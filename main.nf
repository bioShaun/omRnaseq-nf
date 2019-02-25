#!/usr/bin/env nextflow

// TODO LIST
// - different libtype choice
// - optional analysis flag: rmats / snp / ...

def helpMessage() {
    log.info """

    A RNAseq analysis pipeline

    Usage:

    ==================================================================

    References
      --fasta                       Path to Fasta reference
      --gtf                         Path to reference GTF file
      --bed12                       Path to reference bed12 file
      --chr_file                    Chromosome list to split bedfile
      --hisat2_index                Hisat2 index directory
      --go_anno                     GO annotation file
      --gene_length                 Gene length file
      --kegg_anno                   KEGG blasttab annotation
      --kegg_abbr                   KEGG abbr
      --anno_file                   Gene annotation file
      --pathway_db                  KEGG pathway database
    
    Optional ref files
      --transfactor                 Transfactor annotation
      --split_bed                   Split bed file for snp calling

    Mandatory arguments:
      --proj_name                   Project name for output name
      --reads                       Path to input data (must be surrounded with quotes)
      --sample_group                Sample vs Group file

    Other options:
      --lib_type                    RNAseq library type, default is "unstranded"
      --exp_diff_pval               Expression differential analysis p.adjust cutoff
      --exp_lgfc                    Expression differential analysis log2foldchange cutoff
      --contrast                    Differential analysis compare
      --outdir                      Path to store index files, default is directory store fasta
      --snp                         Apply snp analysis
      --snp_quality                 Quality cutoff for snp analysis
      --snp_depth                   Depth cutoff for snp analysis
      --ppi                         PPI link file
      --proten2gene                 PPI protein id <-> gene id map file
      --ref_flat                    RefFlat format annotation
      --kegg_pathway                Run kegg_pathway analysis

    """.stripIndent()
}

// workflow internal path&files
script_dir = file("$baseDir/script/")

 // Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


def check_ref_exist = {file_path, file_type ->
    if (file_path) {
        file_path = file(file_path)
        if( !file_path.exists() ) exit 1, "${file_type} file not found: ${file_path}"
        return file_path
    } else {
        exit 1, "${file_type} is required!"
    }
}
/*
 * SET UP CONFIGURATION VARIABLES
*/
// reference
params.fasta = false
params.hisat2_index = false
params.gtf = false
params.chr_file = false
params.bed12 = false
params.kegg_abbr = false
params.kegg_backgroud = false
params.transfactor = false
params.snp = false
// 30 for real data
params.snp_quality = 30
// 5 for read data
params.snp_depth = 5
params.snpEff = '/public/software/snpEff/snpEffv4.3T/'
params.pathway_db = '/public/database/kegg_html/'
params.snpEff_db = false
params.anno_file = false
params.ref_flat = false
params.kegg_pathway = false
params.proj_name = 'OMS_test'

if (!params.kegg_backgroud) {
    kegg_backgroud = params.kegg_abbr
} else {
    kegg_backgroud = params.kegg_backgroud
}

gtf = check_ref_exist(params.gtf, 'gtf')
bed12 = check_ref_exist(params.bed12, 'bed12')
fasta = check_ref_exist(params.fasta, 'fasta')
fai_idx = file("${fasta}.fai")
fai_idx = check_ref_exist(fai_idx, 'fasta fai')
fasta_path = fasta.getParent() 
go_anno = check_ref_exist(params.go_anno, 'GO')
gene_length = check_ref_exist(params.gene_length, 'Gene Length')
kegg_anno = check_ref_exist(params.kegg_anno, 'KEGG')
anno_file = check_ref_exist(params.anno_file, 'gene annotation')

// hisat index
hisat2_index = Channel
                    .fromPath("${params.hisat2_index}/*.ht2")
                    .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }

alignment_splicesites = file("${params.hisat2_index}/${gtf.baseName}.hisat2_splice_sites.txt")

// 
params.outdir = false
params.reads = false
params.insert_size = 350
params.exp_diff_pval = 0.05
params.exp_lgfc = 1
params.lib_type = ""
params.contrast = false
params.sample_group = false
sample_group = check_ref_exist(params.sample_group, 'Sample vs Group file')


// Prepare analysis fastq files
Channel
    .fromFilePairs("${params.reads}")
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n!" }
    .set { raw_fq_files }

// contrast file
if (params.contrast) {
    contrast_file = file(params.contrast)
} else {
    process generate_contrast {
        publishDir "${params.outdir}/${params.proj_name}/configure" , mode: 'copy'

        input:
        file sample_group from sample_group
        
        output:
        file 'contrast.ini' into contrast_file
        
        script:
        """
        pipenv run python ${script_dir}/utils/make_analysis_compare.py \\
            generate_contrast \\
            --sample-group ${sample_group} \\
            --contrast-file contrast.ini
        """
    }
}


/*
 * Fastp
 */

process fastp {

    tag "${name}"

    module "fastp/0.19.5"

    publishDir "${params.outdir}/${params.proj_name}/data/fastp_trimmed_reads/${name}", mode: 'copy'

    input:
    set name, file(reads) from raw_fq_files

    output:
    file "*trimmed.R*.fq.gz" into trimmed_reads, kallisto_reads
    file "${name}.json" into fastp_json
    file "${name}.html"

    cpus = 4

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${name}.trimmed.R1.fq.gz \\
        --out2 ${name}.trimmed.R2.fq.gz \\
        --json ${name}.json \\
        --html ${name}.html    
    """
}  

/*
* fq quality control
*/
process reads_qc_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/reads_qc/", mode: 'copy'

    input:
    file 'fastp_json/*' from fastp_json.collect()
    
    output:
    file "data.summary.csv"
    file "reads_gc"
    file "reads_quality"
    file "reads_filter"
    
    script:
    """
    pipenv run python ${script_dir}/reads_qc/extract_fastp_info.py \\
        --fastp-dir fastp_json \\
        --outdir .

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/reads_qc/gc_plot.R \\
        --gc_dir reads_gc

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/reads_qc/reads_quality_plot.R \\
        --rq_dir reads_quality

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/reads_qc/filter_stats.R \\
        --filter_dir reads_filter   
    """
}


 /*
 * Hisat2 mapping
 */
process hisat2_align {

    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/data/bam/${name}", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".hisat2_summary.txt") > 0 ? "$filename" : null}

    input:
    file reads from trimmed_reads
    file index from hisat2_index.collect()
    file gtf from gtf
    file alignment_splicesites from alignment_splicesites

    output:
    file "${name}.bam" into hisat2_bam, rnaseq_matrix_bam
    file "${name}.hisat2_summary.txt" into alignment_logs

    cpus = 16

    script:
    name = reads[0].toString() - '.trimmed.R1.fq.gz'
    index_base = index[0].toString() - ~/.\d.ht2/
    """
    hisat2 -x ${index_base} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --known-splicesite-infile ${alignment_splicesites} \\
            --no-mixed \\
            --no-discordant \\
            -p ${task.cpus} \\
            --met-stderr \\
            --new-summary \\
            --summary-file ${name}.hisat2_summary.txt \\
            --rg-id ${name} \\
            --rg SM:${name}\\
            --rg LB:${name} \\
            | samtools view -bS -F 4 -F 8 -F 256 - > ${name}.bam
    """
}

/*
* mapping stats
*/
if (params.ref_flat) {
    ref_flat = check_ref_exist(params.ref_flat, 'Ref flat')
    process RnaSeqMetrics {
        tag "${name}"

        publishDir "${params.outdir}/${params.proj_name}/data/bam/${name}", mode: 'copy'

        input:
        file bam from rnaseq_matrix_bam
        file ref_flat from ref_flat
        
        output:
        file "${name}.RNA_Metrics" into rnaseq_matrix
        
        script:
        name = bam.baseName
        """
        java -jar /public/software/picard/2.17.3/picard.jar CollectRnaSeqMetrics \\
            I=${bam} \\
            O=${name}.RNA_Metrics \\
            REF_FLAT=${ref_flat} \\
            STRAND=NONE
        """
    }
}

/*
* mapping summary
*/
process mapping_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/mapping/", mode: 'copy'

    input:
    file 'rnaseq_matrix/*' from rnaseq_matrix.collect()
    file 'hisat2_summary/*' from alignment_logs.collect()
    
    output:
    file "overall_mapping_stats.txt"
    file "mapping_portion.p*"
    
    script:
    """
    pipenv run python ${script_dir}/mapping/hisat_mapping_rate.py \\
        --summary-dir hisat2_summary \\
        --summary-file mapping_rate.txt

    pipenv run python ${script_dir}/mapping/rnaseq_matrix_summary.py \\
        --analysis_dir rnaseq_matrix \\
        --summary_file mapping_genome_portion.txt

    pipenv run python ${script_dir}/utils/merge_files.py \\
        mapping_rate.txt \\
        mapping_genome_portion.txt \\
        overall_mapping_stats.txt \\
        --by_colname

    Rscript ${script_dir}/mapping/mapping_portion.R \\
        --mapping_file overall_mapping_stats.txt \\
        --out_prefix mapping_portion
    """
}


/*
* Hisat2 bam sort
*/
process hisat2_sortOutput {
    tag "${name}"

    module "samtools/1.9"

    publishDir "${params.outdir}/${params.proj_name}/data/bam/${name}", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".sorted.bam") > 0 ? "$filename" : null}

    input:
    file hisat2_bam

    output:
    file "${name}.sorted.bam" into bam_scallop, bam_for_genebody, bam_saturation, bam_rmats, bam_snp
    file "${name}.fixmate.bam"

    cpus = 8

    script:
    name = hisat2_bam.baseName
    """
    samtools fixmate --threads ${task.cpus} \\
        -m ${hisat2_bam} ${name}.fixmate.bam

    samtools sort \\
        ${name}.fixmate.bam \\
        -m 2400M\\
        --threads ${task.cpus} \\
        -o ${name}.sorted.bam
    """
}


/*
* Scallop assembly
*/
process scallop {

    tag "${name}"

    publishDir "${params.outdir}/${params.proj_name}/data/assembly/${name}", mode: 'copy'

    input:
    file bam from bam_scallop
    
    output:
    file "*gtf" into assembly_gtf

    script:
    name = bam.baseName - '.sorted'
    """
    scallop \\
        -i ${bam} \\
        -o ${name}.gtf
    """
}


/*
* Stringtie merge scallop output
*/
process stringtie_merge {

    publishDir "${params.outdir}/${params.proj_name}", mode: 'copy',
            saveAs: {filename ->
                if (filename == "novel.gtf") "result/assembly/$filename"
                else if (filename == "assembly.gtf") "result/assembly/$filename"
                else if (filename == "transcripts.gtf") "data/ref/$filename"
                else null
            }

    input:
    file "gtf/*" from assembly_gtf.collect()
    file gtf from gtf
    
    output:
    file "transcripts.gtf" into kallisto_idx_gtf, quant_gtf, rmats_gtf
    file "novel.gtf" into novel_gtf, novel_anno_gtf
    file "assembly.gtf"
    file "gene_trans.map" into gene2tr_quant, gene2tr_anno

    script:
    """
    ls gtf/* > gtf.list

    stringtie --merge \\
        -G ${gtf} \\
        -m 200 \\
        -o assembly.gtf \\
        -T 1 \\
        -f 0.1 \\
        gtf.list

    gffcompare \\
        -r ${gtf} \\
        -R assembly.gtf \\
        -o cmp2ref

    pipenv run python ${script_dir}/assembly/novel_gtf_from_gffcompare.py \\
        --compare-gtf cmp2ref.annotated.gtf \\
        --outfile novel.gtf

    cat ${gtf} novel.gtf > transcripts.gtf

    pipenv run python ${script_dir}/assembly/get_gene_to_trans.py \\
        --gff transcripts.gtf \\
        --out_dir .
    """
}

/*
* Make kallisto index for quantification
*/
process mk_kallisto_index {
    tag "KALLISTO index"

    module "kallisto/0.45.0"

    publishDir "${params.outdir}/${params.proj_name}", mode: 'copy',
        saveAs: {filename ->
            if (filename == "novel.fa") "result/assembly/$filename"
            else if (filename == "transcripts.fa.kallisto_idx") "data/ref/$filename"
            else if (filename == "transcripts.fa") "data/ref/$filename"
            else null
        }

    input:
    file gtf from kallisto_idx_gtf
    file novel_gtf from novel_gtf
    file fasta from fasta
    
    output:
    file "transcripts.fa" into merged_fa
    file "novel.fa" into novel_fa
    file "transcripts.fa.kallisto_idx" into kallisto_idx
    
    script:
    """
    gffread ${gtf} \\
        -g ${fasta} \\
        -w transcripts.fa

    gffread ${novel_gtf} \\
        -g ${fasta} \\
        -w novel.fa    

    kallisto index \\
        -i transcripts.fa.kallisto_idx \\
        transcripts.fa
    """
}

/*
* novel gene annotation
*/
process novel_gene_annotation {

    publishDir "${params.outdir}/${params.proj_name}/result/assembly/", mode: 'copy'

    input:
    file novel_fa from novel_fa
    file novel_gtf from novel_anno_gtf
    file gene2tr from gene2tr_anno
    
    output:
    file "novel.gene.annotation.csv"
    
    script:
    """
    pipenv run python ${script_dir}/novel_gene/gene_location.py \\
        --gtf-file ${novel_gtf} \\
        --outfile novel.loc.txt

    pipenv run python ${script_dir}/novel_gene/annotate_novel_gene.py \\
        --input-file ${novel_fa} \\
        --gene2tr ${gene2tr}
    
    pipenv run python ${script_dir}/utils/merge_files.py \\
        novel.loc.txt \\
        novel.annotation.detail.txt \\
        novel.gene.annotation.csv \\
        --na_rep '--' \\
        --method left \\
        --out_sep ','
    """
}


/*
* Kallisto quant
*/
process kallisto {

    tag "${name}"

    module "kallisto/0.45.0"

    publishDir "${params.outdir}/${params.proj_name}/data/kallisto/", mode: 'copy'

    input:
    file reads from kallisto_reads
    file index from kallisto_idx
    
    output:
    file name into kallisto_out

    cpus 4
    
    script:
    name = reads[0].toString() - '.trimmed.R1.fq.gz'
    """
    kallisto quant \
        --threads ${task.cpus} \
        -i ${index} \
        --output-dir=${name} \
        ${reads}   
    """
}

// for differential analysis
/*
* Load Kallisto results
*/

process rmats_cfg {

    publishDir "${params.outdir}/${params.proj_name}/configure" , mode: 'copy'

    input:
    file sample_group from sample_group
    file contrast_file from contrast_file
    file 'bam/*' from bam_rmats.collect()
    
    output:
    file 'rmats_compare/*' into rmats_compare
    
    script:
    """
    pipenv run python ${script_dir}/utils/make_analysis_compare.py \\
        rmats-sample-files \\
        --bam-dir ./bam \\
        --sample-group ${sample_group} \\
        --out-dir rmats_compare \\
        --contrast ${contrast_file}
    """
}

rmats_compare
    .flatMap()
    .into { flat_rmats_compare; diff_compare }

process load_kallisto_results {

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/", mode: 'copy',
        saveAs: {filename -> filename == "deg_input.RData" ? null : "$filename"}

    input:
    file 'kallisto/*' from kallisto_out.collect()
    file sample_group from sample_group
    file gene2tr from gene2tr_anno
    
    output:
    file 'expression_summary' into expression_summary
    file 'deg_input.RData' into deg_obj
    
    script:
    """   
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/kallisto_to_table.R \\
        --kallisto_dir kallisto \\
        --sample_inf ${sample_group} \\
        --gene2tr ${gene2tr} \\
        --out_dir expression_summary
    """
}

/*
* Differential analysis
*/
process diff_analysis {
    tag "${compare}"

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/differential_analysis/", mode: 'copy'

    input:
    file compare from diff_compare
    file deg_obj from deg_obj
    file sample_group from sample_group
    
    output:
    file compare into diff_out_go, diff_out_kegg, diff_out_all, diff_tf, diff_out_ppi, diff_out_enrich_plot, pathway_compare, diff_summary
    
    script:
    """
    mv ${compare} ${compare}_compare

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/diff_edgeR.R \\
        --deg_rdata ${deg_obj} \\
        --compare ${compare} \\
        --sample_inf ${sample_group} \\
        --out_dir ${compare} \\
        --qvalue ${params.exp_diff_pval} \\
        --logfc ${params.exp_lgfc}
    """
}


/*
* diff summary
*/
process diff_exp_summary {

    publishDir "${params.outdir}/${params.proj_name}/result/quantification/expression_summary", mode: 'copy'

    input:
    file expression_summary from expression_summary
    file "diff/*" from diff_summary.collect()
    file sample_group from sample_group

    output:
    file "Diff.gene*"
    file "cluster_data"
    
    script:
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/quant/quant_report.R \\
        --exp_dir ${expression_summary} \\
        --diff_dir diff \\
        --sample_inf ${sample_group} \\
        --out_dir .
    """
}

/*
* enrichment analysis
*/
// go enrichment

def compare2reg = {compare ->
    compare_name = compare.baseName
    (up_name, down_name) = "${compare.baseName}".split('_vs_')
    return ["${compare_name};ALL", "${compare_name};${up_name}-UP", "${compare_name};${down_name}-UP"]
}

diff_out_all
    .flatMap(compare2reg)
    .into { go_compare; kegg_compare; ppi_compare }


process go_analysis {
    tag "${go_compare}"

    publishDir "${params.outdir}/${params.proj_name}/result/enrichment/go/", mode: 'copy'

    input:
    file "diff/*" from diff_out_go.collect()
    val go_compare from go_compare
    file go_anno from go_anno
    file gene_length from gene_length
    
    output:
    file "${compare}/${compare}.${reg}.go.enrichment.txt" into go_out
    file "${compare}/DAG/*" into dag_plot
    file "${compare}/${compare}.${reg}.go.enrichment.barplot*"
    
    script:
    (compare, reg) = go_compare.split(';')
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/enrichment/go_analysis.R \\
        --name ${compare}.${reg} \\
        --gene_list diff/${compare}/${compare}.${reg}.edgeR.DE_results.diffgenes.txt \\
        --go_anno ${go_anno} \\
        --gene_length ${gene_length} \\
        --out_dir ${compare}

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/enrichment/enrich_bar.R \\
        --enrich_file ${compare}/${compare}.${reg}.go.enrichment.txt
    """
}


//KEGG enrichment
process kegg_analysis {

    tag "${compare_reg}"

    publishDir "${params.outdir}/${params.proj_name}/result/enrichment/kegg/", mode: 'copy',
        saveAs: {filename -> filename.indexOf("kegg.enrichment") > 0 ? "$filename" : null}   

    input:
    val compare_reg from kegg_compare
    file "diff/*" from diff_out_kegg.collect()
    file kegg_anno from kegg_anno
    
    output:
    file "${compare}/${compare}.${reg}.kegg.enrichment.txt" into kegg_out, kegg_pathway_input
    file "${compare}.${reg}.blasttab" into kegg_pathway_blast
    file "${compare}/${compare}.${reg}.kegg.enrichment.barplot*"
    
    script:
    (compare, reg) = compare_reg.split(';')
    """
    #!/bin/bash

    pipenv run python ${script_dir}/utils/extract_info_by_id.py \\
        --id diff/${compare}/${compare}.${reg}.edgeR.DE_results.diffgenes.txt \\
        --table ${kegg_anno} \\
        --output ${compare}.${reg}.blasttab 

    mkdir -p ${compare}
    
    source /usr/bin/virtualenvwrapper.sh
    workon oms_pub

    run_kobas.py \\
        -i ${compare}.${reg}.blasttab \\
        -t blastout:tab \\
        -s ${params.kegg_abbr} \\
        -b ${kegg_backgroud} \\
        -d K \\
        -o ${compare}/${compare}.${reg}.kegg.enrichment.txt

    python ${script_dir}/enrichment/treat_kegg_table.py \\
        --kegg ${compare}/${compare}.${reg}.kegg.enrichment.txt

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/enrichment/enrich_bar.R \\
        --enrich_file ${compare}/${compare}.${reg}.kegg.enrichment.txt
    """
}


/*
* KEGG pathway
*/
if (params.kegg_pathway) {
    process kegg_pathway {
        tag "${compare_reg}"

        publishDir "${params.outdir}/${params.proj_name}/result/enrichment/kegg/${compare}", mode: 'copy'

        input:
        file kegg_out from kegg_pathway_input
        file blast_out from kegg_anno
        file "diff/*" from pathway_compare.collect()
        
        output:
        file "${compare_reg}.pathway"
        
        script:
        compare_reg = kegg_out.baseName - '.kegg.enrichment'
        compare = kegg_out.baseName - ~/.(\w+-UP|ALL).kegg.enrichment$/
        """
        mkdir ${compare_reg}.pathway
        pipenv run python ${script_dir}/enrichment/pathway_html.py \\
            --pathway-db ${params.pathway_db} \\
            --blasttab ${blast_out} \\
            --kegg-table ${kegg_out} \\
            --diff-table diff/${compare}/${compare}.edgeR.DE_results.txt \\
            --outdir ${compare_reg}.pathway
        """
    }
}

/*
* Diff gene ppi analysis
*/
if (params.ppi && params.proten2gene) {
    ppi_anno = check_ref_exist(params.ppi, 'PPI')
    proten2gene = check_ref_exist(params.proten2gene, 'proten2gene id map')
    process ppi_analysis {

        tag "${compare}"

        publishDir "${params.outdir}/${params.proj_name}/result/ppi/", mode: 'copy'

        input:
        file compare from diff_out_ppi
        file ppi_anno from ppi_anno
        file proten2gene from proten2gene
        file anno_file from anno_file
        
        output:
        file "${compare}.ppi.csv"
        
        script:
        """
        pipenv run python ${script_dir}/ppi/diff_reg_list.py \\
            --diff-dir ${compare} \\
            --outfile ${compare}.deff_reg.txt 

        pipenv run python ${script_dir}/ppi/diff_gene_ppi.py \\
            --ppi-link ${ppi_anno} \\
            --proten2gene ${proten2gene} \\
            --diff-gene-list ${compare}.deff_reg.txt \\
            --outfile ${compare}.ppi.csv \\
            --anno-file ${anno_file}
        """
    }
}

/*
* transfactor in diff genes
*/
if (params.transfactor) {
    process transfactor_anno {

        publishDir "${params.outdir}/${params.proj_name}/result/transfactor/", mode: 'copy'

        input:
        file 'diff/*' from diff_tf.collect()
        file tf_file from file(params.transfactor)
        
        output:
        file "Diff_TF_number.stats.csv"
        file "TF_in_diff_compare.stats.csv"
        
        script:
        """
        pipenv run python ${script_dir}/transfactor/diff_tf.py \\
            --tf-file ${tf_file} \\
            --diff-dir ./diff \\
            --outdir .
        """
    }
}



/*
* Rseqc create BigWig coverage
*/
process createBigWig {
    tag "${name}"

    input:
    file bam from bam_for_genebody

    output:
    file "*.bigwig" into bigwig_for_genebody

    cpus 8

    script:
    name = bam.baseName - '.sorted'
    """
    samtools index ${bam}
    pipenv run bamCoverage -b ${bam} -p ${task.cpus} -o ${name}.bigwig
    """
}

/*
 * Rseqc genebody coverage
 */
process genebody_coverage {
    tag "${name}"

    input:
    file bigwig from bigwig_for_genebody
    file bed12 from bed12

    output:
    file "${name}.geneBodyCoverage.txt" into genebody_coverage_results

    cpus 8

    script:
    name = bigwig.baseName
    """
    pipenv run geneBody_coverage2.py \\
        -i ${bigwig} \\
        -o ${name} \\
        -r ${bed12}
    """
}


/*
 * Rseqc genebody coverage plot
 */
process genebody_coverage_plot {

    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc" , mode: 'copy'

    input:
    file "genebody_coverage/*" from genebody_coverage_results.collect()

    output:
    file "genebody_coverage"


    script:
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/rnaseq_qc/genebody_cov.R \\
        --cov_dir genebody_coverage
    """
}

/*
 * Rseqc gene exp saturation analysis
 */
process saturation {
    tag "${name}"

    input:
    file bam from bam_saturation
    file bed12 from bed12

    output:
    file "${name}.rawCount.xls" into count_saturation

    cpus 12

    script:
    name = bam.baseName - '.sorted'
    """
    pipenv run RPKM_saturation.py \\
        -i ${bam} \\
        -o ${name} \\
        -r ${bed12}
    """
}

/*
 * Rseqc gene exp saturation analysis
 */
process saturation_plot {

    publishDir "${params.outdir}/${params.proj_name}/result/rnaseq_qc/" , mode: 'copy'

    input:
    file 'sequencing_saturation/*' from count_saturation.collect()

    output:
    file "sequencing_saturation"

    script:
    """
    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/rnaseq_qc/seq_saturation.R \\
        --saturation_count sequencing_saturation
    """
}

/*
* rMATS for alternative analysis
*/
process rmats {
    tag "${each_compare}"

    publishDir "${params.outdir}/${params.proj_name}/result/rmats" , mode: 'copy'

    input:
    file each_compare from flat_rmats_compare
    file gtf from rmats_gtf

    output:
    file "${each_compare}/*MATS.JC.txt" into rmats_jc_out
    file "${each_compare}/*MATS.JCEC.txt" into rmats_jcec_out
    file "${each_compare}/*as_summary.p*" into rmats_plot

    cpus 12

    script:
    """
    mv ${each_compare} ${each_compare}_compare

    python /public/software/rMATS/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py \\
        --b1 ${each_compare}_compare/b1.txt \\
        --b2 ${each_compare}_compare/b2.txt \\
        -t paired \\
        --readLength 150 \\
        --gtf ${gtf} \\
        --od ${each_compare} \\
        --libType fr-unstranded

    /public/software/R/R-3.5.1/executable/bin/Rscript ${script_dir}/rmats/rmats_plot.R \\
        --as_dir ${each_compare}
    """
}


/*
* snp
*/

if (params.snp) {
    // check fasta index
    fa_dict = file("${fasta_path}/${fasta.baseName}.dict")
    fa_dict = check_ref_exist(fa_dict, 'fasta dict')

    // prepare split bed files
    Channel
        .fromPath("${params.split_bed}/genome/*bed")
        .ifEmpty { exit 1, "Cannot find any bed file in directory: ${params.split_bed}/genome\n!" }
        .into { haplotype_beds; combine_gvcf_beds }

    /*
    * Reads mark duplication and recalibration
    */
    process bam_remove_duplicate {
        tag "${sample_name}"

        publishDir "${params.outdir}/${params.proj_name}/data/bam/${sample_name}", mode: 'copy'

        input:
        file bam from bam_snp
    
        output:
        file "${sample_name}.rmdup.bam" into br_rmdup_bam
    
        cpus = 8

        script:
        sample_name = bam.baseName - '.sorted'
        """
        samtools markdup -r -s \\
            --threads ${task.cpus} \\
            ${bam} \\
            ${sample_name}.rmdup.bam
        """
    }

    process SplitNCigarReads {

        tag "${sample_name}"
            
        input:
        file bam from br_rmdup_bam
        file fasta from fasta
        file fa_dict from fa_dict
        file fai_idx from fai_idx
        
        output:
        file "${sample_name}.SplitN.bam" into splitn_bam

        cpus = 8
        
        script:
        sample_name = bam.baseName - '.rmdup'
        """
        gatk SplitNCigarReads \\
            --reference ${fasta} \\
            --input ${bam} \\
            --output ${sample_name}.SplitN.bam \\
            --max-reads-in-memory 1000000  \\
        """
    }

    /*
    *  GATK HaplotypeCaller
    */
    process gatk_HaplotypeCaller {
        tag "${sample_name}|${chr_name}"

        input:
        file bam from splitn_bam
        each file(bed) from haplotype_beds
        file refer from fasta
        file refer_fai from fai_idx
        file refer_dict from fa_dict    

        output:
        file "${sample_name}.${chr_name}.hc.g.vcf.gz" into sample_gvcf
        file "${sample_name}.${chr_name}.hc.g.vcf.gz.tbi" into sample_gvcf_index
        
        cpus = 8

        script:
        sample_name = bam.baseName - '.SplitN'
        chr_name = bed.baseName
        """
        gatk HaplotypeCaller  \\
            --input ${bam} \\
            --output ${sample_name}.${chr_name}.hc.g.vcf.gz \\
            --reference ${refer} \\
            --intervals ${bed} \\
            --emit-ref-confidence GVCF \\
            --read-filter AmbiguousBaseReadFilter \\
            --read-filter MappingQualityReadFilter \\
            --read-filter NonZeroReferenceLengthAlignmentReadFilter \\
            --read-filter ProperlyPairedReadFilter \\
            --minimum-mapping-quality 30         
        """
    }

    /*
    * GATK CombineGVCFs
    */
    process gatk_CombineGVCFs {
        tag "Chrom: ${chr_name}"

        publishDir "${params.outdir}/${params.proj_name}/data/snp/gvcf/", mode: 'copy'  

        input:
        file ('gvcf/*') from sample_gvcf.collect()
        file ('gvcf/*') from sample_gvcf_index.collect()
        each file(bed) from combine_gvcf_beds
        file refer from fasta
        file refer_fai from fai_idx
        file refer_dict from fa_dict    
        
        output:
        file "all_sample.${chr_name}.g.vcf.gz" into merged_sample_gvcf
        file "all_sample.${chr_name}.g.vcf.gz.tbi" into merged_sample_gvcf_index
        
        script:
        chr_name = bed.baseName
        """
        ls gvcf/*.${chr_name}.hc.g.vcf.gz > ${chr_name}.gvcf.list

        gatk CombineGVCFs \\
            --output all_sample.${chr_name}.g.vcf.gz \\
            --reference ${refer} \\
            --variant ${chr_name}.gvcf.list
        """
    }

    /*
    * GATK GenotypeGVCFs
    */
    process gatk_GenotypeGVCFs {
        tag "Chrom: ${chr_name}"

        publishDir "${params.outdir}/${params.proj_name}/data/snp/vcf/", mode: 'copy'

        input:
        file gvcf from merged_sample_gvcf
        file gvcf_index from merged_sample_gvcf_index
        file refer from fasta
        file refer_fai from fai_idx
        file refer_dict from fa_dict    
        
        output:
        file "${vcf_prefix}.vcf.gz" into merged_sample_vcf
        file "${vcf_prefix}.vcf.gz.tbi" into merged_sample_vcf_index
        
        script:
        vcf_prefix = gvcf.baseName - '.g.vcf'
        chr_name = vcf_prefix - 'all_sample.'
        """
        gatk GenotypeGVCFs \\
            --reference ${refer} \\
            --variant ${gvcf} \\
            --output ${vcf_prefix}.vcf.gz
        """
    }

    /*
    * Concat vcf
    */
    process concat_vcf {

        publishDir "${params.outdir}/${params.proj_name}/result/snp/", mode: 'copy'   

        input:
        file ('vcf/*') from merged_sample_vcf.collect()
        file ('vcf/*') from merged_sample_vcf_index.collect()

        output:
        file "all_sample.raw.vcf.gz" into all_sample_raw_vcf
        file "all_sample.raw.vcf.gz.tbi" into all_sample_raw_vcf_idx
        
        script:
        """
        bcftools concat \\
            vcf/*.vcf.gz | \\
            bgzip > all_sample.raw.vcf.gz

        tabix -p vcf all_sample.raw.vcf.gz
        """
    }

    /*
    * Basic quality filter
    */
    process vcf_base_qual_filter {

        publishDir "${params.outdir}/${params.proj_name}/result/snp/", mode: 'copy'   

        input:
        file raw_vcf from all_sample_raw_vcf
        file raw_vcf_idx from all_sample_raw_vcf_idx
        
        output:
        file "all_sample.hq.vcf.gz" into all_hq_vcf
        file "all_sample.hq.vcf.gz.tbi" into all_hq_vcf_idx
        
        script:
        """
        bcftools filter -s LowQual -e '%QUAL<${params.snp_quality} || INFO/DP<${params.snp_depth}' \\
            all_sample.raw.vcf.gz | \\
            grep -v LowQual | bgzip > all_sample.hq.vcf.gz

        tabix -p vcf all_sample.hq.vcf.gz
        """
    }


    /*
    * snpeff for combined vcf
    */
    process snpEff_for_all {

        publishDir "${params.outdir}/${params.proj_name}/result/snp/", mode: 'copy'

        when:
        params.snpEff_db

        input:
        file vcf from all_hq_vcf
        file vcf_idx from all_hq_vcf_idx
        
        output:
        file "all_sample.hq.vcf.stat.csv"
        file "all_sample.hq.vcf.stat.html" 
        file "all_sample.ann.vcf.gz" into all_sample_anno_vcf
        file "all_sample.ann.vcf.gz.tbi" into all_sample_anno_vcf_idx
        
        script:
        """
        java -Xmx10g -jar ${params.snpEff}/snpEff.jar \\
            -c ${params.snpEff}/snpEff.config \\
            -csvStats all_sample.hq.vcf.stat.csv \\
            -htmlStats all_sample.hq.vcf.stat.html \\
            -v ${params.snpEff_db}  \\
            ${vcf} \\
            | bgzip > all_sample.ann.vcf.gz
        
        tabix -p vcf all_sample.ann.vcf.gz
        """
    }

    /*
    * SNP Table
    */
    process snp_table {

        publishDir "${params.outdir}/${params.proj_name}/result/snp/", mode: 'copy'

        input:
        file vcf from all_sample_anno_vcf
        file vcf_idx from all_sample_anno_vcf_idx
        file refer from fasta
        file refer_fai from fai_idx
        file refer_dict from fa_dict  

        output:
        file "all_sample.vcf.table.txt" into om_vcf_table

        script:
        """
        gunzip -c ${vcf} > ${vcf.baseName}

        python /public/scripts/Reseq/omtools/extractTableFromsnpEff.py \\
            -v ${vcf.baseName} \\
            -o all_sample.vcf.table.txt
        """
    }

}

