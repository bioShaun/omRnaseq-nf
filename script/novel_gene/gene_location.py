import pandas as pd
import gtfparse
import fire


def gene_loc_inf(gtf_file, outfile):
    gtf_df = gtfparse.read_gtf(gtf_file)
    gene_chrom = gtf_df.groupby(['gene_id'])['seqname'].first()
    gene_chrom.name = 'chrom'
    gene_start = gtf_df.groupby(['gene_id'])['start'].min()
    gene_end = gtf_df.groupby(['gene_id'])['end'].max()
    gene_strand = gtf_df.groupby(['gene_id'])['strand'].first()
    gene_loc_df = pd.concat([gene_chrom, gene_start,
                             gene_end, gene_strand], axis=1)
    gene_loc_df.to_csv(outfile, sep='\t')


if __name__ == '__main__':
    fire.Fire(gene_loc_inf)
