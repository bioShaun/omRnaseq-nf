import click
import pandas as pd
import os
import glob
import re


def get_rna_matrix_inf(m_file):
    with open(m_file) as file_inf:
        flag = 0
        for eachline in file_inf:
            if 'PF_BASES' in eachline:
                flag = 1
                continue
            if flag:
                eachline_inf = eachline.strip().split('\t')
                return eachline_inf


OUT_HEADER = ['Exon',
              'Intron',
              'Intergenic']


@click.command()
@click.option('-d', '--analysis_dir', type=click.Path(exists=True),
              help='rnaseq matrix output directory.', required=True)
@click.option('-f', '--summary_file',
              help='rnaseq matrix summary output file.', required=True)
def main(analysis_dir, summary_file):
    summary_file_inf = open(summary_file, 'w')
    summary_file_inf.write('Sample_id\t{h}\n'.format(h='\t'.join(OUT_HEADER)))
    rnaseq_matrix_files = glob.glob(f'{analysis_dir}/*.RNA_Metrics')
    for each_file in rnaseq_matrix_files:
        each_file_name = os.path.basename(each_file)
        sample_id = re.match('(\S+).RNA_Metrics', each_file_name).groups()[0]
        each_m_inf = get_rna_matrix_inf(each_file)
        summary_file_inf.write('{0}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n'.format(
            sample_id, float(each_m_inf[20]),
            float(each_m_inf[18]),
            float(each_m_inf[19])
        ))
    summary_file_inf.close()


if __name__ == '__main__':
    main()
