import sys
import os
import re
import glob
import fire


def hisat2_summary(summary_dir, summary_file):
    hisat_files = glob.glob(f'{summary_dir}/*.hisat2_summary.txt')
    with open(summary_file, 'w') as summ_inf:
        summ_inf.write('Sample_id\tMapping_rate\n')
        for each_file in hisat_files:
            with open(each_file) as each_file_inf:
                for eachline in each_file_inf:
                    if 'Overall alignment rate' in eachline:
                        mapping_rate = re.search(
                            r'Overall alignment rate: (\S+)%', eachline).groups()[0]
                        mapping_rate = float(mapping_rate) / 100
                        each_file_name = os.path.basename(each_file)
                        sample_id = re.match(
                            '(\S+).hisat2_summary.txt', each_file_name).groups()[0]
                        summ_inf.write(f'{sample_id}\t{mapping_rate:.3f}\n')


if __name__ == '__main__':
    fire.Fire(hisat2_summary)
