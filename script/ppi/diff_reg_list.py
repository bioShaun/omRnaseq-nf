import pandas as pd
import fire
from pathlib import PurePath

REG_LABEL = ('UP', 'DOWN')


def diff_reg_table(diff_dir, outfile):
    diff_dir = PurePath(diff_dir)
    compare_name = diff_dir.name

    def read_reg_file(reg, label):
        nonlocal diff_dir
        nonlocal compare_name
        up_reg_list = diff_dir / f'{compare_name}.{reg}.edgeR.DE_results.diffgenes.txt'
        up_reg_df = pd.read_table(up_reg_list, header=None,
                                  index_col=0)
        up_reg_df.loc[:, 'regulation'] = label
        return up_reg_df

    group_names = compare_name.split('_vs_')
    diff_dfs = list()
    for n, group in enumerate(group_names):
        reg = f'{group}-UP'
        label = REG_LABEL[n]
        each_diff_df = read_reg_file(reg, label)
        diff_dfs.append(each_diff_df)
    diff_df = pd.concat(diff_dfs)
    diff_df.to_csv(outfile, sep='\t', header=False)


if __name__ == '__main__':
    fire.Fire(diff_reg_table)
