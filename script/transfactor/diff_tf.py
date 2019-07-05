import pandas as pd
import os
import fire
from pathlib import Path


REG_LABEL = ('UP', 'DOWN')


def tf_diff_stats(tf_file, diff_dir, outdir):
    compares = os.listdir(diff_dir)
    outdir = Path(outdir)
    diff_dir = Path(diff_dir)
    diff_dfs = list()
    for each_comp in compares:
        comp_groups = each_comp.split('_vs_')
        comp_diff_dfs = []
        for n, each_reg in enumerate(comp_groups):
            name = f'{each_comp}.{each_reg}-UP.edgeR.DE_results.diffgenes.txt'
            diff_file = diff_dir / each_comp / name
            diff_df = pd.read_table(diff_file, header=None,
                                    index_col=0)
            diff_df.loc[:, each_comp] = REG_LABEL[n]
            comp_diff_dfs.append(diff_df)
        comp_diff_df = pd.concat(comp_diff_dfs)
        diff_dfs.append(comp_diff_df)
    diff_df = pd.concat(diff_dfs, axis=1)
    diff_df_cols = diff_df.columns
    tf_df = pd.read_table(tf_file, index_col=2)
    tf_diff_df = diff_df.reindex(tf_df.index)
    # tf number stats
    melt_tf_diff_df = tf_diff_df.melt().dropna()
    diff_tf_summary = melt_tf_diff_df.groupby(
        ['variable', 'value']).size().unstack(1)
    diff_tf_summary.columns.name = ''
    diff_tf_summary.index.name = 'Group_id'
    diff_tf_summary_file = outdir / 'Diff_TF_number.stats.csv'
    diff_tf_summary.to_csv(diff_tf_summary_file)
    # tf up/down regulation stats
    merged_tf_df = tf_df.merge(diff_df, left_index=True,
                               right_index=True,
                               how='left')
    merged_tf_df = merged_tf_df.sort_values(list(diff_df_cols))
    merged_tf_file = outdir / 'TF_in_diff_compare.stats.csv'
    merged_tf_df.to_csv(merged_tf_file, na_rep='None')


if __name__ == '__main__':
    fire.Fire(tf_diff_stats)
