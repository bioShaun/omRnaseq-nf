import pandas as pd
import fire


def diff_gene_ppi(ppi_link, proten2gene, diff_gene_list, outfile,
                  anno_file=None):
    lk_df = pd.read_table(ppi_link, compression='gzip', delim_whitespace=True)
    lk_df['ordered-cols'] = lk_df.apply(lambda x: '-'.join(
        sorted([x['protein1'], x['protein2']])), axis=1)
    lk_df.drop_duplicates(['ordered-cols'], inplace=True)
    lk_df.drop(['ordered-cols'], axis=1, inplace=True)
    old_cols = lk_df.columns
    pep_gene_map_df = pd.read_table(
        proten2gene, header=None, index_col=0, names=['gene_id'])
    lk_df.loc[:, 'protein1'] = [each.split('.')[-1] for each in lk_df.protein1]
    lk_df.loc[:, 'protein2'] = [each.split('.')[-1] for each in lk_df.protein2]
    lk_df.loc[:, 'gene1'] = [
        pep_gene_map_df.loc[each].gene_id
        if each in pep_gene_map_df.index else None
        for each in lk_df.protein1]
    lk_df.loc[:, 'gene2'] = [
        pep_gene_map_df.loc[each].gene_id
        if each in pep_gene_map_df.index else None
        for each in lk_df.protein2]
    new_cols = old_cols.insert(0, 'gene2')
    new_cols = new_cols.insert(0, 'gene1')
    lk_df.dropna(inplace=True)
    diff_gene_df = pd.read_table(diff_gene_list, header=None,
                                 index_col=0, names=['regulation'])
    diff_gene1_df = diff_gene_df.copy()
    diff_gene1_df.columns = ['gene1_regulation']
    diff_lk_df = lk_df.merge(diff_gene1_df, left_on='gene1',
                             right_index=True)
    diff_gene2_df = diff_gene_df.copy()
    diff_gene2_df.columns = ['gene2_regulation']
    diff_lk_df = diff_lk_df.merge(diff_gene2_df, left_on='gene2',
                                  right_index=True)
    new_cols = new_cols.append(
        pd.Index(['gene1_regulation', 'gene2_regulation']))
    if anno_file is not None:
        anno_df = pd.read_table(anno_file, index_col=0)
        anno_cols = anno_df.columns
        col1 = pd.Index([f'gene1_{each}' for each in anno_cols])
        col2 = pd.Index([f'gene2_{each}' for each in anno_cols])
        anno_df.columns = col1
        diff_lk_df = diff_lk_df.merge(anno_df, left_on='gene1',
                                      right_index=True,
                                      how='left')
        anno_df.columns = col2
        diff_lk_df = diff_lk_df.merge(anno_df, left_on='gene2',
                                      right_index=True,
                                      how='left')
        diff_lk_df.fillna('--', inplace=True)
        new_cols = new_cols.append(col1)
        new_cols = new_cols.append(col2)

    diff_lk_df.to_csv(outfile, index=False,
                      columns=new_cols)


if __name__ == '__main__':
    fire.Fire(diff_gene_ppi)
