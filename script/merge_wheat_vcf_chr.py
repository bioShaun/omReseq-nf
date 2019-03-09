import fire
import pandas as pd
from pathlib import PurePath
from Bio import bgzf


def merge_vcf_split_chr(vcf_file, split_chr_inf):
    split_bed_df = pd.read_table(split_chr_inf, index_col=3,
                                 header=None,
                                 names=['chrom', 'start', 'end'])
    vcf_file = PurePath(vcf_file)
    # a.xxx.[vcf|table].[gz|txt] -> a.xxx.catChr.[vcf|table].[gz|txt]
    merge_chr_vcf_name = vcf_file.name
    name_prefix = '.'.join(merge_chr_vcf_name.split('.')[:-2])
    name_suffix = '.'.join(merge_chr_vcf_name.split('.')[-2:])
    merge_chr_vcf_name = '{0}.catChr.{1}'.format(name_prefix, name_suffix)
    merge_chr_vcf_file = str(vcf_file.with_name(merge_chr_vcf_name))

    is_gz_file = vcf_file.suffix == '.gz'

    if is_gz_file:
        cat_vcf_inf = bgzf.BgzfWriter(merge_chr_vcf_file)
        split_vcf_inf = bgzf.BgzfReader(vcf_file)
    else:
        cat_vcf_inf = open(merge_chr_vcf_file, 'w')
        split_vcf_inf = open(vcf_file)

    for eachline in split_vcf_inf:
        eachline_inf = eachline.split('\t')
        chrom = eachline_inf[0]
        if chrom in split_bed_df.index:
            new_chrom, offset = split_bed_df.loc[chrom, [
                'chrom', 'start']]
            eachline_inf[0] = new_chrom
            eachline_inf[1] = str(offset + int(eachline_inf[1]))
            eachline = '\t'.join(eachline_inf)
        cat_vcf_inf.write(eachline)

    cat_vcf_inf.close()
    split_vcf_inf.close()


if __name__ == '__main__':
    fire.Fire(merge_vcf_split_chr)
