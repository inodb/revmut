"""
Find mutations from those found in exome sequenced samples using samtools mpileup and bcftools
"""
import argparse
import pandas as pd
import sys

import revmut

from sufam.mpileup_parser import run_and_get_mutations
from sufam.mutation import get_mutation


def filter_out_mutations_in_normal_deprecated(tumordf, normaldf):
    df = tumordf.merge(normaldf, on=["chrom","pos"], suffixes=("_T", "_N"))

    # filters
    common_al = (df.most_common_al_count_T == df.most_common_count_T) & (df.most_common_al_T == df.most_common_al_N)
    common_indel = (df.most_common_indel_count_T == df.most_common_count_T) & (df.most_common_indel_T == df.most_common_indel_N)
    normal_criteria = ((df.most_common_count_N >= 20) & df.most_common_maf_N > 0.2) | ((df.most_common_count_N < 20) & df.most_common_count_N > 1)
    df = df[~(common_al | common_indel) & normal_criteria]

    # restore column names of tumor
    for c in df.columns:
        if c.endswith("_N"):
            del df[c]
    df.columns = [c[:-2] if c.endswith("_T") else c for c in df.columns]

    return df


def select_only_revertant_mutations_deprecated(bpdf, pos, snv=None, ins=None, dlt=None):
    """
    Selects only mutations that revert the given mutations in a single event.
    """
    if sum([bool(snv), bool(ins), bool(dlt)]) != 1:
        raise(Exception("Should be either snv, ins or del".format(snv)))

    if bool(snv):
        if snv not in ["A","C","G","T"]:
            raise(Exception("snv {} should be A, C, G or T".format(snv)))
        rv = bpdf[(bpdf.pos.astype(float) == pos) & (bpdf.most_common_al == snv) & (bpdf.most_common_al_count == bpdf.most_common_count)]
    elif bool(ins):
        rv = bpdf[((bpdf.most_common_indel.apply(lambda x: len(x) + len(ins) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "+") & (bpdf.most_common_count == bpdf.most_common_indel_count)) |
                ((bpdf.most_common_indel.apply(lambda x: len(ins) - len(x) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "-") & (bpdf.most_common_count == bpdf.most_common_indel_count))]
    elif bool(dlt):
        rv = bpdf[((bpdf.most_common_indel.apply(lambda x: len(x) - len(dlt) % 3 if x else None) == 0) & (bpdf.most_common_indel_type == "+") & (bpdf.most_common_count == bpdf.most_common_indel_count)) |
                ((bpdf.most_common_indel.apply(lambda x: -len(dlt) - len(x) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "-") & (bpdf.most_common_count == bpdf.most_common_indel_count))]
    else:
        # should never happen
        raise(Exception("No mutation given?"))

    # find mutations that delete the mutation partially or completely
    del_mut_df = bpdf[(bpdf.most_common_indel_type == "-") & (bpdf.pos.astype(float) <= pos) & (bpdf.pos.astype(float) + bpdf.most_common_indel.apply(lambda x: len(x) if x else None) >= pos)]
    return pd.concat([rv, del_mut_df], axis=0)


def select_only_revertant_mutations(to_be_reverted_mut, mutations_at_pos):
    muts = []

    dpos = to_be_reverted_mut.pos - mutations_at_pos.pos

    if to_be_reverted_mut.type == ".":
        for k in mutations_at_pos.deletions:
            if dpos >= 0 and len(k) - dpos >= 0:
                muts += [mutations_at_pos.deletions[k]]
    elif to_be_reverted_mut.type == "-":
        for k in mutations_at_pos.deletions:
            if dpos >= 0 and len(k) - dpos >= 0 or \
                    ((len(k) + len(to_be_reverted_mut.change)) % 3 == 0):
                muts += [mutations_at_pos.deletions[k]]
        for k in mutations_at_pos.insertions:
            if (len(k) - len(to_be_reverted_mut.change)) % 3 == 0:
                muts += [mutations_at_pos.insertions[k]]
    elif to_be_reverted_mut.type == "+":
        for k in mutations_at_pos.deletions:
            if dpos >= 0 and len(k) - dpos >= 0 or \
                    ((len(k) - len(to_be_reverted_mut.change)) % 3 == 0):
                muts += [mutations_at_pos.deletions[k]]
        for k in mutations_at_pos.insertions:
            if (len(k) + len(to_be_reverted_mut.change)) % 3 == 0:
                muts += [mutations_at_pos.insertions[k]]
    else:
        raise(Exception("Unkown mutation type"))

    return muts


def find_revertant_mutations(reffa, mutations_tsv, search_bam, normal_bam,
                             outfile):
    """Find revertant mutations. Expecting CHROM, POS, REF, ALT,
    SEARCH_START, SEARCH_END in mutations_tsv"""
    mutsdf = pd.read_csv(mutations_tsv, sep="\t", dtype={"CHROM": str})
    for i, row in mutsdf.iterrows():
        m = get_mutation(row.CHROM, int(row.POS), row.REF, row.ALT)
        search_start = int(row.SEARCH_START)
        search_end = int(row.SEARCH_END)

        muts = run_and_get_mutations(search_bam, m.chrom, search_start,
                                          search_end, reffa)
        nmuts = run_and_get_mutations(normal_bam, m.chrom, search_start,
                                          search_end, reffa)
        nmuts = {muts_at_pos.pos: muts_at_pos for muts_at_pos in nmuts}

        filt_muts = []
        for muts_at_pos in muts:
            try:
                filt_muts += [muts_at_pos.filter_against_normal(nmuts[muts_at_pos.pos])]
            except KeyError:
                # no normal coverage at that position
                filt_muts += [muts_at_pos]

        revmuts = []
        for muts_at_pos in filt_muts:
            revmuts += select_only_revertant_mutations(m, muts_at_pos)

        outfile.write("CHROM\tPOS\tID\tREF\tALT\tCOV\tCOUNT\tMAF\n")
        for revm in revmuts:
            outfile.write(revm.to_tsv() + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reffa", type=str, help="Reference genome (fasta)")
    parser.add_argument("mutations_tsv", type=str, help="tsv with mutation(s) to be reverted")
    parser.add_argument("search_bam", type=str, help="Find revertant mutations in these bams")
    parser.add_argument("normal_bam", type=str, help="Normal used to filter "
                        "mutations.")
    parser.add_argument("--version", action='version', version=revmut.__version__)
    args = parser.parse_args()
    find_revertant_mutations(args.reffa, args.mutations_tsv, args.search_bam,
                             args.normal_bam, sys.stdout)


if __name__ == "__main__":
    main()
