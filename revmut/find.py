"""
Find mutations from those found in exome sequenced samples using samtools mpileup and bcftools
"""
import argparse
import pandas as pd
import sys

import revmut

from sufam.mpileup_parser import run_and_get_mutations
from sufam.mutation import get_mutation


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
