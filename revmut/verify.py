"""
Check if a mutation is reverted
"""
import argparse
import sys
import revmut
import logging
from revmut import oncotator
from revmut.utils import deprecated

import pyhgvs as hgvs
from Bio import SeqIO


def minimumEditDistance(s1, s2):
    """Levenshtein distance between two strings from
    http://rosettacode.org/wiki/Levenshtein_distance"""
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1) + 1)
    for index2, char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1, char1 in enumerate(s1):
            if char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
        distances = newDistances
    return distances[-1]


class RevertantMutationsInfo(object):
    def __init__(self, transcript_record, hgvs_mut):
        self.record = transcript_record
        self.mut = hgvs_mut
        self.normal_p = self.record.seq.translate(to_stop=True)
        self.mut_p = apply_hgvs(self.record.seq, hgvs_mut).translate(to_stop=True)
        self.mut_p_distance = minimumEditDistance(self.normal_p, self.mut_p)
        self.revmuts = []
        self.revmuts_pos_adj = []
        self.revmuts_p = []
        self.revmut_p_distances = []

    def add_revmut(self, hgvs_revmut):
        assert(len(self.revmuts) == len(self.revmuts_p))
        self.revmuts += [hgvs_revmut]
        hgvs_revmut_pos_adj = alter_coords_hgvs_sequential(self.mut, hgvs_revmut)
        self.revmuts_pos_adj += [hgvs_revmut_pos_adj]
        self.revmuts_p += [apply_hgvs(apply_hgvs(self.record.seq, self.mut), hgvs_revmut_pos_adj).translate(to_stop=True)]
        self.revmut_p_distances += [minimumEditDistance(self.normal_p, self.revmuts_p[-1])]

    def to_tsv(self, header=True):
        assert(len(self.revmuts) == len(self.revmuts_p))
        rv = ""
        if header:
            rv += "mut\trevmut\trevmut_pos_adj\ttranscript\tnormal_protein_length\tmut_protein_length\trevmut_protein_length\tmut_p_distance\trevmut_p_distance\n"
        rv += "\n".join(["\t".join(["{}"]*9).format(
                         self.mut.name.split(":")[1],
                         self.revmuts[i].name.split(":")[1],
                         self.revmuts_pos_adj[i].name.split(":")[1],
                         self.record.id,
                         len(self.normal_p), len(self.mut_p),
                         len(self.revmuts_p[i]),
                         self.mut_p_distance,
                         self.revmut_p_distances[i]
                         ) for i in
                         range(len(self.revmuts))]) + "\n"
        return rv

    def __len__(self):
        return len(self.revmuts)


def alter_coords_hgvs_sequential(h1, h2):
    """Change HGVS coords of h2 after applying h1"""
    if h1.kind == "c" and h2.kind == "c":
        if h1.mutation_type == ">":
            h3 = hgvs.HGVSName(h2.name)
        elif h1.mutation_type == "del":
            if h1.cdna_start.coord > h2.cdna_end.coord:
                h3 = hgvs.HGVSName(h2.name)
            elif h1.cdna_end.coord < h2.cdna_start.coord:
                h3 = hgvs.HGVSName(h2.name)
                h3.cdna_start = hgvs.CDNACoord(coord=h3.cdna_start.coord-len(h1.ref_allele))
                h3.cdna_end = hgvs.CDNACoord(coord=h3.cdna_end.coord-len(h1.ref_allele))
            else:
                raise(Exception("Overlapping del not implemented"))
        elif h1.mutation_type == "ins":
            if h1.cdna_start.coord > h2.cdna_end.coord:
                h3 = hgvs.HGVSName(h2.name)
            elif h1.cdna_end.coord < h2.cdna_start.coord:
                h3 = hgvs.HGVSName(h2.name)
                h3.cdna_start = hgvs.CDNACoord(coord=h3.cdna_start.coord+len(h1.alt_allele))
                h3.cdna_end = hgvs.CDNACoord(coord=h3.cdna_end.coord+len(h1.alt_allele))
            else:
                raise(Exception("Overlapping ins not implemented"))
        elif h1.mutation_type == "dup":
            if h1.cdna_start.coord > h2.cdna_end.coord:
                h3 = hgvs.HGVSName(h2.name)
            elif h1.cdna_end.coord < h2.cdna_start.coord:
                h3 = hgvs.HGVSName(h2.name)
                h3.cdna_start = hgvs.CDNACoord(coord=h3.cdna_start.coord+len(h1.alt_allele)-len(h1.ref_allele))
                h3.cdna_end = hgvs.CDNACoord(coord=h3.cdna_end.coord+len(h1.alt_allele)-len(h1.ref_allele))
            else:
                raise(Exception("Overlapping dup not implemented"))
        elif h1.mutation_type == "delins":
            if h1.cdna_start.coord > h2.cdna_end.coord:
                h3 = hgvs.HGVSName(h2.name)
            elif h1.cdna_end.coord < h2.cdna_start.coord:
                h3 = hgvs.HGVSName(h2.name)
                h3.cdna_start = hgvs.CDNACoord(coord=h3.cdna_start.coord-len(h1.ref_allele)+len(h1.alt_allele))
                h3.cdna_end = hgvs.CDNACoord(coord=h3.cdna_end.coord-len(h1.ref_allele)+len(h1.alt_allele))
            else:
                raise(Exception("Overlapping delins not implemented"))
        else:
            raise(Exception("Unexpected mutation_type {}".format(h1.mutation_type)))
        return h3
    else:
        raise(Exception("Only cDNA mutations have been implemented"))


def apply_hgvs(seq, h):
    """Apply 1-based HGVS mutation to 0-based biopython Sequence"""
    if h.kind == "c":
        start = h.cdna_start.coord - 1
        end = h.cdna_end.coord - 1
        if h.mutation_type == ">":
            assert(seq[start] == h.ref_allele)
            new_seq = seq[:start] + h.alt_allele + seq[start+1:]
        elif h.mutation_type == "del":
            assert(seq[start:end+1] == h.ref_allele)
            new_seq = seq[:start] + seq[end+1:]
            assert(len(seq) == len(new_seq) + len(h.ref_allele))
        elif h.mutation_type == "ins":
            new_seq = seq[:start+1] + h.alt_allele + seq[end:]
            assert(len(seq) + len(h.alt_allele) == len(new_seq))
        elif h.mutation_type == "dup":
            assert(seq[start:end+1] == h.ref_allele)
            new_seq = seq[:start] + h.alt_allele + seq[end+1:]
            assert(len(seq) - len(h.ref_allele) + len(h.alt_allele) == len(new_seq))
        elif h.mutation_type == "delins":
            assert(seq[start:end+1] == h.ref_allele)
            new_seq = seq[:start] + h.alt_allele + seq[end+1:]
            assert(len(seq) - len(h.ref_allele) + len(h.alt_allele) == len(new_seq))
        else:
            raise(Exception("Unexpected mutation_type {}".format(h.mutation_type)))
        return new_seq
    else:
        raise(Exception("Only cDNA mutations have been implemented"))


@deprecated
def is_revertant(record, hgvs_mut, hgvs_rev_mut):
    normal_p = record.seq.translate(to_stop=True)
    mut_p = apply_hgvs(record.seq, hgvs_mut).translate(to_stop=True)
    revmut_p = apply_hgvs(apply_hgvs(record.seq, hgvs_mut), hgvs_rev_mut).translate(to_stop=True)
    logging.info("---")
    logging.info("mutation: {}".format(hgvs_mut))
    logging.info("putative revertant mutation: {}".format(hgvs_rev_mut))
    logging.info("transcript: {}".format(record.id))
    logging.info("normal length: {}".format(len(normal_p)))
    logging.info("mutation length: {}".format(len(mut_p)))
    logging.info("putative revertant length: {}".format(len(revmut_p)))
    logging.info("---")
    return len(normal_p) == len(revmut_p)


def parse_hgvs_muts_file(hgvs_muts_file, raise_exception=True):
    hgvs_muts = []
    for m in [l.rstrip('\n') for l in open(hgvs_muts_file).readlines()]:
        if raise_exception:
            hgvs_mut = hgvs.HGVSName(m)
            hgvs_muts += [hgvs_mut]
        else:
            try:
                hgvs_mut = hgvs.HGVSName(m)
                hgvs_muts += [hgvs_mut]
            except hgvs.InvalidHGVSName:
                sys.stderr.write("Invalid HGVS found: {}\n".format(m))
                pass
    return hgvs_muts


def print_revertant_mutations_info(hgvs_muts_file, revmuts_file, fasta, revmuts_file_format='hgvs', outfile=sys.stdout):
    assert(revmuts_file_format in ['hgvs', 'oncotator'])

    hgvs_muts = parse_hgvs_muts_file(hgvs_muts_file)
    if revmuts_file_format == 'hgvs':
        hgvs_revmuts = parse_hgvs_muts_file(revmuts_file)
    elif revmuts_file_format == 'oncotator':
        ot = oncotator.Oncotator(revmuts_file)
    transcripts = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    print_header = True
    for m in hgvs_muts:
        rmi = RevertantMutationsInfo(transcripts[m.transcript], m)

        # get only mutations for given transcript, version id postfix is ignored
        if revmuts_file_format == 'oncotator':
            hgvs_revmuts_t = ot.get_hgvs_mutations(m.transcript)
        else:
            hgvs_revmuts_t = [rm for rm in hgvs_revmuts
                if m.transcript == rm.transcript or
                   m.transcript == ".".join(rm.transcript.split(".")[:-1])]

        for rm in hgvs_revmuts_t:
            rmi.add_revmut(rm)

        if len(rmi) > 0:
            outfile.write(rmi.to_tsv(header=print_header))
            print_header = False


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        # format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
    )
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hgvs_muts_file", type=str, help="File with one mutation per line in HGVS format TRANSCRIPT_ID:c.HGVS_MUT")
    parser.add_argument("revmuts_file", type=str, help="File with one putative revertant mutation per line in HGVS format TRANSCRIPT_ID:c.HGVS_MUT")
    parser.add_argument("fasta", type=str, help="Fasta file with transcripts")
    parser.add_argument("--version", action='version', version=revmut.__version__)
    parser.add_argument("--revmuts_file_format", type=str, choices=["hgvs", "oncotator"], default="hgvs", help="Set revmuts_file format")
    args = parser.parse_args()
    print_revertant_mutations_info(args.hgvs_muts_file, args.revmuts_file, args.fasta, revmuts_file_format=args.revmuts_file_format, outfile=sys.stdout)


if __name__ == "__main__":
    main()
