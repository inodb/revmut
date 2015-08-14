from cStringIO import StringIO
import sys
import os
from nose.tools import ok_, assert_equals
import numpy as np
from os.path import join as ospj

FILE_PATH = os.path.realpath(__file__)
TEST_DIR_PATH = os.path.dirname(FILE_PATH)
DATA_PATH = os.path.abspath(ospj(TEST_DIR_PATH, "test_data"))
TMP_DIR_PATH = ospj(TEST_DIR_PATH, "nose_tmp_output")
TMP_BASENAME_DIR = ospj(TMP_DIR_PATH, "validation")
PKG_PATH = ospj(TEST_DIR_PATH, '..')

sys.path.append(PKG_PATH)
from revmut import revertant_mutation_checker
from revmut import utils

from Bio import SeqIO
import pyhgvs as hgvs


class TestRevertantMutation(object):
    def setUp(self):
        """Delete temporary dir if it exists then create it"""
        self.tearDown()
        utils.mkdir_p(TMP_BASENAME_DIR)

    def tearDown(self):
        """remove temp output files"""
        utils.rm_rf(TMP_DIR_PATH)

    def test_apply_hgvs(self):
        # p.Leu1303Phefs == c.3908dupT
        transcripts = \
            SeqIO.to_dict(SeqIO.parse("tests/test_data/BRCA1_transcripts.fa",
            "fasta"))
        brca1_mut = hgvs.HGVSName("ENST00000357654:c.3908dupT")
        normal_p = transcripts["ENST00000357654"].seq.translate()
        assert_equals("L", normal_p[1302])
        mut_c = revertant_mutation_checker.apply_hgvs(transcripts["ENST00000357654"].seq, brca1_mut)
        assert_equals("TT", mut_c[3907:3909])
        mut_p = mut_c.translate()
        assert_equals("F", mut_p[1302])

    def test_revertant_mutation_checker_oncotator_del(self):
        out = StringIO()
        revertant_mutation_checker.print_revertant_mutations_info(
            ospj(DATA_PATH, "to_be_reverted_mutations.txt"),
            ospj(DATA_PATH, "oncotator.del.maf.txt"),
            ospj(DATA_PATH, "BRCA_transcripts.fa"),
            revmuts_file_format='oncotator',
            outfile=out
        )
        assert_equals(open(ospj(DATA_PATH, "output", "oncotator.del.maf.out.tsv")).read(), out.getvalue())

    def test_revertant_mutation_checker_oncotator_ins(self):
        out = StringIO()
        revertant_mutation_checker.print_revertant_mutations_info(
            ospj(DATA_PATH, "to_be_reverted_mutations.txt"),
            ospj(DATA_PATH, "oncotator.ins.maf.txt"),
            ospj(DATA_PATH, "BRCA_transcripts.fa"),
            revmuts_file_format='oncotator',
            outfile=out
        )
        assert_equals(open(ospj(DATA_PATH, "output", "oncotator.ins.maf.out.tsv")).read(), out.getvalue())

    def test_revertant_mutation_checker_ins(self):
        out = StringIO()
        revertant_mutation_checker.print_revertant_mutations_info(
            ospj(DATA_PATH, "to_be_reverted_mutations.txt"),
            ospj(DATA_PATH, "oncotator.ins.txt"),
            ospj(DATA_PATH, "BRCA_transcripts.fa"),
            revmuts_file_format='hgvs',
            outfile=out
        )
        assert_equals(open(ospj(DATA_PATH, "output", "oncotator.ins.maf.out.tsv")).read(), out.getvalue())

    def test_revertant_mutation_checker_del(self):
        out = StringIO()
        revertant_mutation_checker.print_revertant_mutations_info(
            ospj(DATA_PATH, "to_be_reverted_mutations.txt"),
            ospj(DATA_PATH, "oncotator.del.txt"),
            ospj(DATA_PATH, "BRCA_transcripts.fa"),
            revmuts_file_format='hgvs',
            outfile=out
        )
        assert_equals(open(ospj(DATA_PATH, "output", "oncotator.del.maf.out.tsv")).read(), out.getvalue())
