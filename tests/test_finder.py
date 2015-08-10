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
from revmut.find_mutations import find_revertant_mutations
from revmut import utils


class TestRevertantMutation(object):
    def setUp(self):
        """Delete temporary dir if it exists then create it"""
        self.tearDown()
        utils.mkdir_p(TMP_BASENAME_DIR)

    def tearDown(self):
        """remove temp output files"""
        utils.rm_rf(TMP_DIR_PATH)

    def test_revertant_mutation_finder(self):
        out = StringIO()

        reffa = ospj(DATA_PATH, "human_g1k_v37_chr17.fa")
        mutations_tsv = ospj(DATA_PATH, "germline_mutations", "T1_test_mutation.tsv")
        search_bam = ospj(DATA_PATH, "T1.bam")
        normal_bam = ospj(DATA_PATH, "N1.bam")

        find_revertant_mutations(reffa, mutations_tsv, search_bam, normal_bam, out)
        assert_equals(open(ospj(DATA_PATH, "output", "T1_test.tsv")).read(), out.getvalue())
