# -*- coding: utf-8 -*-
"""Nose tests for Alignment and its subclasses.
"""
from bseq.sequence import NuclSequence
from bseq.alignment import Alignment


class TestAlignmentEmpty:
    def setup(self):
        self.aln = Alignment('test')

    def test_add_sequence_obj(self):
        seq_obj = NuclSequence('test', 'ATGCATGCATGCAAA')
        self.aln.add_sequence_obj(seq_obj)
        assert len(self.aln._records) == 1  # pylint: disable=W0212
        assert len(self.aln._records_lookup_d) == 1  # pylint: disable=W0212
        assert len(self.aln._aln_matrix) == 1  # pylint: disable=W0212

    def test_add_sequence(self):
        self.aln.add_sequence('test', 'ATGCATGCATGCAAA', 'nucleotide')
        assert len(self.aln._records) == 1  # pylint: disable=W0212
        assert len(self.aln._records_lookup_d) == 1  # pylint: disable=W0212
        assert len(self.aln._aln_matrix) == 1  # pylint: disable=W0212

    def test_len(self):
        assert len(self.aln) == 0  # pylint: disable=C1801


class TestAlignment:
    def setup(self):
        self.aln = Alignment('test')
        self.aln.add_sequence('seq1', 'ATGCATGCATGCAAA', 'nucleotide')
        self.aln.add_sequence('seq2', 'ATGTATGCATGCAAA', 'nucleotide')
        self.aln.add_sequence('seq3', 'ATGCATGCATGCATA', 'nucleotide')
        self.aln.add_sequence('seq4', 'ATGCATGCATGCAAG', 'nucleotide')

    def test_len(self):
        assert len(self.aln) == 15
