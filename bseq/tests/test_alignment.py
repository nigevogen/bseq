# -*- coding: utf-8 -*-
"""Nose tests for Alignment and its subclasses.
"""
from bseq.sequence import NuclSequence
from bseq.alignment import Alignment
from bseq.marker import Marker
import numpy as np


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

    def test_add_markers(self):
        # Add two sequence, polymorphic on 4th site
        self.aln.add_sequence('test1', 'ATGCATGCATGCAAA', 'nucleotide')
        self.aln.add_sequence('test2', 'ATGGATGCATGCAAA', 'nucleotide')
        self.aln.add_markers(
            Marker('test_marker', {'O':'conserved', 'X':'polymorphic'},
                   'OOOXOOOOOOOOOOO')
        )
        assert len(self.aln.markers) == 1
        assert 'test_marker' in self.aln.markers.keys()
        assert self.aln.markers['test_marker'].sequence == 'OOOXOOOOOOOOOOO'


class TestAlignment:
    def setup(self):
        self.aln = Alignment('test')
        self.aln.add_sequence('seq1', 'ATGCATGCATGCAAA', 'nucleotide')
        self.aln.add_sequence('seq2', 'ATGTATGCATGCAAA', 'nucleotide')
        self.aln.add_sequence('seq3', 'ATGCATGCATGCATA', 'nucleotide')
        self.aln.add_sequence('seq4', 'ATGCATGCATGCAAG', 'nucleotide')

    def test_len(self):
        assert len(self.aln) == 15

    def test_getitem(self):
        assert list(self.aln[0]) == list('AAAA')
        assert list(self.aln[0:1]) == list('AAAA')
        assert list(map(list, self.aln[0:2])) == [['A', 'T'], ['A', 'T'],
                                                  ['A', 'T'], ['A', 'T']]
        assert list(self.aln['seq2']) == list('ATGTATGCATGCAAA')
