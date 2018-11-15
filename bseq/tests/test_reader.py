# -*- coding: utf-8 -*-
"""Nose tests for reader functions.
"""
import tempfile
from bseq.reader import read_fasta_file

class TestReadFastaFile:
    """Unit tests for read_fasta_file
    """
    def setup(self):
        self.fp = tempfile.NamedTemporaryFile(mode='w+')
        lines = ['>Test1 This is a description',
                 'ATGCATGCATGCAAA',
                 'ATGCATGCATGCAAA',
                 '>Test2 This is a another description',
                 'CATGCATGCAAATTT',
                 'CATGCATGCAAATTT',
                 '']
        self.fp.write('\n'.join(lines))
        self.fp.seek(0)
        self.path = self.fp.name

    def teardown(self):
        self.fp.close()

    def test_read_fasta_file(self):
        seq_list = read_fasta_file(self.path)
        assert len(seq_list) == 2


