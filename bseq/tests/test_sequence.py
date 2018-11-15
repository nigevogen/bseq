# -*- coding: utf-8 -*-
"""Nose tests for Sequence and its subclasses.
"""
from bseq.sequence import Sequence, NuclSequence, CodonSequence


class TestSequence:
    """Unit tests for Sequence.
    """
    def setup(self):
        self.seq = Sequence('test', 'ATGCATGCATGCAAA')

    def teardown(self):
        pass

    def test_fasta_format(self):
        # No line width specified
        fasta_string = self.seq.fasta_format()
        assert fasta_string == '>test\nATGCATGCATGCAAA\n'

        # Line width is shorter than the sequence.
        # Sequence should wrap according to the line width while the ID
        # line remains unaffected.
        fasta_string = self.seq.fasta_format(line_width=6)
        assert fasta_string == '>test\nATGCAT\nGCATGC\nAAA\n'

        # Line width is longer than the sequence.
        # No wrapper is expected.
        fasta_string = self.seq.fasta_format(line_width=20)
        assert fasta_string == '>test\nATGCATGCATGCAAA\n'

        self.seq.description = 'test description'
        fasta_string = self.seq.fasta_format()
        assert fasta_string == '>test test description\nATGCATGCATGCAAA\n'

    def test_len(self):
        assert len(self.seq) == 15

    def test_getitem(self):
        assert self.seq[0] == 'A'
        assert self.seq[10] == 'G'
        assert self.seq[14] == 'A'
        assert self.seq[-1] == 'A'

        assert self.seq[0:1] == 'A'
        assert self.seq[0:3] == 'ATG'
        assert self.seq[12:] == 'AAA'
        assert self.seq[-3:] == 'AAA'
        assert self.seq[-6:-3] == 'TGC'

    def test_contains(self):
        assert 'ATGCAT' in self.seq
        assert 'ATGCATGCATGCAAA' in self.seq
        assert 'ATGCGG' not in self.seq

    def test_str(self):
        assert str(self.seq) == 'ATGCATGCATGCAAA'

    def test_count(self):
        assert self.seq.count('A') == 6
        assert self.seq.count('T') == 3
        assert self.seq.count('C') == 3
        assert self.seq.count('G') == 3

    def test_count_all(self):
        counter = self.seq.count_all()
        assert len(
            set(counter.keys()).symmetric_difference(set(list('ATCG')))
        ) == 0
        assert counter['A'] == 6
        assert counter['T'] == 3
        assert counter['C'] == 3
        assert counter['G'] == 3


class TestNuclSequence:
    """Unit test for NuclSequence.

    NuclSequence is a subclass of Sequence.
    """
    def setup(self):
        self.seq = NuclSequence('test', 'ATGCATGCATGCAAA')
    
    def test_len(self):
        assert len(self.seq) == 15

    def test_getitem(self):
        assert self.seq[0] == 'A'
        assert self.seq[1] == 'T'
        assert self.seq[4] == 'A'
        assert self.seq[-1] == 'A'

        assert self.seq[0:1] == 'A'
        assert self.seq[0:2] == 'AT'
        assert self.seq[3:4] == 'C'
        assert self.seq[3:] == 'CATGCATGCAAA'
        assert self.seq[-2:] == 'AA'
        assert self.seq[-4:-1] == 'CAA'

    def test_contains(self):
        assert 'ATGCAT' in self.seq
        assert 'ATGCATGCATGCAAA' in self.seq
        assert 'ATGCGG' not in self.seq

    def test_str(self):
        assert str(self.seq) == 'ATGCATGCATGCAAA', print(str(self.seq))

    def test_count(self):
        assert self.seq.count('A') == 6
        assert self.seq.count('T') == 3
        assert self.seq.count('C') == 3
        assert self.seq.count('G') == 3

    def test_count_all(self):
        counter = self.seq.count_all()
        assert len(
            set(counter.keys()).symmetric_difference(set(list('ATCG')))
        ) == 0
        assert counter['A'] == 6
        assert counter['T'] == 3
        assert counter['C'] == 3
        assert counter['G'] == 3

class TestProtSequence:
    """Unit test for ProtSequence.

    ProtSequence is a subclass of Sequence.
    """
    def setup(self):
        self.seq = NuclSequence('test', 'VCWMMYDCGVVEIDC')
    
    def test_len(self):
        assert len(self.seq) == 15

    def test_getitem(self):
        assert self.seq[0] == 'V'
        assert self.seq[1] == 'C'
        assert self.seq[4] == 'M'
        assert self.seq[-1] == 'C'

        assert self.seq[0:1] == 'V'
        assert self.seq[0:2] == 'VC'
        assert self.seq[3:4] == 'M'
        assert self.seq[3:] == 'MMYDCGVVEIDC'
        assert self.seq[-2:] == 'DC'
        assert self.seq[-4:-1] == 'EID'

    def test_contains(self):
        assert 'VCWMMY' in self.seq
        assert 'VCWMMYDCGVVEIDC' in self.seq
        assert 'VCWMQY' not in self.seq

    def test_str(self):
        assert str(self.seq) == 'VCWMMYDCGVVEIDC', print(str(self.seq))

    def test_count(self):
        assert self.seq.count('V') == 3
        assert self.seq.count('C') == 3
        assert self.seq.count('W') == 1
        assert self.seq.count('M') == 2
        assert self.seq.count('Y') == 1
        assert self.seq.count('D') == 2
        assert self.seq.count('G') == 1
        assert self.seq.count('E') == 1
        assert self.seq.count('I') == 1


    def test_count_all(self):
        counter = self.seq.count_all()
        assert len(
            set(counter.keys()).symmetric_difference(set(list('VCWMYDGEI')))
        ) == 0
        assert counter['V'] == 3
        assert counter['C'] == 3
        assert counter['W'] == 1
        assert counter['M'] == 2
        assert counter['Y'] == 1
        assert counter['D'] == 2
        assert counter['G'] == 1
        assert counter['E'] == 1
        assert counter['I'] == 1

class TestCodonSequence:
    """Unit test for CodonSequence.

    CodonSequence is a subclass of Sequence.
    """
    def setup(self):
        self.seq = CodonSequence('test', 'ATGCATGCATGCAAA')

    def test_len(self):
        assert len(self.seq) == 5

    def test_getitem(self):
        assert self.seq[0] == 'ATG'
        assert self.seq[1] == 'CAT'
        assert self.seq[4] == 'AAA'
        assert self.seq[-1] == 'AAA'

        assert self.seq[0:1] == 'ATG'
        assert self.seq[0:2] == 'ATGCAT'
        assert self.seq[3:4] == 'TGC'
        assert self.seq[3:] == 'TGCAAA'
        assert self.seq[-2:] == 'TGCAAA'
        assert self.seq[-3:-1] == 'GCATGC'

    def test_contains(self):
        assert 'ATGCAT' in self.seq
        assert 'ATGCATGCATGCAAA' in self.seq
        assert 'ATGCGG' not in self.seq

    def test_str(self):
        assert str(self.seq) == 'ATG CAT GCA TGC AAA', print(str(self.seq))
