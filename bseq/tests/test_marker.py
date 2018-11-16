# -*- coding: utf-8 -*-
"""Nose tests for Alignment and its subclasses.
"""
from bseq.marker import Marker


class TestMarker:
    def setup(self):
        self.marker = Marker('test', {'O': 'keep', 'X': 'remove'},
                             'OOOXOOOOOOXXXOO')

    def test_sequence(self):
        assert self.marker.sequence == 'OOOXOOOOOOXXXOO'

    def test_encoded_sequence(self):
        assert self.marker.encoded_sequence == '0O3X4O10X13O15'

    def test_encode(self):
        self.marker._encode('XOOXOOOOOOXXXOX')
        assert self.marker._char_list == ('X', 'O', 'X', 'O', 'X', 'O', 'X')
        assert self.marker._pos_list == (
            (0, 1), (1, 3), (3, 4), (4, 10), (10, 13), (13, 14), (14, 15)
        )

    def test_coords(self):
        assert self.marker.coords('X') == [3, 10, 11, 12]
        assert self.marker.coords('O') == [0, 1, 2, 4, 5, 6, 7, 8, 9, 13, 14]
        assert self.marker.coords('X', inverse=True) == [0, 1, 2, 4, 5, 6,
                                                         7, 8, 9, 13, 14]
        assert self.marker.coords('X', explicit=False) == [[3, 4], [10, 13]]
        assert self.marker.coords('O', explicit=False) == [[0, 3], [4, 10],
                                                           [13, 15]]

    def test_filter(self):
        sequence = 'ATTCAATATACCCAT'
        assert self.marker.filter(sequence, 'X') == 'ATTAATATAAT'
        assert self.marker.filter(sequence, 'O') == 'CCCC'
        assert self.marker.filter(sequence, 'X', output_coords=True) == \
            [0, 1, 2, 4, 5, 6, 7, 8, 9, 13, 14]
        assert self.marker.filter(sequence, 'O', output_coords=True) == \
            [3, 10, 11, 12]

    def test_mask(self):
        sequence = 'ATTCAATATACCCAT'
        assert self.marker.mask(sequence, 'X') == 'ATT_AATATA___AT'
        assert self.marker.mask(sequence, 'O') == '___C______CCC__'
        assert self.marker.mask(sequence, 'X', mask_char='.') == \
        'ATT.AATATA...AT'
