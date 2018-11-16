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

    def test_coords(self):
        assert self.marker.coords('X') == [3, 10, 11, 12]
        assert self.marker.coords('O') == [0, 1, 2, 4, 5, 6, 7, 8, 9, 13, 14]
        assert self.marker.coords('X', inverse=True) == [0, 1, 2, 4, 5, 6,
                                                         7, 8, 9, 13, 14]
        assert self.marker.coords('X', explicit=False) == [[0, 3], [4, 10],
                                                           [13, 15]]
