# -*- coding: utf-8 -*-
"""Data structure models for biological sequences.

This module contains the classes that describe representations
of nucleic acid, protein, and codon sequences.

A sequence is a series of characters that represent a biological
sequence of bases or amino acids. The Sequence object is meant to
represent these biological molecules. A Sequence object has a name
to identify the sequence, the string representation of the sequence,
a sequence type to identify what king of biological sequence the
string sequence is representing, and optionally, a description that
contains related information regarding the sequence. However, the
Sequence object is not meant to be used directly. Instead, it is
encouraged to use the subclasses of Sequence: NuclSequence for
nucleotide sequences, ProtSequence for protein/amino acid sequences,
and CodonSequence for coding sequences.

"""
from collections import Counter
import numpy as np
from bseq.formatter import fasta_formatted_string


class Sequence(object):
    """Represents a biological sequence of characters.

    This object can be used to represent a singleton sequence,
    a sequence that is part of a multi-sequence file, or
    one sequence in an alignment.

    Attributes
    ----------
    name : str
        Name of the sequence.
        This corresponds to the ID in the FASTA format.
    description : str
        Description and other information about the sequence.
        This corresponds to the text after the ID, separated by a whitespace, in the FASTA format.
    seq_type : str
        nucleotide, protein, or codon
    sequence : str
        Biological sequence

    See also
    --------
    NuclSequence
    ProtSequence
    CodonSequence

    Notes
    -----
    The Sequence object is indexable using [] and can be
    sliced like a string or a list.

    """
    def __init__(self, name, sequence, description=None, seq_type=None):
        """Creates a new instance of the Sequence object.

        Parameters
        ----------
        name : str
            Name of the sequence. This corresponds to the text found in the
            identifier line of a FASTA-formatted file.
        sequence : str
            Biological sequence.
        description : str, optional
            Description or other miscellenous information about the sequence.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.
        seq_type: str, optional
            Indicates whether the sequence is a nucleotide, protein
            (amino acid), or codon sequence.

        """
        self.name = name
        self.description = description
        self.seq_type = seq_type
        self._sequence = sequence

    @property
    def i(self):
        """Always returns a read-only array of single characters,
        even for codon sequences.
        """
        return np.array(self._sequence)

    @property
    def sequence(self):
        """Returns a read-only string representation of the
        sequence.
        """
        return self._sequence

    def fasta_format(self, line_width=None):
        """Output the sequence as a FASTA-formatted string.

        Parameters
        ----------
        line_width : int
            Number of characters per line.

        Returns
        -------
        str
            FASTA formatted string such that `>{self.name}\n{self.seq}`.

        """
        return fasta_formatted_string(self.name, self._sequence,
                                      description=self.description,
                                      line_width=line_width)

    def count(self, char):
        """Counts the number of times a given character occurs in the sequence.

        Parameters
        ----------
        char : str
            Character to count in the sequence

        Returns
        -------
        int
            Number of occurrences in the sequence

        """
        return self._sequence.count(char)

    def count_all(self):
        """Counts the number of times each character occurs in the sequence.

        Returns
        -------
        Counter
            Keys are the characters present in the sequence and
            values are the corresponding number of occurrences in the sequence

        """
        return Counter(self._sequence)

    def __len__(self):
        return len(self._sequence)

    def __getitem__(self, i):
        return self._sequence[i]

    def __iter__(self):
        return iter(self._sequence)

    def __contains__(self, x):
        return x in self._sequence

    def __str__(self):
        return self._sequence


class NuclSequence(Sequence):
    """Represents a nucleotide sequence.

    Attributes
    ----------
    name : str
        Name of the sequence.
        This corresponds to the ID in the FASTA format.
    description : str
        Description and other information about the sequence.
        This corresponds to the text after the ID, separated by a whitespace, in the FASTA format.
    seq_type : str
        nucleotide, protein, or codon
    sequence : str
        Biological sequence

    See also
    --------
    Sequence
    ProtSequence
    CodonSequence

    Notes
    -----
    The Sequence object is indexable using [] and can be
    sliced like a string or a list.

    """
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description,
                         seq_type='nucleotide')

class ProtSequence(Sequence):
    """Represents a protein or amino acid sequence.

    Attributes
    ----------
    name : str
        Name of the sequence.
        This corresponds to the ID in the FASTA format.
    description : str
        Description and other information about the sequence.
        This corresponds to the text after the ID, separated by a whitespace, in the FASTA format.
    seq_type : str
        nucleotide, protein, or codon
    sequence : str
        Biological sequence

    See also
    --------
    Sequence
    NuclSequence
    CodonSequence

    Notes
    -----
    The Sequence object is indexable using [] and can be
    sliced like a string or a list.

    """
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description,
                         seq_type='protein')

class CodonSequence(Sequence):
    """Represents a codon-based coding sequence.

    Attributes
    ----------
    name : str
        Name of the sequence.
        This corresponds to the ID in the FASTA format.
    description : str
        Description and other information about the sequence.
        This corresponds to the text after the ID, separated by a whitespace, in the FASTA format.
    seq_type : str
        nucleotide, protein, or codon
    sequence : str
        Biological sequence

    See also
    --------
    Sequence
    NuclSequence
    ProtSequence

    Notes
    -----
    Unlike the base class Sequence, length-related values of CodonSequence are
    counted by the number of codons and not by the number of nucleotides in the
    sequence. This means that for a sequence `s` of 10 codons or 30 nucleotides,
    `len(s)` returns 10 instead of 30. Similarly, indexing is via codon
    and notby nucleotide position.

    """
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description,
                         seq_type='codon')

    def count_codon(self, codon):
        """Counts the number of times a given codon occurs in the sequence.

        Parameters
        ----------
        codon : str
            Codon to count in the sequence

        Returns
        -------
        int
            Number of occurrences in the sequence

        """
        return sum([1 for c in self if c == codon])

    def count_codon_all(self):
        """Counts the number of times each character occurs in the sequence.

        Returns
        -------
        Counter
            Keys are the characters present in the sequence and
            values are the corresponding number of occurrences in the sequence

        """
        return Counter(list(self))

    def __len__(self):
        return int(len(self._sequence) / 3)

    def __getitem__(self, i):
        if isinstance(i, int):
            if i >= 0:
                i = i * 3
                return self._sequence[i:i+3]
            i = len(self._sequence) + 1 - 1
            return self._sequence[i-3:i]
        elif isinstance(i, slice):
            start = i.start * 3
            end = (i.stop * 3) if isinstance(i.stop, int) \
                  else len(self._sequence)
            return self._sequence[start:end]
        return IndexError()

    def __iter__(self):
        return (self._sequence[i:i+3] for i in range(0, len(self._sequence), 3))

    def __str__(self):
        return ' '.join(list(self))
