# -*- coding: utf-8 -*-
"""Data structure models for biological alignments.

This module contains the classes that describe representations
of nucleic acid, protein, and codon alignments.

An Alignment is a list of sequences that have been aligned such that
characters in the same alignment column are evolutionary related to
each other. Similar to the Sequence object, an Alignment has a name
to identify the alignment, a sequence type to identify what king of
biological sequence is in the alignment, and optionally, a description
that contains related information regarding the alignment. 

The sequences in the Alignment object are not added at the time of
creating the object. Instead, sequences are added one by one through a
method. When reading from a FASTA file, a helper function called
`fasta_to_alignment` can be used to create the Alignment object
and add in the sequences automatically.

"""
from collections import namedtuple
from copy import deepcopy
import numpy as np
from bseq.sequence import Sequence, NuclSequence, ProtSequence, CodonSequence
from bseq.formatter import fasta_formatted_string


class Alignment(object):
    """Represents an alignment of biological sequences.

    Attributes
    ----------
    name : str
        Name of the alignment
    description : str
        Description of the alignment
    aln_type : str
        nucleotide, protein, or codon
    filters : list
        Contains the list of Filter objects associated to the
        Alignment instance.
    sequences
    i

    Notes
    -----
    The Alignment object is indexable using [] like a string or a list.
    If the key used is an object, it acts like a column-wise positional index.
    But if the key used is a string, it retrieves like a dictionary key and
    looks up the list of sequence names for a match.

    """
    def __init__(self, name, description=None, aln_type=None):
        """Creates a new alignment.

        Sequences cannot be added directly at instantiation.
        Use the `add_sequence_obj` or `add_sequence` methods to add
        sequences into the alignment.

        Parameters
        ----------
        name : str
            Name of the alignment.
        description : str, optional
            Description or other information about the alignment.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.
        aln_type : str, optional
            Alignment type indicates whether the alignment is a nucleotide,
            protein (amino acid), or codon alignment.

        """
        self.name = name
        self.description = description
        self.aln_type = aln_type
        self._records = []  # list of SequenceAnnotation objects
        self._records_lookup_d = dict()
        self._aln_matrix = np.array([])
        self.filters = dict()

    @property
    def sequences(self):
        """Returns a read-only copy of the alignment as a list of
        strings.
        """
        return [s for s in self]

    @property
    def i(self):
        """Returns a read-only copy of the alignment matrix.
        """
        return self._aln_matrix

    def add_sequence_obj(self, sequence_obj):
        """Adds a Sequence object containing a single aligned sequence
        to the alignment.

        Assumes that the sequence to be added is a Sequence object.
        If adding a sequence string, use the `add_sequence` method
        instead.

        Parameters
        ----------
        sequence_obj : Sequence

        See also
        --------
        add_sequence

        """
        assert sequence_obj.name not in self._records_lookup_d.keys()
        seq_annot = SequenceAnnotation(sequence_obj.name,
                                       sequence_obj.description,
                                       sequence_obj.seq_type)
        if not self.aln_type:
            self.aln_type = seq_annot.seq_type
        assert self.aln_type == seq_annot.seq_type
        self._records_lookup_d[seq_annot.name] = len(self._records)
        self._records.append(seq_annot)
        seq_array = np.array(list(sequence_obj.sequence))
        if self._aln_matrix.shape[-1] == 0:
            self._aln_matrix = np.array([seq_array,])
        else:
            self._aln_matrix = np.vstack((self._aln_matrix, seq_array))

    def add_sequence(self, name, sequence, seq_type, description=None):
        """Adds a single alignmed sequence to the alignment.

        In this method, the sequence is a string and not yet wrapped
        as a Sequence object. For adding sequence objects, use the
        `add_sequence_obj` method instead.

        Parameters
        ----------
        name : str
            Name of the sequence. This corresponds to the text found in the
            identifier line of a FASTA-formatted file.
        sequence : str
            Biological sequence.
        seq_type: str
            Indicates whether the sequence is a nucleotide, protein
            (amino acid), or codon sequence.
        description : str, optional
            Description or other miscellenous information about the sequence.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.

        See aslo
        --------
        add_sequence_obj

        """
        seq = Sequence(name, sequence,
                       description=description, seq_type=seq_type)
        self.add_sequence_obj(seq)

    def filter_sites(self, *filter_names, exclude_char='X'):
        """Filters out sites in the alignmnet using a given list of filters.

        Assumes that all filter tracks use the same excluding character to mark
        sites. The default excluding character is "X".

        Parameters
        ----------
        filter_names: str or Filter object
            Names of registered filters or Filter objects to be used
            to filter the alignment.
        exclude_char : str
            Filter track character that marks sites to be excluded
            from the returned alignment.

        Returns
        -------
        Alignment
            New alignment object with excluded sites removed from the
            alignment. The number of entries in the alignment does not change.

        See also
        --------
        use_all_filters
        filter_sequences

        """
        keep_coords = set(range(self._aln_matrix.shape[0]))
        for filt in filter_names:
            coords = set()
            if isinstance(filt, str) and filt in self.filters.keys():
                coords = set(self.filters[filt].coords(
                    exclude_char, inverse=True))
            elif isinstance(filt, Marker):
                coords = set(filt.coords(exclude_char, inverse=True))
            else:
                raise ValueError()
            keep_coords = keep_coords.intersection(coords)
        new_aln = deepcopy(self)
        new_aln._aln_matrix = self._aln_matrix[:, coords]  # pylint: disable=protected-access
        return new_aln

    def filter_sequences(self, *sequence_names):
        """Filters the alignmnet by returning only the specified entries.

        This does not change the number of sites in the alignment.

        Parameters
        ----------
        sequence_names
            Filter track character that marks sites to be excluded
            from the returned alignment

        Returns
        -------
        Alignment
            New alignment object with excluded sites removed from the
            alignment. The number of entries in the alignment does not change.

        See also
        --------
        filter_sites
        use_all_filters

        """
        positions = []
        new_records = []
        for name in sequence_names:
            if name in self._records_lookup_d.keys():
                i = self._records_lookup_d[name]
                positions.append(i)
                new_records.append(self._records[i])

        new_aln = deepcopy(self)
        new_aln._aln_matrix = self._aln_matrix[positions]  # pylint: disable=protected-access
        new_aln._records = new_records  # pylint: disable=protected-access
        new_aln._records_lookup_d = {k:v for k, v in self._records_lookup_d # pylint: disable=protected-access
                                     if k in sequence_names}
        return new_aln

    def use_all_filters(self, exclude_char='X'):
        """Filters the alignment using all filters associated with the
        current alignment.

        Parameters
        ----------
        exclude_char : str
            Character indicating the positions to be exclude from the
            final alignment.

        Returns
        -------
        Alignment
            New alignment object with excluded sites removed from the
            alignment. The number of entries in the alignment does not change.

        See also
        --------
        filter_sites
        filter_sequences

        """
        return self.filter_sites(*list(self.filters.keys()),
                                 exclude_char=exclude_char)

    def fasta_format(self, line_width=None):
        """Output the alignment as a FASTA-formatted string.

        Parameters
        ----------
        line_width : int
            Number of characters per line.

        Returns
        -------
        str
            FASTA formatted string such that
            `>{self.name}\n{self.seq}\n>{self.name}\n{self.seq}`.

        """
        fasta_string = ''
        for i, record in enumerate(self._records):
            fasta_string += \
                fasta_formatted_string(record.name,
                                       self._aln_matrix[i],
                                       description=record.description,
                                       line_width=line_width)
        return fasta_string

    def __len__(self):
        return self._aln_matrix.shape[-1]

    def __getitem__(self, i):
        if isinstance(i, int):
            return self._aln_matrix[:, i]
        elif isinstance(i, slice):
            return self._aln_matrix[:, i]
        elif isinstance(i, str):
            if i in self._records_lookup_d.keys():
                pos = self._records_lookup_d[i]
                return self._aln_matrix[pos]
        return IndexError()

    def __iter__(self):
        return iter(self._aln_matrix.transpose())

    # def __str__(self):
    #     pass
    #
    # def __repr__(self):
    #     pass


class NuclAlignment(Alignment):
    """Alignment composed of nucleotide sequences.

    Attributes
    ----------
    name : str
        Name of the alignment
    description : str
        Description of the alignment
    aln_type : str
        nucleotide, protein, or codon
    filters : list
        Contains the list of Filter objects associated to the
        Alignment instance.
    sequences
    i

    Notes
    -----
    The Alignment object is indexable using [] like a string or a list.
    If the key used is an object, it acts like a column-wise positional index.
    But if the key used is a string, it retrieves like a dictionary key and
    looks up the list of sequence names for a match.

    """
    def __init__(self, name, description=None):
        """Creates a new nucleotide alignment.

        Sequences cannot be added directly at instantiation.
        Use the `add_sequence_obj` or `add_sequence` methods to add
        sequences into the alignment.

        Parameters
        ----------
        name : str
            Name of the alignment.
        description : str, optional
            Description or other information about the alignment.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.

        """
        super().__init__(name, description=description, aln_type='nucleotide')

    def add_sequence_obj(self, sequence_obj: NuclSequence):
        """Adds a NuclSequence object containing a single aligned sequence
        to the alignment.

        Asserts that the sequence to be added is a NuclSequence object.
        If adding a sequence string, use the `add_sequence` method
        instead.

        Parameters
        ----------
        sequence_obj : NuclSequence

        See also
        --------
        add_sequence
        Alignment.add_sequence_obj

        """
        assert sequence_obj.seq_type == 'nucleotide'
        super().add_sequence_obj(sequence_obj)

    def add_sequence(self, name, sequence, seq_type, description=None):
        """Adds a single aligned sequence to the alignment.

        In this method, the sequence is a string and not yet wrapped
        as a Sequence object. For adding sequence objects, use the
        `add_sequence_obj` method instead.

        Parameters
        ----------
        name : str
            Name of the sequence. This corresponds to the text found in the
            identifier line of a FASTA-formatted file.
        sequence : str
            Biological sequence.
        seq_type: str
            Indicates whether the sequence is a nucleotide, protein
            (amino acid), or codon sequence.
        description : str, optional
            Description or other miscellenous information about the sequence.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.

        See aslo
        --------
        add_sequence_obj
        Alignment.add_sequence

        """
        assert seq_type == 'nucleotide'
        super().add_sequence(name, sequence, seq_type, description=description)


class ProtAlignment(Alignment):
    """Alignment composed of protein or amino acid sequences.

    Attributes
    ----------
    name : str
        Name of the alignment
    description : str
        Description of the alignment
    aln_type : str
        nucleotide, protein, or codon
    filters : list
        Contains the list of Filter objects associated to the
        Alignment instance.
    sequences
    i

    Notes
    -----
    The Alignment object is indexable using [] like a string or a list.
    If the key used is an object, it acts like a column-wise positional index.
    But if the key used is a string, it retrieves like a dictionary key and
    looks up the list of sequence names for a match.

    """
    def __init__(self, name, description=None):
        """Creates a new protein or amino acid alignment.

        Sequences cannot be added directly at instantiation.
        Use the `add_sequence_obj` or `add_sequence` methods to add
        sequences into the alignment.

        Parameters
        ----------
        name : str
            Name of the alignment.
        description : str, optional
            Description or other information about the alignment.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.

        """
        super().__init__(name, description=description, aln_type='protein')

    def add_sequence_obj(self, sequence_obj: NuclSequence):
        """Adds a ProtSequence object containing a single aligned sequence
        to the alignment.

        Asserts that the sequence to be added is a ProtSequence object.
        If adding a sequence string, use the `add_sequence` method
        instead.

        Parameters
        ----------
        sequence_obj : ProtSequence

        See also
        --------
        add_sequence
        Alignment.add_sequence_obj

        """
        assert sequence_obj.seq_type == 'protein'
        super().add_sequence_obj(sequence_obj)

    def add_sequence(self, name, sequence, seq_type, description=None):
        """Adds a single aligned sequence to the alignment.

        In this method, the sequence is a string and not yet wrapped
        as a Sequence object. For adding sequence objects, use the
        `add_sequence_obj` method instead.

        Parameters
        ----------
        name : str
            Name of the sequence. This corresponds to the text found in the
            identifier line of a FASTA-formatted file.
        sequence : str
            Biological sequence.
        seq_type: str
            Indicates whether the sequence is a nucleotide, protein
            (amino acid), or codon sequence.
        description : str, optional
            Description or other miscellenous information about the sequence.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.

        See aslo
        --------
        add_sequence_obj
        Alignment.add_sequence

        """
        assert seq_type == 'protein'
        super().add_sequence(name, sequence, seq_type, description=description)


class CodonAlignment(Alignment):
    """Alignment composed of coding sequences.

    Assumes that alignment is aligned by codon sites and not by nucleic acid.
    This means that the alignment length should be multiples of 3 and
    that gaps also in multiples of 3.

    Attributes
    ----------
    name : str
        Name of the alignment
    description : str
        Description of the alignment
    aln_type : str
        nucleotide, protein, or codon
    filters : list
        Contains the list of Filter objects associated to the
        Alignment instance.
    sequences
    i

    Notes
    -----
    The Alignment object is indexable using [] like a string or a list.
    If the key used is an object, it acts like a column-wise positional index.
    But if the key used is a string, it retrieves like a dictionary key and
    looks up the list of sequence names for a match.

    However, unlike the Alignment base class, length-related values of
    CodonAlignment are counted by the number of codons and not by the
    number of nucleotides. This means that for a alignment `a` that is
    10 codons or 30 nucleotides long, `len(s)` returns 10 instead of 30.
    Similarly, indexing is via codon and not by nucleotide position.

    """
    def __init__(self, name, description=None):
        """Creates a new codon alignment.

        Sequences cannot be added directly at instantiation.
        Use the `add_sequence_obj` or `add_sequence` methods to add
        sequences into the alignment.

        Parameters
        ----------
        name : str
            Name of the alignment.
        description : str, optional
            Description or other information about the alignment.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.

        """
        super().__init__(name, description=description, aln_type='codon')

    def add_sequence_obj(self, sequence_obj: CodonSequence):
        """Adds a CodonSequence object containing a single aligned sequence
        to the alignment.

        Asserts that the sequence to be added is a CodonSequence object.
        If adding a sequence string, use the `add_sequence` method
        instead.

        Parameters
        ----------
        sequence_obj : CodonSequence

        See also
        --------
        add_sequence
        Alignment.add_sequence_obj

        """
        assert sequence_obj.seq_type == 'codon'
        super().add_sequence_obj(sequence_obj)

    def add_sequence(self, name, sequence, seq_type, description=None):
        """Adds a single aligned sequence to the alignment.

        In this method, the sequence is a string and not yet wrapped
        as a Sequence object. For adding sequence objects, use the
        `add_sequence_obj` method instead.

        Parameters
        ----------
        name : str
            Name of the sequence. This corresponds to the text found in the
            identifier line of a FASTA-formatted file.
        sequence : str
            Biological sequence.
        seq_type: str
            Indicates whether the sequence is a nucleotide, protein
            (amino acid), or codon sequence.
        description : str, optional
            Description or other miscellenous information about the sequence.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.

        See aslo
        --------
        add_sequence_obj
        Alignment.add_sequence

        """
        assert seq_type == 'codon'
        super().add_sequence(name, sequence, seq_type, description=description)

    def __len__(self):
        return int(len(self._aln_matrix.shape[-1]) / 3)

    def __getitem__(self, i):
        if isinstance(i, int):  # self[0] returns the first alignment column
            x = int(i/3)
            return self._aln_matrix[:, x:x+3]
        elif isinstance(i, slice):  # self[0:2] returns the first 2 columns
            start = i.start * 3
            end = (i.stop * 3) + 3
            return self._aln_matrix[:, start:end]
        elif isinstance(i, str):  # self['test'] returns the sample's sequence
            if i in self._records_lookup_d.keys():
                pos = self._records_lookup_d[i]
                return self._aln_matrix[pos]
        raise TypeError('Alignment can only be indexed by integer, slice, \
                        or string.')

    def __iter__(self):
        return iter(self._aln_matrix)
