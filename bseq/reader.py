# -*- coding: utf-8 -*-
"""Helper functions for reading files.
"""
from bseq.sequence import NuclSequence, ProtSequence, CodonSequence
from bseq.alignment import NuclAlignment, ProtAlignment, CodonAlignment


def read_fasta_file(path, seq_type='nucleotide'):
    """Reads the FASTA file and stores the contents
    as a list of Sequence objects.

    Parameters
    ----------
    path : str
        Path to the FASTA file.
    seq_type : str, optional
        Type of sequence expected. Choices are 'nucleotide',
        'protein', or 'codon'.

    Returns
    -------
    list of Sequence objects

    """
    name = ''
    description = ''
    seq = []
    sequence_list = []
    with open(path, 'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                if seq:
                    seq_string = ''.join(seq)
                    seq_obj = None
                    if seq_type == 'nucleotide':
                        seq_obj = NuclSequence(name, seq_string, description)
                    elif seq_type == 'protein':
                        seq_obj = ProtSequence(name, seq_string, description)
                    elif seq_type == 'codon':
                        seq_obj = CodonSequence(name, seq_string, description)
                    else:
                        raise ValueError('seq_type must be \
                                          "nucleotide", "protein", \
                                          or "codon".')
                    sequence_list.append(seq_obj)
                    seq = []
                name, description = line[1:].rstrip().split(' ', 1)
            else:
                seq += [line.rstrip()]

        if seq:
            seq_string = ''.join(seq)
            seq_obj = None
            if seq_type == 'nucleotide':
                seq_obj = NuclSequence(name, seq_string, description)
            elif seq_type == 'protein':
                seq_obj = ProtSequence(name, seq_string, description)
            elif seq_type == 'codon':
                seq_obj = CodonSequence(name, seq_string, description)
            else:
                raise ValueError('seq_type must be \
                                    "nucleotide", "protein", \
                                    or "codon".')
            sequence_list.append(seq_obj)
    return sequence_list

def read_fasta_alignment(path, seq_type='nucleotide',
                         name=None, description=None):
    """Reads the FASTA alignment and stores the contents
    as an Alignment object.

    Parameters
    ----------
    path : str
        Path to the FASTA file.
    seq_type : str, optional
        Type of sequence expected. Choices are 'nucleotide',
        'protein', or 'codon'.

    Returns
    -------
    Alignment

    """
    alignment = None
    if seq_type == 'nucleotide':
        alignment = NuclAlignment(name, description)
    elif seq_type == 'protein':
        alignment = ProtAlignment(name, description)
    elif seq_type == 'codon':
        alignment = CodonAlignment(name, description)

    for seq_obj in read_fasta_file(path, seq_type=seq_type):
        alignment.add_sequence_obj(seq_obj)

    return alignment
