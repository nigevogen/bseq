# -*- coding: utf-8 -*-
"""Helper functions for reading files.
"""
from bseq.sequence import NuclSequence, ProtSequence, CodonSequence


def read_fasta_file(path, seq_type='nucleotide'):
    """Reads the FASTA file and stores the contents
    as a list of Sequence objects

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
                if not seq:
                    seq_string = ''.join(seq)
                    seq_obj = None
                    if seq_type == 'nucleotide':
                        seq_obj = NuclSequence(name, seq_string, description)
                    elif seq_string == 'protein':
                        seq_obj = ProtSequence(name, seq_string, description)
                    elif seq_string == 'codon':
                        seq_obj = CodonSequence(name, seq_string, description)
                    else:
                        raise ValueError('seq_type must be \
                                          "nucleotide", "protein", \
                                          or "codon".')
                    sequence_list.append(seq_obj)
                name, description = line[1:].rstrip().split(' ', 1)
                seq = []
            else:
                seq += [line.rstrip()]

        if not seq:
            seq_string = ''.join(seq)
            seq_obj = None
            if seq_type == 'nucleotide':
                seq_obj = NuclSequence(name, seq_string, description)
            elif seq_string == 'protein':
                seq_obj = ProtSequence(name, seq_string, description)
            elif seq_string == 'codon':
                seq_obj = CodonSequence(name, seq_string, description)
            else:
                raise ValueError('seq_type must be \
                                    "nucleotide", "protein", \
                                    or "codon".')
    return sequence_list
