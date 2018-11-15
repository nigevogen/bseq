# -*- coding: utf-8 -*-
"""Helper functions for formatting files.
"""
def fasta_formatted_string(name, sequence, description=None, line_width=None):
    """Returns the name and character sequence in the FASTA format.

    Parameters
    ----------
    name : str
        Name describing the sequence. Usually the ID in the sequence in the a FASTA file.
    sequence : str
        Characters describing nucleotide sequences in a FASTA file.
    description : str, optional
        Longer description and other notes about the sequence.
    line_width : int, optional
        The number of characters in a line.

    Returns
    -------
    str
        String in FASTA format consisting of 2 lines, first line is the `name` and
        the second line is `sequence`.

    """
    string = '>' + name
    if description:
        string += ' ' + description
    string += '\n'
    if line_width:
        last = 0
        for i in range(line_width, len(sequence), line_width):
            string += sequence[i-line_width:i] + '\n'
            last = i
        string += sequence[last:]
        return string + '\n'
    string += sequence
    return string + '\n'
