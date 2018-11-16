# -*- coding: utf-8 -*-
"""Data structure models for marker sequences in alignments.
"""
from collections import Iterable
import numpy as np


class Marker(object):
    """Annotates sites in the alignment using a marker sequence
    to label various types of sites.

    Attributes
    ----------
    name : str
        Name of the marker.
    description : str
        Description of the marker.
    char_description : dict
        Dictionary that descripts what each marker character means.
        Keys are allowed characters and values are the corresponding
        description of what a particular marker character means.
    sequence
    encoded_sequence

    """
    def __init__(self, name, char_description: dict, marker_sequence,
                 description=None):
        """Creates a new Marker object.

        Parameters
        ----------
        name : str
            Name of the marker.
        char_desction : dict
            Keys are allowed characters and values are description of what
            the particular marker character means.
        marker_sequence : str
            Marker sequence as a string.
        description : str, optional
            Description for this marker sequence.

        """
        self.name = name
        self.description = description
        self.char_description = char_description
        self._pos_list = ()
        self._char_list = ()
        # Encode sequences into a list of positions and its
        # corresponding characters.
        self._encode(marker_sequence)
        # Check if the characters in the sequence are in the
        # character descriptions.
        self.check_sequence(self._char_list, list(self.char_description.keys()))

    @property
    def sequence(self):
        """Returns the string representation of the marker sequence.
        """
        assert len(self._pos_list) == len(self._char_list)
        len_list = [tup[-1]-tup[0] for tup in self._pos_list]
        return ''.join(
            [char * length for length, char in zip(len_list, self._char_list)]
        )

    @property
    def encoded_sequence(self):
        """Returns the encoded representation of the marker sequence.

        For instance, 0C10N11C20 means from the start to the 10th character is
        C, the 11th character is N, and from the 12th to the last position,
        the marker is C. The encoding 0C10N11C20 is equal to
        CCCCCCCCCCNCCCCCCCCC
        """
        num_list = list(map(str, (tup[0] for tup in self._pos_list)))
        combined_list = [None] * (len(num_list) + len(self._char_list))
        combined_list[::2] = num_list
        combined_list[1::2] = self._char_list
        last = str(self._pos_list[-1][-1])
        return ''.join(list(map(str, combined_list))) + last

    def coords(self, char, inverse=False, explicit=True):
        """Returns the all the positions of the target character, or the all the
        positions where the target character is not found.

        Parameters
        ----------
        char : str
            Marker character of interest
        inverse : bool, optional
            Indicates what set of positions to return with regards to the target
            character given.
            If True, positions where the target character is not found will be
            returned. If False (default), it returns the positions that matched
            the character of interest.
        explicit : bool, optional
            Indicates whether to return a list of numbered positions, or to
            return a list of position intervals.
            If True, returns a complete listing positions,
            ie. [0,1,2,3,4,6,7,8]. If False, returns a list of intervals,
            lower bound inclusive and upper bound exclusive, ie.
            [(0,5), (6,9)]

        Returns
        -------
        list
            List of positions if explicit == True, or a list of intervals if False.
        """
        coords_list = []
        for pos_tup, c in zip(self._pos_list, self._char_list):
            if (char == c and not inverse) or (char != c and inverse):
                if explicit:
                    coords_list += list(range(*pos_tup))
                else:
                    coords_list.append(list(pos_tup))
        return coords_list

    def filter(self, sequence, *exclude_chars, inverse=False,
               output_coords=False):
        """Filters out alignment columns based on the set of marker
        characters to be excluded.

        Parameters
        ----------
        sequence : str
            Sequence to undergo filtering.
        exclude_chars : list of str
            Marker characters to be excluded. The positions of these
            characters will indicate what sites to excluded in the output.
        output_coords : bool, optional
            Outputs the list of positions that passed the filter/s if True,
            otherwise outputs the filtered sequence. By default, `output_coords`
            is False.

        Returns
        -------
        str
            Sequence after removing sites containing the marked characters.
            This means the output sequence must be the same length or shorter
            compared to the input sequence.

        """
        assert len(sequence) == len(self)
        seq_array = np.array(list(sequence))
        coords = []
        for c in exclude_chars:
            coords += self.coords(c, inverse=True if not inverse else False,
                                  explicit=True)
        coords = sorted(set(coords))
        if output_coords:
            return list(coords)
        return ''.join(seq_array[coords])

    def mask(self, sequence, *exclude_chars, mask_char='_'):
        """Masks positions based on the set of marker characters.

        Parameters
        ----------
        sequence : str
            Sequence to be masked.
        exclude_chars : list of str
            Marker characters to be excluded by masking the positions
            where they occur.
        mask_char : str, optional
            Character to replace sites to be excluded.

        Returns
        -------
        str
            Masked sequence that masks positions where the specified characters
            occur. This means the output sequence will be the same length
            as the input sequence.

        """
        assert len(sequence) == len(self)
        seq_array = np.array(list(sequence))
        coords = []
        for c in exclude_chars:
            coords += self.coords(c, inverse=False, explicit=True)
        coords = sorted(set(coords))
        seq_array[coords] = mask_char
        return ''.join(seq_array)

    @staticmethod
    def check_sequence(sequence, allowed_chars):
        """Checks if the sequence only contains the set of allowed characters.

        Paramters
        ---------
        sequence : iterable
            Sequence of characters.
        allowed_chars : iterable
            Set of allowed characters.

        Raises
        ------
        ValueError
            If at least one character in `sequence` is not found
            in `allowed_chars`

        """
        for c in set(sequence):
            if c not in allowed_chars:
                raise ValueError(
                    'character {} is not one of the allowed ' \
                    'characters in the sequence ({})'.format(
                        c, ', '.join(allowed_chars))
                    )

    def _encode(self, sequence):
        """Encodes the marker sequence string into a list of coordinates.

        Parameters
        ----------
        sequence : str
            Marker sequence string to be encoded.

        """
        pos_list = []
        char_list = []
        for i, c in enumerate(sequence):
            if not char_list or char_list[-1] != c:
                pos_list.append(i)
                char_list.append(c)
        pos_list.append(len(sequence))
        self._pos_list = tuple(
            (pos_list[i], pos_list[i+1]) for i in range(len(pos_list)-1)
        )
        self._char_list = tuple(char_list)

    def __repr__(self):
        return self.name + ':' + self.encoded_sequence

    def __len__(self):
        return self._pos_list[-1][-1]

class ConsAlignMarker(Marker):
    """Marker subclass specifically for ConsAlign.

    """
    def __init__(self, marker_sequence,  # pylint: disable=W0102
                 name='ConsAlign_marker_sequence',
                 description=None,
                 char_description={  # pylint doesnt like a dict as a param
                     'C': 'consistent site',
                     'N': 'inconsistent site'
                 }):
        """Creates a new ConsAlignMarker object.

        Parameters
        ----------
        marker_sequence : str
            Marker sequence string generated by the ConsAlign program.
        name : str, optional
            Name of the ConsAlignMarker. By default, the name is
            'ConsAlign_marker_sequence'.
        description : str, optional
            Description for the ConsAlign marker.
        char_description : dict, optional
            Keys are the allowed marker characters and values are the
            description of what the particular marker character means.
            By default, 'C' indicates 'consistent site' while
            'N' indicates 'inconsistent site'

        """
        super().__init__(name, char_description, marker_sequence,
                         description=description)

    def consistent_sites(self, aligned_sequence, marker_char='C'):
        """Returns the aligned sequence containing only the consistent sites
        according to the ConsAlign marker sequence.

        Parameters
        ----------
        aligned_sequence : str
            Aligned sequence
        marker_char : str, optional
            ConsAlign marker character representing a consistently aligned site.
            By default this character is 'C'.

        Returns
        -------
        str
            Aligned sequence that contains only consistent sites. This means
            the output sequence must be equal to or less than the length of the
            input sequence.

        """
        return self.filter(aligned_sequence, marker_char, inverse=True)

    def inconsistent_sites(self, aligned_sequence, marker_char='N'):
        """Returns the aligned sequence containing only the inconsistent sites
        according to the ConsAlign marker sequence.

        Parameters
        ----------
        aligned_sequence : str
            Aligned sequence
        marker_char : str, optional
            ConsAlign marker character representing an inconsistently
            aligned site. By default this character is 'N'.

        Returns
        -------
        str
            Aligned sequence that contains only inconsistent sites. This means
            the output sequence must be equal to or less than the length of the
            input sequence.

        """
        return self.filter(aligned_sequence, marker_char, inverse=True)

    def mask_consistent_sites(self, aligned_sequence,
                              marker_char='C', mask_char='_'):
        """Returns the aligned sequence masking the consistent sites
        according to the ConsAlign marker sequence.

        Parameters
        ----------
        aligned_sequence : str
            Aligned sequence
        marker_char : str, optional
            ConsAlign marker character representing a consistently
            aligned site. By default this character is 'C'.
        mask_char : str, optional
            Character used to mask sites.

        Returns
        -------
        str
            Aligned sequence where consistent sites are masked by the given
            mask character. This means the output sequence must be equal to
            the length of the input sequence.

        """
        return self.mask(aligned_sequence, marker_char, mask_char=mask_char)

    def mask_inconsistent_sites(self, aligned_sequence,
                                marker_char='N', mask_char='_'):
        """Returns the aligned sequence masking the inconsistent sites
        according to the ConsAlign marker sequence.

        Parameters
        ----------
        aligned_sequence : str
            Aligned sequence
        marker_char : str, optional
            ConsAlign marker character representing an inconsistently
            aligned site. By default this character is 'N'.
        mask_char : str, optional
            Character used to mask sites.

        Returns
        -------
        str
            Aligned sequence where inconsistent sites are masked by the given
            mask character. This means the output sequence must be equal to
            the length of the input sequence.

        """
        return self.mask(aligned_sequence, marker_char, mask_char=mask_char)

    def consistent_site_coords(self, marker_char='C'):
        """Returns the positions of consistent sites in the
        ConsAlign marker sequence.

        Parameters
        ----------
        marker_char : str, optional
            ConsAlign marker character representing a consistently
            aligned site. By default this character is 'C'.

        Returns
        -------
        list of int
            Positions of consistent sites according to the
            ConsAlign marker sequence. Positions count from zero.

        """
        return self.coords(marker_char)

    def inconsistent_site_coords(self, marker_char='N'):
        """Returns the positions of inconsistent sites in the
        ConsAlign marker sequence.

        Parameters
        ----------
        marker_char : str, optional
            ConsAlign marker character representing an inconsistently
            aligned site. By default this character is 'N'.

        Returns
        -------
        list of int
            Positions of inconsistent sites according to the
            ConsAlign marker sequence. Positions count from zero.

        """
        return self.coords(marker_char)

class GapMarker(Marker):
    """Marker subclass specifically to mark gaps.
    """
    def __init__(self, marker_sequence,  # pylint: disable=W0102
                 name='Gap_marker_sequence',
                 description=None,
                 char_description={  # pylint doesnt like a dict as a param
                     'O': 'ungapped site',
                     'X': 'site has at least one gap present'
                 }):
        """Creates a new GapMarker object.

        Parameters
        ----------
        marker_sequence : str
            Marker sequence string indicating the presence of gaps
            in the alignment.
        name : str, optional
            Name of the GapMarker. By default, the name is
            'Gap_marker_sequence'.
        description : str, optional
        char_description : dict, optional
            Keys are the allowed marker characters and values are the
            description of what the particular marker character means.
            By default, 'O' indicates 'ungapped site' while
            'X' indicates 'gapped site'

        """
        super().__init__(name, char_description, marker_sequence,
                         description=description)

    def remove_gaps(self, aligned_sequence, marker_char='X'):
        """Returns the aligned sequence containing only the sites
        that do not have gaps based on the gap marker sequence.

        Parameters
        ----------
        aligned_sequence : str
            Aligned sequence
        marker_char : str, optional
            Marker character representing a gapped site.
            By default this character is 'X'.

        Returns
        -------
        str
            Aligned sequence that contains only sites without gaps. This means
            the output sequence must be equal to or less than the length of the
            input sequence.

        """
        return self.filter(aligned_sequence, marker_char)
