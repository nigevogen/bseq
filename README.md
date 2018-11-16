# bseq
Python package for handling biological sequences and alignments

# Quickstart

## Install bseq
```bash
pip install bseq
```

## Nucleotide sequence into FASTA format
```python
from bseq.sequence import NuclSequence
seq = NuclSequence('test_seq_name', 'GATTACA', description='Test nucleotide sequence')
print(seq.fasta_format())
```

## Import alignment from file
```python
from bseq.reader import read_fasta_alignment
aln = read_fasta_alignment('/Users/Me/Downloads/seq.aln.fasta', 
                           seq_type='nucleotide', name='test_alignment', 
                           description='Test alignment for bseq')
```

# Overview
bseq helps biologists import and manipulate biological sequences such as nucleotides,
amino acids, and codons. Unlike other bioinformatics packages in Python, bseq solely focuses
on analyzing and manipulating sequences and alignments.

# Dependencies
bseq uses [numpy](http://www.numpy.org) for many alignment operations, and [nose](https://nose.readthedocs.io/en/latest/) to test its source code. When using `pip` to install `bseq`, pip will automatically install `numpy` and `nose`.
