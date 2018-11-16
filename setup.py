import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bseq",
    version="1.0.1",
    author="Kent Kawashima",
    author_email="kentkawashima@nig.ac.jp",
    description="Python package for handling biological sequences and alignments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/nigevogen/bseq",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    install_requires=[
        'numpy',
        'nose',
    ],
)
