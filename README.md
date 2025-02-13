# sequence-alignment

A python implementation of the Needleman-Wunsch algorithm for global sequence alignment, as described in _Biological sequence analysis: Probabilistic models of proteins and nucleic acids_ by Durbin et al.

## Dependencies

This script requires **numpy** and **pandas** in order to run. They can be installed using:
```
$ pip install -r requirements.txt
```
## Usage

To run the script, use:
```
$ python src/alignment.py input_fasta substitution_matrix [-o output]
```

```
positional arguments:
    input_fasta           The path of the file containing the two sequences to be aligned in FASTA format.
    substitution_matrix   The path of the subsitution matrix file with space-separated values. The order of the rows and columns is: "ARNDCQEGHILKMFPSTWYV".

options:
    -o output, --output output
                        The path of the alignment file to be outputted. If not set, will simply print the results
```
