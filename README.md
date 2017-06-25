# Protcolor
Script to color protein structure residues, using as input protein sequences.

Useful, yet simple tool to color-code and export a protein structure file, according to selected amino acids in a FASTA sequence.

The input is:

- A file in FASTA format, with the amino acids which should be colored in lower-case letters, and the remaining in upper-case letters.
- A prefix to rename the output file.
- A PDB structure file.

The program proceeds by aligning the FASTA sequence to the sequence from the protein structure, and color coding the selected residues.

The command line to start the program is:

pymol protcolor.py --- H3HA_2YPG.pdb HK68 hk68.fasta

To run the program you will need to have installed:

- Biopython
- Muscle
- Pymol

Please have a look at the source code for the library dependencies, and contact me if you have doubts or suggestions for improvements.



