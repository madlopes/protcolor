#!/usr/bin/python

#############################################################################
#
#   This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	Author: Tiago J.S. Lopes - University of Wisconsin, Madison, USA
#	Contact: dasilvalopes@wisc.edu 
#
#############################################################################

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from pymol import cmd
from Bio.PDB import PDBList
from re import sub
import sys

###
## Arguments:
## Structure file to open, e.g. 1RUZ.pdb
## Prefix to use, e.g. H1HA
## Conserved FASTA file: the file with the sequence conservation information
###

## Function to calculate the similarity
## of two sequences in an alignment

def get_similarity(align):
	
	conserved_seq = str(list(align)[0].seq)
	chain_seq = str(list(align)[1].seq)
		
	count = -1
	
	identical_aa = 0;
	different_aa = 0;
	
	for i in range(0, len(chain_seq)):
		if chain_seq[i] != '-':
			count = count + 1
	
			if chain_seq[i] == conserved_seq[i]:
				identical_aa = identical_aa + 1
			else:
				different_aa = different_aa + 1
	
	similarity = (identical_aa * 100) / count

	return similarity;



#############################################################################
##
## Main program
##
#############################################################################

pdb_file = sys.argv[2]
prefix_outputfile = sys.argv[3]
conserved_fasta_file = sys.argv[4]

print("%s %s %s" % (sys.argv[2], sys.argv[3], sys.argv[4]))

## Prepare the structures

## First for the sequence analysis

p = PDBParser(PERMISSIVE=1)
structure = p.get_structure(pdb_file, pdb_file)

## Next the structure for colouring

cmd.load(pdb_file)
cmd.show_as("surface")
cmd.color("white")

for model in structure:
	for chain in model:
		print("Working on chain %s." % (chain.get_id()))
		seq = list()
		position = list()
		
		for residue in chain:
			res_id = residue.get_id()
			if res_id[0] == ' ':
				seq.append(three_to_one(residue.get_resname()))
				position.append(res_id[1])
				#print("%s -- %d" % (three_to_one(residue.get_resname()), res_id[1]))


		my_prot = Seq(str(''.join(seq)), IUPAC.protein)

		tmp_seq1 = SeqIO.read(conserved_fasta_file, "fasta", IUPAC.IUPACProtein())
		
		my_x_seq = Seq(sub('[a-z]', 'X', str(tmp_seq1.seq)))

		seq1 = SeqRecord(my_x_seq, id=tmp_seq1.id, name="", description="")
		
		#print("%s" % seq1)
		
		seq2 = SeqRecord(my_prot, id="Chain_" + chain.get_id(), name="", description="")
		
		myseqs = [seq1, seq2]
		
		fasta_filename = prefix + "_" + chain.get_id() + ".fasta"
		align_filename = prefix + "_" + chain.get_id() + ".align"

		SeqIO.write(myseqs, fasta_filename, "fasta")

		cline = MuscleCommandline(input=fasta_filename, out=align_filename, clw="1", clwstrict="1")

		cline()

		aln = AlignIO.read(align_filename, "clustal")

		seq_similarity = get_similarity(aln)
		
		print("\n\t ** Chain %s --- Similarity: %3.2f%%" % (chain.get_id(), seq_similarity))
		
		if seq_similarity < 30:
			print("\t ** Chain ignored\n")
			continue
		else:
			print("\t ** Minimum similarity achieved. Using it further\n")
			
			conserved_seq = str(list(aln)[0].seq)
			chain_seq = str(list(aln)[1].seq)
		
			count = -1
			residues_to_colour = list()
			
			for i in range(0, len(chain_seq)):
				if chain_seq[i] != '-':
					count = count + 1


				# This option will colour in red only the conserved residues
				this_residue = conserved_seq[i]
				if chain_seq[i] != '-' and conserved_seq[i] == 'X':
				
					output_filename = prefix_outputfile + ".pse"
					
					selection_string = "resi " + format(position[count]) + " and chain " + chain.get_id()
					cmd.select("residuedata", selection_string)
					cmd.color ('red', "residuedata")
					cmd.delete("residuedata")
					cmd.save(output_filename, format="pse")
					


cmd.quit()			
				
				
				
								
				
				
				
				
