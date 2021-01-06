# Generation of a Phylogenetic Tree
This code (phylogeny.py) aligns sequences which can then be used to generate a basic phylogenetic tree which is exported it in to PhyloXML format so the tree can be visualised using other phylogenetic software. 

To run:
1. Collate all gene sequences of interest into a single fasta file, for sample data the CRN.fasta file which contains the protein sequences of CRN homologues from various wheat species and Arabidopsis and maize is located in the repository and can be downloaded
2. Create a single working directory with your fasta file containing the protein sequences of interest and the MUSCLE binary which can be downloaded from https://www.drive5.com/muscle/downloads.htm
3. Align the protein sequences you want to investigate:

Make sure you specify the location of the muscle binary and the input sequences, as well as specifying where you want your aligned output to be saved. 
```
# Alignment of sequences
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
from Bio import AlignIO

muscle_exe = "/Users/kj16227/TreeGeneration/muscle3" #specify the location of your muscle exe file
input_sequences = "/Users/kj16227/TreeGeneration/CRN.fasta" #specify the location of the sequences you want to investigate
output_fas = "/Users/kj16227/TreeGeneration/CRN_aligned.fasta" #specify where you want your output file to be saved and its name

def musclealign (Fasta): 
    muscle_cline = MuscleCommandline(muscle_exe, input=Fasta, out=output_fas)
    stdout, stderr = muscle_cline()
    MultipleSeqAlignment = AlignIO.read(output_fas, "fasta") #output file format specified as FASTA
    print(MultipleSeqAlignment) 
#the output provides a visualisation of the alignment so confirms that it has been carried out correctly

musclealign(input_sequences)
``` 
4. Convert aligned sequence output from FASTA to phylip:

This is a necessary step as the aligned output in FASTA format is incompatible with the tree generation software from BioPython
```
#Converting aligned sequence output from FASTA to phylip
from Bio import AlignIO
output_phy= "/Users/kj16227/TreeGeneration/CRN_aligned.phy" #specify the location of where you want the phylip aligned output to be saved
alignments = AlignIO.parse(output_fas, "fasta")
AlignIO.write(alignments, output_phy, "phylip") #specify the format you want the file to be saved to, in this case phylip

aln = AlignIO.read(output_phy, 'phylip')
print(aln) #visualises the output in phylip format to confirm that conversion was successful
```
5. Generate a basic phylogenetic tree:

This constructs a phylogenetic tree and allows visualisation in a basic dendogram format. 
```
#Generation of a basic phylogenetic tree
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
aln = AlignIO.read(output_phy, 'phylip')
constructor = DistanceTreeConstructor()
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print(dm) #prints the distance matrix in plain text
njtree = constructor.nj(dm) #constructs tree using neighbour joining methods
print(njtree) #prints the nj tree in text format
Phylo.draw_ascii(njtree) #print the nj tree in basic dendrogram format
```
6. Generate a tree plot and export it to phyloXML format:

This prints the tree in a plot, however due it is hard to edit the plot using BioPython so this code also enables the export of the tree into PhyloXML format which can then be visualised and edited in other software. 
```
#Generation of a tree plot and export to Phyloxmyl format
Phylo.draw(njtree) #prints the nj tree in a plot
import sys
Phylo.write(njtree, sys.stdout, "phyloxml") #prints the tree in phyloxml format
tree_out= "/Users/kj16227/TreeGeneration/tree_out.xml" #speciify where you want the tree in phyloxml format to be saved
Phylo.write(njtree, tree_out, "phyloxml") #saves the phylogenetic tree in phyloxml output for further use
```
