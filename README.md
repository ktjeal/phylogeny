# Phylogeny
This code aligns protein sequences which can then be used to generate a basic phylogenetic tree which is exported it in to phyloXML format so the tree can be visualised using other phylogenetic software. 

Example data is provided and contains the protein sequences of CRN homologues from various wheat species and Arabidopsis and maize which act as an outgroup. 

To run:
1. Create a single directory with your protein sequences of interest and the MUSCLE binary which can be downloaded from (https://www.drive5.com/muscle/downloads.htm)
2. For sample data download FASTA file into this directory
3. Align the protein sequences you want to investigate:
```
# code block
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

