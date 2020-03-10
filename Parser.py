# Robert Bennett - rbennet8@uncc.edu - 2/12/2020

# Sequence - Should hold sequence of each object and a method to print them
# DNASequence - Holds methods to transcribe (moved from the Sequence class) and translate
    # Can only contain the characters A, T, C, and G
# ProteinSequence - Holds method to display structure
    # Can only contain the characters M, F, L, C, Y, W, P, H, Q, R, I, T, N, K, S, V, A, D, E, G, AND *


#ENTER PATHS TO FASTA SEQUENCES AT THE LINES CONTAINING ################################################################



# Parent class to DNA and Protein sequences
class Sequence:
    def __init__(self, seq):
        self.seq = seq

    def __repr__(self):
        return self.seq



# Child class of Sequence and stores the nucleotide sequence
# Also has methods to transcribe and translate the sequence to a protein
class DNASequence(Sequence):
    def __init__(self, seq):
        # If secondary sequence check passes, then the constructor calls the super class to create the object
        if self.seqCheck(seq):
            super().__init__(seq)
        # Else, a boolean is returned, which prevents the program from continuing
        else:
            return False

    # Second check to make sure whatever is being passed is a DNA sequence, in case it isn't submitted via parse method
    def seqCheck(self, seq):
        # Checks all characters in seq to make sure they match all characters in ATCG and returns boolean
        bool = all(x in seq for x in "ATCG")
        return bool

    # Method that replaces every T in the string with a U and returns string
    # Method from previous lab
    def transcribe(self):
        mrna = self.seq.replace('T', 'U')
        return mrna

    # Method that takes transcribed sequence and translates it, then saves the sequence as a protein object
    def translate(self):
        aa_dict = {'M': ['ATG'],
                   'F': ['TTT', 'TTC'],
                   'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                   'C': ['TGT', 'TGC'],
                   'Y': ['TAC', 'TAT'],
                   'W': ['TGG'],
                   'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                   'H': ['CAT', 'CAC'],
                   'Q': ['CAA', 'CAG'],
                   'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                   'I': ['ATT', 'ATC', 'ATA'],
                   'T': ['ACT', 'ACC', 'ACA', 'ACG'],
                   'N': ['AAT', 'AAC'],
                   'K': ['AAA', 'AAG'],
                   'S': ['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'],
                   'V': ['GTT', 'GTC', 'GTA', 'GTG'],
                   'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                   'D': ['GAT', 'GAC'],
                   'E': ['GAA', 'GAG'],
                   'G': ['GGT', 'GGC', 'GGA', 'GGG'],
                   '*': ['TAA', 'TAG', 'TGA']}

        # Declaring variables used for while loop, triplet selection, and protein sequence
        x = 0
        y = 3
        protein = ""
        # While the last triplet position is less than the length of the sequence, continue; prevents while loop from searching
        # for bases in an out of bounds index
        while y < len(self.seq):
            # Setting the triplet at the beginning of each loop
            triplet = self.seq[x:y]
            # For the key and value in amino acid disctionary, if triplet matches the value, append key to protein sequence
            for key, val in aa_dict.items():
                if triplet in val:
                    protein += key
            # Incrementing for triplets
            x += 3
            y += 3
        # Returns protein object
        return ProteinSequence(protein)



# Child class of Sequence and stores the protein translation of the sequence
class ProteinSequence(Sequence):
    def __init__(self, seq):
        # If secondary sequence check passes, then the constructor calls the super class to create the object
        if self.seqCheck(seq) == True:
            super().__init__(seq)
        # Else, a boolean is returned, which prevents the program from continuing
        else:
            return False

    # Second check to make sure whatever is being passed is a DNA sequence, in case it isn't submitted via parse method
    def seqCheck(self, seq):
        # Checks all characters in seq to make sure they match all characters in MFLCYWPHQRITNKSVADEG* and returns boolean
        bool = all(x in seq for x in "MFLCYWPHQRITNKSVADEG")
        return bool

    # This method would search protein databases for similar sequences and return the structure of the closest match
    def displayStructure(self):
        pass


# Takes in a "label", which is the name of the sequence, and a sequence object, where the sequence is stored.
class SequenceRecord:
    def __init__(self, label, seqObj):
        # Checking to make sure a Sequence object is being passed, by making sure it is an instance of the parent class Sequence,
        # then storing the value
        if isinstance(seqObj, Sequence):
            self.seqObj = seqObj
            self.label = label
        else:
            return False

    # Method to return an output for the class
    def __repr__(self):
        return self.label + "\n" + self.seqObj.seq



# Function that takes in a FASTA file and separates information into separate variables
def parse(path):
    label = None
    seq = ""
    file = open(path)
    # For loop that increments through each line of the FASTA file
    for line in file:
        # Handles line that begins with >
        if line.startswith(">"):
            if label:
                # Checking sequence and handling it accordingly; if 1, create DNA object; if 2, create Protein object
                if checkSeq(seq) == 1:
                    seqRecord = SequenceRecord(label, DNASequence(seq))
                    yield seqRecord
                elif checkSeq(seq) == 2:
                    seqRecord = SequenceRecord(label, ProteinSequence(seq))
                    yield seqRecord
                label = None
                seq = ""
            label = line.rstrip().lstrip(">")
        # If no > then concats line to sequence
        else:
            seq += line.rstrip()
    # Checking sequence and handling it accordingly; if 1, create DNA object; if 2, create Protein object
    if checkSeq(seq) == 1:
        seqRecord = SequenceRecord(label, DNASequence(seq))
        yield (seqRecord)
    elif checkSeq(seq) == 2:
        seqRecord = SequenceRecord(label, ProteinSequence(seq))
    yield seqRecord

# Checks sequence to see if it's protein or DNA, then returns an int or boolean
def checkSeq(sequence):
    # if sequence is DNA, return 1
    if all(x in sequence for x in "ATCG"):
        return 1
    # If sequence is Protein, return 2
    elif all(x in sequence for x in "MFLCYWPHQRITNKSVADEG"):
        return 2
    # If sequence is neither, return boolean and break program to prevent it from conitnuing
    else:
        return False

# Takes in path for FASTA file and passes it to parse method
path = r'DNA FASTA PATH' ###############################################################################################
tert = []
for seq in parse(path):
    tert.append(seq)
# Prints out information before and after calling transcription method
print("Printing name and sequence of DNA FASTA:\n", tert[0], "\n")

# Calling transcribe method and printing the result
rna = tert[0].seqObj.transcribe()
print("Printing mRNA sequence:\n", rna, "\n")

# Prints just the sequence of the obj in the sequence record
print("Printing just the sequence from the Parent Sequence class:\n", tert[0].seqObj, "\n")

# Creating protein variable to hold the object being returned from translate method in DNASequence class and printing it
protein = tert[0].seqObj.translate()
print("Printing the translation of DNA sequence:\n", protein, "\n")

# Taking in another path, this time to a protein FASTA, and outputting the information; making sure the program can handle
# both types of FASTA files
path2 = r'PROTEIN FASTA PATH' ##########################################################################################
tertP = []
for seq in parse(path2):
    tertP.append(seq)
print("Printing name and sequence of Protein FASTA:\n", tertP[0])