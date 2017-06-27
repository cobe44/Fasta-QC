import os
import sys
import math
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna
import argparse

parser=argparse.ArgumentParser(description='Create a concatenated alignment for phylogenomic analysis from subpar baits capture data.')

parser.add_argument('fastafile', metavar='fasta infile,', help="Path to fasta file")
parser.add_argument('outfile', metavar='Outfile name,', help="Name of output file")
parser.add_argument('exclude_num', metavar='num gaps,', default=0, help='Number of sequences in alignment allowed to have a gap at a given position before that position is eliminated from the alignment.')
parser.add_argument('already_cleaned', metavar='Already cleaned,', nargs='?', default=None, help='So you dont have to do that malarkey again')

args=parser.parse_args()
fastapath = str(args.fastafile)
outfile = str(args.outfile)
exclude = int(args.exclude_num)
if(args.already_cleaned != None):
    already_cleaned = True
else:
    already_cleaned = None


def fastaparse(fasta):
    #It parses a fasta. Very exciting stuff.
    seq_store=[]
    names=[]
    for record in SeqIO.parse(fasta, "fasta"):
        names.append(record.id)
        seq_store.append(record.seq)

    return seq_store, names

def removal_team(sequences, i):
    #Remove a position from MSA given sequences and position index
    new_sequences = []
    if(i == 0):
        for sequence in sequences:
            new_sequences.append(sequence[1:])
    if(i == len(sequences[0])-1):
        for sequence in sequences:
            new_sequences.append(sequence[:(i)])
    else:
        for sequence in sequences:
            new_sequences.append(sequence[:(i)] + sequence[(i+1):])
    return new_sequences

def allgap_detector(sequences, i):
    for sequence in sequences:
        if(sequence[i] != "-"):
            return False
    return True

def detector_elim(sequences, detector):
    for i in range(len(sequences[0])-1, 0, -1):
        if(detector[i] == 1):
            sequences = removal_team(sequences, i)
    return sequences

def initial_cleanup(sequences):
    #Remove all positions in MSA comprised of only gap characters
    detector = np.zeros(len(sequences[0]))
    for i in range(len(sequences[0])):
        if(allgap_detector(sequences, i) == True):
            detector[i] = 1
        elif(allgap_detector(sequences, i) == False):
            pass
    sequences = detector_elim(sequences, detector)
    print("Initial cleanup complete")
    return sequences

def gap_detector(sequences, i):
    for sequence in sequences:
        if(sequence[i] == "-"):
            return False
    return True

def threshold_detector(sequences, i, threshold):
    count = 0
    for sequence in sequences:
        if(sequence[i] == "-"):
            count += 1
    if(count >= threshold):
        return False
    else:
        return True

def take_out_trash(sequences, exclude):
    #Remove positions in MSA based upon number of gap characters included; wobble bases are OK
    if(exclude == 0):
        #For if you don't want any gap characters because they clash with your personal aesthetic
        detector = np.zeros(len(sequences[0]))
        for i in range(len(sequences[0])):
            if(gap_detector(sequences, i) == False):
                detector[i] = 1
    else:
        #If you want to keep some more biologically relevant data in your alignment even though it'll make your bootstrap values go down
        detector = np.zeros(len(sequences[0]))
        for i in range(len(sequences[0])):
            if(threshold_detector(sequences, i, exclude) == False):
                detector[i] = 1
    sequences = detector_elim(sequences, detector)

    return sequences

def writer(names, seq_store):
    proper_seqs = []
    for i in range(len(names)):
        proper_seqs.append(SeqRecord(seq_store[i], id = names[i], description=""))

    SeqIO.write(proper_seqs, outfile, "fasta")
    print("Data successfully written")

def main():
    sequences, names = fastaparse(fastapath) #Get that fasta file!

    if(already_cleaned == None):
        #Well you'll have to clean it up in that case, won't you?
        cleaned_sequences = initial_cleanup(sequences)
        #Now take out the trash like your mother told you to hours ago, young man
        cleaned_sequences = take_out_trash(cleaned_sequences, exclude)
        #Now write the resulting cleaned sequences to output fasta
        writer(names, cleaned_sequences)
    else:
        #For if you've already cleaned your room but forgot to take out the trash
        cleaned_sequences = take_out_trash(sequences, exclude)
        writer(names, cleaned_sequences)


if __name__ == "__main__":
    main()
