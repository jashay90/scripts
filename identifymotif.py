# Julie Shay
# Given a fasta sequence alignment file with an unknown stx2 protein sequence and at least one known sequence 
# (the known sequence should have sites of interest in the positions that the script is expecting),
# This script will check sites of interest for differentiating between stx2a/c/d,
# and it will print the list of subtypes that would fit the amino acids at those sites.
# I have been using this for multiple sequence alignments:
# https://github.com/bmh-genomics/Stx_subtype_DBs/blob/main/stx2_aa_subunitA-B_joinedDB_20220505.fasta

from Bio import Seq
from Bio import SeqIO
import argparse
import sys

parser = argparse.ArgumentParser(description="Given a fasta sequence alignment file with an unknown stx2 protein sequence and at least one known sequence (the known sequence should have sites of interest in the positions that the script is expecting), this script will check sites of interest for differentiating between stx2a/c/d, and it will print the list of subtypes that would fit the amino acids at those sites. I have been using this for multiple sequence alignments: https://github.com/bmh-genomics/Stx_subtype_DBs/blob/main/stx2_aa_subunitA-B_joinedDB_20220505.fasta")
parser.add_argument('-i', dest='infile', default='muscleout.fasta', help="""Input fasta alignment file""")
parser.add_argument('-s', dest='seqoi', required=True, help="""title of sequence of interest in fasta file""")
parser.add_argument('-k', dest='knownseq', default='stx2a_subunitA:B_CAA30714.1:CAA30715.1', help="""title of reference sequence in fasta file""")
args = parser.parse_args()

# I think the easiest thing to do is count the gaps in one or more known sequences
# to figure out exactly which sites to check in the unknown sequence.
# String matching could be messed up by gaps (I don't think there would be gaps within this region, but just in case)
def getsitelocs(knownseq, sitesoi): # assumes sites of interest are in ascending order
    actuallocs = {}
    currentloc = 0
    gaps = 0
    for site in sitesoi:
        while site > (currentloc - gaps):
            currentloc = site + gaps
            # grab the subsequence from the beginning to where the first base of interest would be
            gaps = knownseq[0:currentloc].seq.count('-')
        # iteratively grab more until getting to the site of interest
        actuallocs[site] = currentloc
    return actuallocs


def getaa(seqoi, reallocs):
    aa = ""
    info = ""
    aas = []
    for locoi in reallocs.keys():
        realloc = reallocs[locoi]
        aa = seqoi[realloc]
        try:
            info = sitesoi[locoi][aa]
        except KeyError:
            info = []

        aas.append((locoi, aa, info)) # tuple with amino acid and info about subtype it's associated with
    return aas


# infile = "muscleout.fasta"
sitesoi = {312: {"S": ["a", "d"], "F": ["a", "c"]}, 318: {"E": ["a", "d"], "K": ["a", "c"]}, 352: {"E": ["a", "c", "d"]}, 353: {"D": ["a"], "N": ["c", "d"]}, 354: {"D": ["a", "c", "d"]}} # remember 0-based index for keys. Is a dataframe better here?
# knownseq = "stx2a_subunitA:B_CAA30714.1:CAA30715.1"
# seqoi = "stx2c_subunitA:B_AAA16362.1:AAA16363.1"
records = SeqIO.to_dict(SeqIO.parse(args.infile, "fasta"))

# should I do this for a, c, and d in case they disagree?
#try:
reallocs = getsitelocs(records[args.knownseq], sitesoi.keys())
#except KeyError:    
#    print("Couldn't find {}. Quitting now.".format(args.knownseq))
#    exit()


results = getaa(records[args.seqoi], reallocs)
resultsummary = set(results[0][2])
for o in results:
    print("Site {} has amino acid {} which supports subtype {}".format(o[0], o[1], o[2]))
    resultsummary = resultsummary.intersection(o[2])

print("Together these support subtype {}".format(list(resultsummary)))
