# Julie Shay
# July 22, 2019
# This script makes two output files from a bam alignment file: a file with gene counts/coverages
# and a file with cluster counts/coverages based on the most abundant gene per cluster
import pysam
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', dest='bam', type=str, help="""bam alignment file""")
parser.add_argument('-c', dest='clusters', type=str, required=True, help="""Tab-delimited file with genes and clusters""")
parser.add_argument('-p', dest='minprop', type=float, default=0, help="""minimum proportion of cluster covered""")
parser.add_argument('-g', dest='geneout', type=str, default="geneout.tsv", help="""output file for gene counts""")
parser.add_argument('-o', dest='clusterout', type=str, default="clusterout.tsv", help="""output file for cluster counts""")
args = parser.parse_args()

def gene_count(gene, bam):
    read_count = bamfile.count(contig=gene, read_callback='all')
    cov_counts = bamfile.count_coverage(contig=gene)
    cov_number = 0
    total_number = 0
    for a, c, g, t in zip(*cov_counts):   # doesn't matter whether the order is right
        total_number = total_number + 1
        if (((a + g) + c) + t) != 0:
            cov_number = cov_number + 1
    prop_covered = float(cov_number) / float(total_number)
    return read_count, prop_covered


clusterdf = pd.read_csv(args.clusters, sep="\t")
genelist = list(clusterdf["ID"])
bamfile = pysam.AlignmentFile(args.bam, 'rb')
# for every gene in genelist, get all the reads that mapped to that gene from bamfile
counts = []
coverages = []
for g in genelist:
    count, coverage = gene_count(g, bamfile)
    counts.append(count)
    coverages.append(coverage)

genedf = pd.DataFrame({"Gene": genelist, "Hits": counts, "Gene Fraction": coverages})
genedf.set_index("Gene").to_csv(args.geneout, sep="\t")
clusterdf = pd.merge(clusterdf, genedf.rename(columns={"Gene": "ID"}), on="ID").drop("ID", axis=1).groupby("Cluster").aggregate(np.max)
clusterdf[clusterdf["Gene Fraction"] >= args.minprop].to_csv(args.clusterout, sep="\t")
