#!/bin/bash

ls *.bam | paste -s - > header.txt
filelist=$(ls *.bam | paste -s -d" " -)

echo "RepBase..."
/shell/bedtools2/bin/bedtools multicov -D -bams $filelist -bed RepBase_clusters_2kb.bed > RepBase_counts_2kb.txt
echo "RefGene..."
#/shell/bedtools2/bin/bedtools multicov -D -bams $filelist -bed RefGene_clustered.bed > RefSeq_counts.txt

echo "RepBase..."
/shell/bedtools2/bin/bedtools multicov -q 2 -bams $filelist -bed RepBase_clusters_2kb.bed > RepBase_u_counts_2kb.txt
echo "RefGene..."
#/shell/bedtools2/bin/bedtools multicov -q 2 -bams $filelist -bed RefGene_clustered.bed > RefSeq_u_counts.txt
