#!/usr/bin/env python3

import csv
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--file1", required=True)
parser.add_argument("-a", "--file2", required=True)
parser.add_argument("-o", "--out", default="result.txt")
args = parser.parse_args()

file2_name = os.path.basename(args.file2).rsplit(".", 1)[0]

anno = {}

with open(args.file2, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        accession = row["Accession"].strip()

        organism = row["Organism"].split(";")[0].strip() if row["Organism"] else ""

        items = []

        if row["Recommended_name"].strip():
            items.append(row["Recommended_name"].strip())

        if row["Gene"].strip():
            items.append("GN=" + row["Gene"].strip())

        if organism:
            items.append("OS=" + organism)

        if row["Entry"].strip():
            items.append("EN=" + row["Entry"].strip())

        if row["Sequence_length"].strip():
            items.append("LEN=" + row["Sequence_length"].strip())

        anno[accession] = " ".join(items)

with open(args.file1, encoding="utf-8") as fin, \
     open(args.out, "w", encoding="utf-8") as fout:

    fout.write("gene\t" + file2_name + "\n")

    for line in fin:
        line = line.rstrip("\n\r")

        if not line:
            continue

        parts = line.split()

        gene = parts[0]
        accessions = parts[1].split(",")

        results = []

        for accession in accessions:
            accession = accession.strip()

            if accession in anno:
                results.append(anno[accession])
            else:
                results.append(accession)

        fout.write(gene + "\t" + ";".join(results) + "\n")
