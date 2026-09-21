#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Extract GO annotations from InterProScan and generate GO annotation files.")
parser.add_argument("-i", "--input", required=True, help="InterProScan result file")
parser.add_argument("-a", "--altid", required=True, help="GOid2altid.txt")
parser.add_argument("-d", "--description", required=True, help="GO_descript.txt")
parser.add_argument("-o", "--go-anno", default="go.anno", help="Output corrected GO annotation file")
parser.add_argument("-b", "--go-bak", default="gene.GOrich.bak", help="Output one-to-one GO annotation file")

args = parser.parse_args()

alt2go = {}
with open(args.altid, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        real_go = fields[0].strip()
        alt2go[real_go] = real_go
        for alt_id in fields[1:]:
            alt_id = alt_id.strip()
            if alt_id:
                alt2go[alt_id] = real_go

go_info = {}
with open(args.description, "r") as f:
    f.readline()
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        go_id = fields[0].strip()
        go_class = fields[1].strip()
        description = fields[2].strip()
        go_info[go_id] = (go_class, description)

gene_go = {}
gene_order = []
seen_gene_go = set()

with open(args.input, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 14:
            continue
        gene = fields[0].strip()
        go_field = fields[13].strip()
        if "GO" not in go_field:
            continue
        go_field = go_field.replace("(InterPro)", "")
        go_field = go_field.replace("(PANTHER)", "")
        go_field = go_field.replace("|", ";")

        for go_id in go_field.split(";"):
            go_id = go_id.strip()
            if not go_id:
                continue
            if not go_id.startswith("GO:"):
                continue
            real_go = alt2go.get(go_id, go_id)
            if real_go not in go_info:
                continue
            key = (gene, real_go)
            if key in seen_gene_go:
                continue
            seen_gene_go.add(key)
            if gene not in gene_go:
                gene_go[gene] = []
                gene_order.append(gene)
            gene_go[gene].append(real_go)

with open(args.go_anno, "w") as fout:
    fout.write("gene\tGO_trems\n")
    for gene in gene_order:
        annotations = []
        for go_id in gene_go[gene]:
            if go_id in go_info:
                description = go_info[go_id][1]
                annotations.append(
                    f"{description} ({go_id})"
                )
            else:
                annotations.append(go_id)
        fout.write(gene + "\t" + "; ".join(annotations) + "\n")

with open(args.go_bak, "w") as fout:
    for gene in gene_order:
        for go_id in gene_go[gene]:
            if go_id not in go_info:
                continue
            go_class, description = go_info[go_id]
            fout.write(f"{gene}\t{go_id}\t{go_class}\t{description}\n")
