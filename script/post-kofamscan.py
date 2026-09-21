#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser(description="Parse KofamScan output and add KEGG pathway annotation")

    parser.add_argument("-i", "--kegg", default="kegg.txt", help="KofamScan raw output, default: kegg.txt")
    parser.add_argument("-k", "--ko-pathway", default="ko2pathway.list", help="KO to pathway mapping")
    parser.add_argument("-d", "--pathway", default="pathway_discrpts.txt", help="Pathway description file")
    parser.add_argument("-o", "--kegg-anno", default="kegg.pathway.anno", help="Output annotated file")
    parser.add_argument("-b", "--pathway-bak", default="gene.pathway.txt", help="Output gene-pathway file")
    args = parser.parse_args()
    
    ko2pathway = {}

    with open(args.ko_pathway, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            ko = parts[0]
            pathways = parts[1].split(",")
            ko2pathway[ko] = pathways

    pathway_info = {}

    with open(args.pathway, "r", encoding="utf-8") as f:
        next(f, None)
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            pathway_id = parts[0].strip()
            description = parts[1].strip() if len(parts) > 1 else ""
            pathway_info[pathway_id] = {
                "description": description,
            }

    gene_pathway_records = []
    with open(args.kegg, "r", encoding="utf-8") as fin, \
         open(args.kegg_anno, "w", encoding="utf-8") as fout:
        fout.write("gene\tKEGG\tPathway\n")
        for line in fin:
            if not line.startswith("*"):
                continue
            line = line[1:].strip()
            if not line:
                continue
            parts = line.split(None, 5)
            if len(parts) < 6:
                continue
            gene = parts[0]
            ko = parts[1]
            definition = parts[5]
            kegg_value = ko + "|" + definition
            pathways = ko2pathway.get(ko, [])
            pathway_strings = []
            for pathway_id in pathways:
                pathway_id = pathway_id.strip()
                if not pathway_id:
                    continue
                info = pathway_info.get(
                    pathway_id,
                    {
                        "description": "",
                    }
                )
                description = info["description"]
                if description:
                    pathway_string = (
                        f"{description} ({pathway_id})"
                    )
                else:
                    pathway_string = pathway_id
                pathway_strings.append(pathway_string)
                gene_pathway_records.append([
                    gene,
                    pathway_id,
                    description,
                ])
            pathway_text = "; ".join(pathway_strings) if pathway_strings else "-"
            fout.write(
                gene + "\t" +
                kegg_value + "\t" +
                pathway_text + "\n"
            )
    with open(args.pathway_bak, "w", encoding="utf-8") as fout:
        for record in gene_pathway_records:
            fout.write("\t".join(record) + "\n")

if __name__ == "__main__":
    main()
