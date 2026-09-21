#!/usr/bin/env python3

import csv
import argparse
import re

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True, help="emapper annotations file")
parser.add_argument("-g", "--go", required=True, help="GO reference file")
parser.add_argument("-k", "--kegg", required=True, help="KEGG reference file")
parser.add_argument("-a", "--go-alt", required=True, help="GOid2altid mapping file")
parser.add_argument("-o", "--output", required=True, help="output annotation file")
args = parser.parse_args()

go_dict = {}

with open(args.go, "r", encoding="utf-8-sig") as f:

    reader = csv.DictReader(f, delimiter="\t")

    if "GO_terms" not in reader.fieldnames: raise ValueError("GO reference 中没有 GO_terms 列")
    if "Class" not in reader.fieldnames: raise ValueError("GO reference 中没有 Class 列")
    if "Description" not in reader.fieldnames: raise ValueError("GO reference 中没有 Description 列")
    for row in reader:
        go_id = row["GO_terms"].strip()
        go_class = row["Class"].strip()
        desc = row["Description"].strip()
        if go_id:
            go_dict[go_id] = (go_class, desc)
print(f"GO terms loaded: {len(go_dict)}")

go_alt_dict = {}

with open(args.go_alt, "r", encoding="utf-8-sig") as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        if not row:
            continue
        if row[0].strip() == "GOid2altid":
            continue
        canonical = row[0].strip()
        if not canonical:
            continue
        go_alt_dict[canonical] = canonical
        for alt_id in row[1:]:
            alt_id = alt_id.strip()
            if alt_id:
                go_alt_dict[alt_id] = canonical


print(f"GO alt ID mappings loaded: " f"{len(go_alt_dict)}")

kegg_dict = {}
with open(args.kegg, "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    required_cols = ["Pathway", "Description", "Class 1", "Class2"]

    for col in required_cols:
        if col not in reader.fieldnames: raise ValueError(f"KEGG reference 中没有 {col} 列")
    for row in reader:
        pathway = row["Pathway"].strip()
        description = row["Description"].strip()
        class1 = row["Class 1"].strip()
        class2 = row["Class2"].strip()

        if not pathway:
            continue
        if not pathway.startswith("ko"):
            pathway = "ko" + pathway
        kegg_dict[pathway] = {
            "Description": description,
            "Class 1": class1,
            "Class2": class2
        }

print(f"KEGG terms loaded: " f"{len(kegg_dict)}")

keep_columns = [ 1, 9, 10, 12, 13, 14, 17, 19, 20, 21 ]
keep_index = [ x - 1 for x in keep_columns]
go_rich_file = "GOrich.bak"
kegg_rich_file = "KEGGrich.bak"
go_out = open(go_rich_file, "w", encoding="utf-8")
kegg_out = open(kegg_rich_file, "w", encoding="utf-8")
go_writer = csv.writer(go_out, delimiter="\t", lineterminator="\n")
kegg_writer = csv.writer(kegg_out, delimiter="\t", lineterminator="\n")
with open(
    args.input,
    "r",
    encoding="utf-8-sig",
    newline=""
) as fin, open(
    args.output,
    "w",
    encoding="utf-8",
    newline=""
) as fout:
    reader = csv.reader(fin, delimiter="\t")
    writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

    for row in reader:

        if not row:
            continue
        if row[0].startswith("##"):
            continue
        if row[0].startswith("#"):
            header = [
                x.lstrip("#").strip()
                for x in row
            ]
            new_header = [
                header[i]
                for i in keep_index
            ]
            writer.writerow(new_header)

            if "GOs" not in new_header: raise ValueError("cut 后的文件中没有找到 GOs 列")
            if "KEGG_Pathway" not in new_header: raise ValueError("cut 后的文件中没有找到 KEGG_Pathway 列")

            go_col = new_header.index("GOs")
            kegg_col = new_header.index("KEGG_Pathway")
            print(f"GOs column          : " f"{go_col + 1}")
            print(f"KEGG_Pathway column : " f"{kegg_col + 1}")
            continue
        if len(row) < 21:
            continue
        row[0] = re.sub(r"\.[0-9]+$", "", row[0])
        gene_id = row[0]
        if row[12] != "-":
            pathways = row[12].split(",")
            new_pathways = []
            for pathway in pathways:
                pathway = pathway.strip()
                if not pathway:
                    continue
                if not pathway.startswith("ko"):
                    pathway = "ko" + pathway
                new_pathways.append(pathway)
            if new_pathways:
                row[12] = ",".join(new_pathways)
            else:
                row[12] = "-"
        if row[16] != "-":
            values = row[16].split(",")
            new_values = []
            for value in values:
                value = value.strip()
                if not value:
                    continue
                if not value.startswith("br"):
                    value = "br" + value
                new_values.append(value)
            if new_values:
                row[16] = ",".join(new_values)
            else:
                row[16] = "-"
        new_row = [
            row[i]
            for i in keep_index
        ]
        kegg = new_row[kegg_col]
        if kegg == "-":
            new_row[kegg_col] = "-"
        else:
            pathways = kegg.split(",")
            new_pathways = []
            for pathway in pathways:
                pathway = pathway.strip()
                if not pathway:
                    continue
                if pathway in kegg_dict:
                    description = kegg_dict[
                        pathway
                    ]["Description"]
                    class1 = kegg_dict[
                        pathway
                    ]["Class 1"]
                    class2 = kegg_dict[
                        pathway
                    ]["Class2"]
                    new_pathways.append(f"{description} ({pathway})")
                    kegg_writer.writerow([gene_id, pathway, description, class1, class2])
                else:
                    continue
            if new_pathways:
                new_row[kegg_col] = ",".join(new_pathways)
            else:
                new_row[kegg_col] = "-"
        gos = new_row[go_col]
        if gos == "-":
            new_row[go_col] = "-"
        else:
            go_list = gos.split(",")
            new_gos = []
            for go in go_list:
                go = go.strip()
                if not go:
                    continue
                canonical_go = go_alt_dict.get(go, go)
                if canonical_go in go_dict:
                    go_class = go_dict[canonical_go][0]
                    description = go_dict[canonical_go][1]
                    new_gos.append(f"{description} ({canonical_go})")
                    go_writer.writerow([gene_id,canonical_go,go_class,description])
                else:
                    continue
                if new_gos:
                    new_row[go_col] = ",".join(new_gos)
                else:
                    new_row[go_col] = "-"
                    writer.writerow(new_row)

go_out.close()
kegg_out.close()

print("Done!")
print(f"Output annotation : {args.output}")
print(f"Output GO rich    : {go_rich_file}")
print(f"Output KEGG rich  : {kegg_rich_file}")