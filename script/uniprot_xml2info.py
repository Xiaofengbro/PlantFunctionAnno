#!/usr/bin/env python3

import argparse
import gzip
import xml.etree.ElementTree as ET


def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rb")
    return open(filename, "rb")


def local_name(tag):
    """
    去掉 XML namespace
    例如：
    {http://uniprot.org/uniprot}entry
    -> entry
    """
    return tag.split("}")[-1]


def get_text(element):
    if element is None:
        return ""

    return "".join(element.itertext()).strip()


def unique_join(items, sep=";"):
    result = []
    seen = set()

    for x in items:
        x = x.strip()

        if x and x not in seen:
            result.append(x)
            seen.add(x)

    return sep.join(result)


def find_children(element, name):

    return [
        x for x in element.iter()
        if local_name(x.tag) == name
    ]


def parse_uniprot(xml_file, fasta_out, annotation_out):

    count = 0

    with open_file(xml_file) as handle, \
         open(fasta_out, "w") as fa, \
         open(annotation_out, "w") as ann:

        header = [
            "Accession",
            "Entry",
            "Recommended_name",
            "Alternative_name",
            "Gene",
            "Organism",
            "TaxID",
            "EC",
            "GO",
            "KEGG",
            "InterPro",
            "Pfam",
            "PANTHER",
            "SMART",
            "SUPFAM",
            "Gene3D",
            "CDD",
            "Sequence_length"
        ]

        ann.write("\t".join(header) + "\n")

        # 流式解析，适合很大的 UniProt XML
        for event, elem in ET.iterparse(handle, events=("end",)):

            if local_name(elem.tag) != "entry":
                continue

            count += 1

            # ==================================================
            # Accession
            # ==================================================

            accessions = []

            for x in find_children(elem, "accession"):
                if x.text:
                    accessions.append(x.text.strip())

            # 每个 UniProt entry 只使用第一个 accession
            primary_accession = accessions[0] if accessions else ""
            accession = primary_accession

            # ==================================================
            # Entry name
            # ==================================================

            entry_name = ""

            for x in elem:

                if local_name(x.tag) == "name":

                    if x.text:
                        entry_name = x.text

                    break

            # ==================================================
            # Protein name
            # ==================================================

            recommended_names = []
            alternative_names = []

            protein_elements = [
                x for x in elem
                if local_name(x.tag) == "protein"
            ]

            if protein_elements:

                protein_elem = protein_elements[0]

                # 只遍历 protein 的直接子节点
                for name_elem in protein_elem:

                    name_type = local_name(name_elem.tag)

                    if name_type not in ("recommendedName", "alternativeName"):
                        continue

                    # 只读取当前 recommendedName / alternativeName
                    # 的直接 fullName
                    for child in name_elem:

                        if local_name(child.tag) == "fullName":

                            if child.text:
                                full_name = child.text.strip()

                                if name_type == "recommendedName":
                                    recommended_names.append(full_name)

                                elif name_type == "alternativeName":
                                    alternative_names.append(full_name)

            recommended_name = unique_join(recommended_names)
            alternative_name = unique_join(alternative_names)

            # ==================================================
            # Gene
            # ==================================================

            genes = []

            for gene_elem in find_children(elem, "gene"):

                for x in gene_elem:

                    if local_name(x.tag) == "name":

                        if x.text:
                            genes.append(x.text)

            gene = unique_join(genes)

            # ==================================================
            # Organism
            # ==================================================

            organism_name = ""

            organisms = find_children(elem, "organism")

            if organisms:

                names = []

                for x in find_children(organisms[0], "name"):

                    if x.text:
                        names.append(x.text)

                organism_name = unique_join(names)

            # ==================================================
            # TaxID
            # ==================================================

            taxids = []

            if organisms:

                for dbref in find_children(
                    organisms[0],
                    "dbReference"
                ):

                    if dbref.attrib.get("type") == "NCBI Taxonomy":

                        taxid = dbref.attrib.get("id")

                        if taxid:
                            taxids.append(taxid)

            taxid = unique_join(taxids)

            # ==================================================
            # EC number
            # ==================================================

            ecs = []

            for x in find_children(elem, "ecNumber"):

                if x.text:
                    ecs.append(x.text)

            # UniProt XML 有些 EC 在 recommendedName 属性中
            for x in find_children(elem, "recommendedName"):

                ec_attr = x.attrib.get("ecNumber")

                if ec_attr:
                    ecs.append(ec_attr)

            ec = unique_join(ecs)

            # ==================================================
            # Database cross references
            # ==================================================

            go = []
            kegg = []
            interpro = []
            pfam = []
            panther = []
            smart = []
            supfam = []
            gene3d = []
            cdd = []

            for dbref in find_children(elem, "dbReference"):

                db = dbref.attrib.get("type")
                dbid = dbref.attrib.get("id")

                if not dbid:
                    continue

                if db == "GO":
                    go.append(dbid)

                elif db == "KEGG":
                    kegg.append(dbid)

                elif db == "InterPro":
                    interpro.append(dbid)

                elif db == "Pfam":
                    pfam.append(dbid)

                elif db == "PANTHER":
                    panther.append(dbid)

                elif db == "SMART":
                    smart.append(dbid)

                elif db == "SUPFAM":
                    supfam.append(dbid)

                elif db == "Gene3D":
                    gene3d.append(dbid)

                elif db == "CDD":
                    cdd.append(dbid)

            # ==================================================
            # Sequence
            # ==================================================

            sequence = ""
            sequence_length = ""

            sequence_elements = find_children(
                elem,
                "sequence"
            )

            if sequence_elements:

                seq_elem = sequence_elements[0]

                sequence = get_text(seq_elem)

                sequence = (
                    sequence
                    .replace(" ", "")
                    .replace("\n", "")
                    .replace("\r", "")
                    .upper()
                )

                sequence_length = seq_elem.attrib.get(
                    "length",
                    ""
                )

            # ==================================================
            # FASTA
            # ==================================================

            if primary_accession and sequence:
                
                fa.write(f">{primary_accession}\n")
                
                for i in range(0, len(sequence), 60):
                    
                                fa.write(sequence[i:i + 60] + "\n")

            # ==================================================
            # Annotation
            # ==================================================

            row = [
                accession,
                entry_name,
                recommended_name,
                alternative_name,
                gene,
                organism_name,
                taxid,
                ec,
                unique_join(go),
                unique_join(kegg),
                unique_join(interpro),
                unique_join(pfam),
                unique_join(panther),
                unique_join(smart),
                unique_join(supfam),
                unique_join(gene3d),
                unique_join(cdd),
                sequence_length
            ]

            ann.write(
                "\t".join(row) +
                "\n"
            )

            # 清除元素，避免超大 XML 占满内存
            elem.clear()

            # 每 1000 条打印一次
            if count % 1000 == 0:

                print(
                    f"Processed {count} entries...",
                    flush=True
                )

    print()
    print(f"Finished.")
    print(f"Total entries: {count}")
    print(f"FASTA: {fasta_out}")
    print(f"Annotation: {annotation_out}")


def main():

    parser = argparse.ArgumentParser(
        description="UniProt XML to FASTA + annotation"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="UniProt XML or XML.gz"
    )

    parser.add_argument(
        "-o",
        "--fasta",
        default="uniprot.fa",
        help="Output FASTA"
    )

    parser.add_argument(
        "-a",
        "--annotation",
        default="uniprot_annotation.tsv",
        help="Output annotation TSV"
    )

    args = parser.parse_args()

    parse_uniprot(
        args.input,
        args.fasta,
        args.annotation
    )


if __name__ == "__main__":
    main()
