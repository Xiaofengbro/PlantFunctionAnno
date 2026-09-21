# 🌱 PlantFunctionAnno
A pipeline for plant genome functional annotation, including plant-specific database construction, annotation transfer, and multi-database result integration.

# ⚠️  Problem Encountered for Me
For functional annotation of annotated genes, the eggNOG Web tool is commonly used. However, when applied to plant genomes, it often generates a large number of annotations that are not plant-specific, particularly in the GO and Pathway categories. These irrelevant annotations can affect downstream functional enrichment analyses, such as those performed for transcriptomic datasets, and may lead to biologically misleading results.

# 🔧 My Approach
This repository provides commands for constructing plant-specific databases from major annotation resources, including NR, eggNOG, UniProt, and KEGG.

In addition, a Bash script for annotation information transfer is provided to help filter and avoid non-plant-specific annotations, thereby improving the reliability of downstream functional annotation and enrichment analyses.

# 🤝 Note❗❗❗❗❗
This workflow  does not provide any form of warranty.
