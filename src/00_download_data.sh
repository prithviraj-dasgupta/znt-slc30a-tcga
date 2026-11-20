#!/bin/bash
# Download TCGA PanCanAtlas expression + survival tables from the UCSC Xena hub
set -e
mkdir -p data/raw && cd data/raw
curl -L -o expression.xena.gz "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
curl -L -o survival.tsv "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp"
