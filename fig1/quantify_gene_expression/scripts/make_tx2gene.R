# ......................................................
# Make tx2gene file
# 10/05/22
# Hugues HERRMANN
# ......................................................
library(tximportData)
library(GenomicFeatures)


# ......................................................
#
#   PARSE ARGUMENT ----
#
# ......................................................
gff <- "/mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/annotation_genome_ref/gencode_hg38v42_annotation.gff3.gz"
out <- "/mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/annotation_genome_ref/gencode_hg38v42_annotation_tx2gene.csv"


# ......................................................
#
#   MAKE TX2G ----
#
# ......................................................
txdb <- makeTxDbFromGFF(file = gff)

keys <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys, "GENEID", "TXNAME")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(tx2gene, out, sep = ",", row.names = FALSE)
