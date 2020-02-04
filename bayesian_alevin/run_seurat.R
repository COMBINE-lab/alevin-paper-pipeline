library(Seurat)
library(tximport)

full_matrix <- tximport("/mnt/scratch5/laraib/alevin2/paperSimAviFinal/em/alevin/quants_mat.gz", type="alevin")
dim(full_matrix$counts)

first_half <- sample(1:4340, 2170, replace=FALSE)
fh.matrix <- full_matrix$counts[, first_half]
dim(fh.matrix)

sh.matrix <- full_matrix$counts[, !(c(1:4340) %in% first_half)]
dim(sh.matrix)

pbmc.fh <- CreateSeuratObject(counts = fh.matrix, project = "fh")
pbmc.sh <- CreateSeuratObject(counts = sh.matrix, project = "sh")

pbmc.fh <- FindVariableFeatures(pbmc.fh)
pbmc.sh <- FindVariableFeatures(pbmc.sh)

transfer.anchors <- FindTransferAnchors(reference = pbmc.fh, 
                                        query = pbmc.sh, 
                                        features = VariableFeatures(object = pbmc.fh), 
                                        reference.assay = "RNA", 
                                        query.assay = "RNA")

anchors <- slot(object = transfer.anchors, name = "anchors")
reference.cells <- slot(object = transfer.anchors, name = "reference.cells")
query.cells <- slot(object = transfer.anchors, name = "query.cells")
anchors <- as.data.frame(x = anchors)

anchors <- anchors[anchors['score'] > 0.5, ]

write.table(file = "./refs.txt", reference.cells[anchors$cell1], quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(file = "./query.txt", query.cells[anchors$cell2], quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(file = "./score.txt", anchors$score, quote = FALSE, row.names = FALSE, col.names = FALSE)
