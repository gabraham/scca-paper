
d <- read.table("relationships_w_pops_121708.txt", header=TRUE, sep="",
   stringsAsFactors=FALSE)
rownames(d) <- d$IID
fam <- read.table("hapmap3_r2_b36_fwd.consensus.qc.poly.fam",
   header=FALSE, sep="", stringsAsFactors=FALSE)
rownames(fam) <- fam$V2

w <- intersect(rownames(d), rownames(fam))
fam <- fam[w, ]
d <- d[w, ]

for(p in unique(d$population)) {
   d2 <- subset(d, population == p, FID:IID)
   write.table(d2, file=paste0("genotyped_samples_", p, ".txt"),
      row.names=FALSE, col.names=FALSE, quote=FALSE)
}

