#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
library(knitr)
library(data.table)
library(ggplot2)
library(GenomicRanges)
# local vs UGER
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir <- "/xchip/beroukhimlab/"
}
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
chromothripsis_overlap <- fread(paste0(workdir, "youyun/complexSV/data/validation/PCAWG_chromothripsisOverlap.txt"))
chromothripsis_overlap_pass_all = chromothripsis_overlap[(Pass_ALL)]
ct_pass_all_granges = GRanges(
    seqnames = gsub('23','X',gsub('^chr|p$|q$', '', chromothripsis_overlap_pass_all$chr)),
    ranges = IRanges(start = chromothripsis_overlap_pass_all$start, end = chromothripsis_overlap_pass_all$end),
    Sample = chromothripsis_overlap_pass_all$samplename
)
#
#
#
#
#
#
#
cluster_bedpe <- fread(paste0(workdir, "youyun/complexSV/data/TCGA/SV/total_bedpe/cluster_bedpe_concat.csv"))
cluster_ranges = rbind(melt(
    cluster_bedpe[,.(seqnames = seqnames_breakend1, start_breakend1, Sample, cluster_ID = paste0(Sample,'_',cluster_ID), SV_ID)],
    id.vars = c('Sample',"cluster_ID", "SV_ID",'seqnames'), variable.name = "breakend", value.name = "value"
),melt(
    cluster_bedpe[,.(seqnames = seqnames_breakend2, start_breakend2, Sample, cluster_ID = paste0(Sample,'_',cluster_ID), SV_ID)],
    id.vars = c('Sample',"cluster_ID", "SV_ID",'seqnames'), variable.name = "breakend", value.name = "value"
))[
    ,.(min = min(value), max = max(value)), by = c('Sample',"cluster_ID", "seqnames")
]
cluster_granges = GRanges(
    seqnames = cluster_ranges$seqnames,
    ranges = IRanges(start = cluster_ranges$min, end = cluster_ranges$max),
    cluster_ID = cluster_ranges$cluster_ID,
    Sample = cluster_ranges$Sample
)
#
#
#
#
#
graph_embedding <- fread(paste0(
    workdir, "youyun/complexSV/data/TCGA/graph_embedding/UGformerV2_UnSup_lr0.0005_bs16_ep30_ff1024_nn5_sd512_do0.5_nl2_nt2_graph_embeddings.csv"
))
dim(graph_embedding)
rsums = apply(graph_embedding, 1, function(x) sum(as.numeric(x[2:89])))
summary(rsums)
#
#
#
#
#
#
#
#| results: hold
cluster_overlaps_w_ct = cluster_ranges[
    ,.(ct_overlap = length(findOverlaps(
        cluster_granges[cluster_granges$cluster_ID == cluster_ID], 
        ct_pass_all_granges[ct_pass_all_granges$Sample == Sample]
    )) > 0),
    .(cluster_ID, Sample)
]
print('Breakdown of clusters that overlap with chromothripsis events:')
print(table(cluster_overlaps_w_ct$ct_overlap))
write.table(
    cluster_overlaps_w_ct,
    paste0(workdir, "youyun/complexSV/data/validation/cluster_overlaps_w_ct.csv"),
    sep = ",", quote = F, row.names = F

)
#
#
#
#
#
#
