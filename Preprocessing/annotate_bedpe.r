suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(GenomicRanges))

# local vs UGER
if (grepl('^/Users',Sys.getenv("HOME"),ignore.case = TRUE)) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir <- "/xchip/beroukhimlab/"
}

# Helper functions -----------------------------------------------------------------------------------------------------------

annotate_around_coordinate = function(
    chr, start, replication_timing, gc_content, gene_ranges, 
    repetitive_elements,telo_granges, cen_granges, TADs, CNAs
){
    # annotate the coordinate -----------------------------------------------------------------------------------------------------------
    # replication timing -- just get the value of the range that the coordinate belongs to
    rep_time = replication_timing[queryHits(findOverlaps(replication_timing, GRanges(chr, IRanges(start, start))))]$score

    # GC content -- proportion per kb 
    gc_value = gc_content[queryHits(findOverlaps(gc_content, GRanges(chr, IRanges(start, start))))]$gc

    # Gene Density -- density within a 1mb window
    # Density in 1Mb windows, sliding every 1kb to centre on the pixel 
    gene_density = length(findOverlaps(gene_ranges, GRanges(chr, IRanges(start-10^6/2, start+10^6/2))))

    # repetitive elements -- log10 (kb distance +1)
    dist_to_line = log10(1 + distanceToNearest(
        GRanges(chr, IRanges(start, start)),
        repetitive_elements[repetitive_elements$repeat.element == "LINE"]
    )@elementMetadata$distance/1000)
    dist_to_sine = log10(1 + distanceToNearest(
        GRanges(chr, IRanges(start, start)),
        repetitive_elements[repetitive_elements$repeat.element == "SINE"]
    )@elementMetadata$distance/1000)
    dist_to_ltr = log10(1 + distanceToNearest(
        GRanges(chr, IRanges(start, start)),
        repetitive_elements[repetitive_elements$repeat.element == "LTR"]
    )@elementMetadata$distance/1000)
    dist_to_dna = log10(1 + distanceToNearest(
        GRanges(chr, IRanges(start, start)),
        repetitive_elements[repetitive_elements$repeat.element == "DNA"]
    )@elementMetadata$distance/1000)
    dist_to_sr = log10(1 + distanceToNearest(
        GRanges(chr, IRanges(start, start)),
        repetitive_elements[repetitive_elements$repeat.element == "Simple_repeat"]
    )@elementMetadata$distance/1000)

    # TAD -- log10 (kb distance +1)
    dist_to_TAD = log10(1 + distanceToNearest(
        GRanges(chr, IRanges(start, start)),
        TADs
    )@elementMetadata$distance/1000)

    # telomere and centromere -- distance to the closest telomere or centromere log10 (mb distance + 1)
    dist_to_telomere = log10(1 + distanceToNearest(
        GRanges(chr, IRanges(start, start)),
        telo_granges
    )@elementMetadata$distance/10^6)
    dist_to_centromere = log10(1 + distanceToNearest(
        GRanges(chr, IRanges(start, start)),
        cen_granges
    )@elementMetadata$distance/10^6)

    # CNAs -- return the total_cn, major_cn, minor_cn of where the coordinate is
    total_cn = CNAs[queryHits(findOverlaps(CNAs, GRanges(chr, IRanges(start, start))))]$total_cn
    major_cn = CNAs[queryHits(findOverlaps(CNAs, GRanges(chr, IRanges(start, start))))]$major_cn
    minor_cn = CNAs[queryHits(findOverlaps(CNAs, GRanges(chr, IRanges(start, start))))]$minor_cn

    return(c(
        ifelse(length(rep_time) == 0, NA, rep_time),
        ifelse(length(gc_value) == 0, NA, gc_value),
        gene_density,
        ifelse(length(dist_to_line) == 0, NA, dist_to_line),
        ifelse(length(dist_to_sine) == 0, NA, dist_to_sine),
        ifelse(length(dist_to_ltr) == 0, NA, dist_to_ltr),
        ifelse(length(dist_to_dna) == 0, NA, dist_to_dna),
        ifelse(length(dist_to_sr) == 0, NA, dist_to_sr),
        ifelse(length(dist_to_TAD) == 0, NA, dist_to_TAD),
        ifelse(length(dist_to_telomere) == 0, NA, dist_to_telomere),
        ifelse(length(dist_to_centromere) == 0, NA, dist_to_centromere),
        ifelse(length(total_cn) == 0, NA, total_cn),
        ifelse(length(major_cn) == 0, NA, major_cn),
        ifelse(length(minor_cn) == 0, NA, minor_cn)
    ))
}


if (!interactive()) {
    option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = "",
        help = "Input SvABA BEDPE file from svaba_vcf2bedpe.r", metavar = "input"
    ),
    make_option(c("-c", "--cluster"),
        type = "character", default = "",
        help = "Input BEDPE file from run_cluster_sv.r", metavar = "clusters"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = "",
        help = "Output BEDPE file path to write to", metavar = "output"
    )
    )

    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    input_path <- opt$input
    cluster_path <- opt$cluster
    output_path <- opt$output
    print(paste0("Input BEDPE file: ", input_path))
    print(paste0("Input cluster file: ", cluster_path))

    # load the BEDPE files -----------------------------------------------------------------------------------------------------------
    # input_path = '~/Documents/HMS/PhD/beroukhimlab/broad_mount/youyun/complexSV/data/TCGA/SV/bedpe/00493087-9d9d-40ca-86d5-936f1b951c93.broad-dRanger_snowman.20150918.somatic.sv.bedpe'
    # cluster_path = '~/Documents/HMS/PhD/beroukhimlab/broad_mount/youyun/complexSV/data/TCGA/SV/bedpe/00493087-9d9d-40ca-86d5-936f1b951c93.broad-dRanger_snowman.20150918.somatic.sv.sv_clusters_and_footprints.tsv'
    input_bedpe = fread(input_path, sep = "\t",header = FALSE)
    cluster_bedpe = fread(cluster_path, sep = "\t",header = FALSE)
    print(paste0("Input bedpe file has ", nrow(input_bedpe), " rows"))

    colnames(input_bedpe) = c(
        'seqnames_breakend1', 'start_breakend1', 'end_breakend1',
        'seqnames_breakend2', 'start_breakend2', 'end_breakend2',
        'Sample', 'SV_ID', 'cnt_type_breakend1', 'cnt_type_breakend2',
        'ins_seq', 'ins_len', 'mh_seq', 'mh_len', 'N_ALT', 'N_ALT_RP', 'N_ALT_SR', 'svclass',
        'REF_breakend1', 'REF_breakend2', 'ALT_breakend1', 'ALT_breakend2', 'ID_breakend1', 'ID_breakend2'
    )
    colnames(cluster_bedpe) = c(
        'seqnames_breakend1', 'start_breakend1', 'end_breakend1',
        'seqnames_breakend2', 'start_breakend2', 'end_breakend2',
        'Sample', 'SV_ID', 'cnt_type_breakend1', 'cnt_type_breakend2',
        'cluster_ID','cluster_size','footprint_ID_low','footprint_ID_high',
        'span_low','span_high','pval'
    )

    # merge the two bedpe files -----------------------------------------------------------------------------------------------------------
    # Not all SVs maybe preserved in the cluster SV output according to the github page:
    #       "SVs with at least one breakpoint within the provided telomere and centromere boundaries are discarded. 
    #           Therefore the output may have fewer SVs than the input."
    merged_bedpe = merge(cluster_bedpe, input_bedpe, by = c(
        'seqnames_breakend1', 'start_breakend1', 'end_breakend1',
        'seqnames_breakend2', 'start_breakend2', 'end_breakend2',
        'Sample', 'SV_ID', 'cnt_type_breakend1', 'cnt_type_breakend2'
    ), all.x = TRUE)
    print(paste0("Merged bedpe file has ", nrow(merged_bedpe), " rows"))

    # load the annotation files -----------------------------------------------------------------------------------------------------------
    replication_timing = readRDS(paste0(workdir, "Jeremiah/GenomeTracks/gr.repl100k.rds"))
    gc_content = fread(paste0(workdir, "Jeremiah/GenomeTracks/GCcontent/hg19.1000.gc5.bed"), header = FALSE)
    gc_content = GRanges(gc_content$V1, IRanges(gc_content$V2, gc_content$V3), strand = "*", gc = gc_content$V4)
    gene_ranges = readRDS(paste0(workdir, "Jeremiah/tracks/gr.allgenes.rds"))
    repetitive_elements = readRDS(paste0(workdir, "Jeremiah/tracks/gr.repeatMasker.rds"))
    telo_centro = fread(
        cmd = paste0("grep -v '^#' ",workdir, "youyun/complexSV/code/ClusterSV/references/hg19_centromere_and_telomere_coords.txt"), 
        header = TRUE, sep = "\t"
    )[, c('qtel','note') := data.table(str_split_fixed(qtel, "\\#",2))]
    telo_granges = with(melt(
            telo_centro[,.(chr,ptel,qtel,note)], id.vars = c('chr','note'), 
            variable.name = 'telomere', value.name = 'start'
        )[, end := start], GRanges(
            chr, IRanges(as.numeric(start), as.numeric(end)), strand = "*",
            note = note, telomere = telomere
        )
    )
    cen_granges = with(melt(
            telo_centro[,.(chr,cen_start,cen_end,note)], id.vars = c('chr','note'), 
            variable.name = 'centromere', value.name = 'start'
        )[, end := start], GRanges(
            chr, IRanges(as.numeric(start), as.numeric(end)), strand = "*",
            note = note, centromere = centromere
        )
    )
    TADs = fread(paste0(workdir, "Jeremiah/tracks/TADannotations.bed"), header = FALSE)[, V1 := gsub("chr", "", V1)]
    TADs = rbind(TADs[,.(V1, V2, V4)], TADs[,.(V1, V2 = V3, V4)])[, V3 := V2]
    TADs = GRanges(TADs$V1, IRanges(TADs$V2, TADs$V3), strand = "*", annotation = TADs$V4)
    CNAs = fread(list.files(
        paste0(workdir,'youyun/complexSV/data/TCGA/CNA'),
        gsub(".*/|\\..*", "", input_path), full.names = TRUE
    ))
    CNAs = GRanges(
        CNAs$chromosome, IRanges(CNAs$start, CNAs$end), strand = "*",
        total_cn = CNAs$total_cn, major_cn = CNAs$major_cn, minor_cn = CNAs$minor_cn
    )

    # What are the values of the annotation? -----------------------------------------------------------------------------------------------------------
    summary(replication_timing$score)
    summary(as.numeric(gc_content$gc))
    sort(table(repetitive_elements$repeat.element))
    
    # Getting started with the annotation now -----------------------------------------------------------------------------------------------------------
    merged_bedpe[
        , c(
            'rep_time_breakend1', 'gc_breakend1', 'gene_density_breakend1',
            'dist_to_line_breakend1', 'dist_to_sine_breakend1', 'dist_to_ltr_breakend1',
            'dist_to_dna_breakend1', 'dist_to_sr_breakend1', 'dist_to_TAD_breakend1',
            'dist_to_telomere_breakend1', 'dist_to_centromere_breakend1',
            'total_cn_breakend1', 'major_cn_breakend1', 'minor_cn_breakend1',
            'rep_time_breakend2', 'gc_breakend2', 'gene_density_breakend2',
            'dist_to_line_breakend2', 'dist_to_sine_breakend2', 'dist_to_ltr_breakend2',
            'dist_to_dna_breakend2', 'dist_to_sr_breakend2', 'dist_to_TAD_breakend2',
            'dist_to_telomere_breakend2', 'dist_to_centromere_breakend2',
            'total_cn_breakend2', 'major_cn_breakend2', 'minor_cn_breakend2'
        ) := data.table(t(c(
            annotate_around_coordinate(
                seqnames_breakend1, start_breakend1, replication_timing, gc_content, gene_ranges, 
                repetitive_elements, telo_granges, cen_granges, TADs, CNAs
            ),
            annotate_around_coordinate(
                seqnames_breakend2, start_breakend2, replication_timing, gc_content, gene_ranges, 
                repetitive_elements, telo_granges, cen_granges, TADs, CNAs
            )
        ))), .(SV_ID)
    ]

    # write the output -----------------------------------------------------------------------------------------------------------
    write.table(merged_bedpe, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
