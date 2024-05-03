suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

# Helper functions -----------------------------------------------------------------------------------------------------------
collapse_VCF_to_bedpe = function(vcf_data){
    # minimal columns needed -- Sample, seqnames, start, cnt_type, ID
    vcf_data[, c("SV_ID", "breakend_num") :=
        list(
            paste0(Sample, "__", gsub(":[0-2]$", "", ID)),
            gsub(".*:", "", ID)
    )]
    print(paste0('There are ',nrow(vcf_data[,.N,SV_ID][N != 2]),' SVs without two breakend.'))
    # there are some SVs without two breakends, we will discard them
    vcf_data = vcf_data[SV_ID %in% vcf_data[,.N,SV_ID][N == 2]$SV_ID]
    # per sv pair, we sort the breakend numbers, sort them, and give them a unified order
    vcf_data[order(breakend_num), breakend_order := c("breakend1", "breakend2"), SV_ID]
    # SV_config_combo
    vcf_data[,'SV_config_combo' := paste0(cnt_type,collapse = ''),SV_ID]
    # dcast into the short format with breakend 1 and 2 identified by SV_ID
    bedpe_data <- dcast(
        # melt into long format
        melt(
            vcf_data,
            id.vars = c(
                "SV_ID", "Sample", "breakend_order", "SV_config_combo",
                'ins_seq', 'ins_len', 'mh_seq', 'mh_len', 
                'N_ALT', 'N_ALT_RP', 'N_ALT_SR'
            )
        )[, variable := paste0(variable, "_", breakend_order)][, breakend_order := NULL],
        SV_ID + Sample + SV_config_combo + ins_seq + ins_len + 
            mh_seq + mh_len + N_ALT + N_ALT_RP + N_ALT_SR ~ variable,
        value.var = "value"
    )
    # there are some bugs from old code that we will ignore here
    if(nrow(bedpe_data[grepl('DUP|DEL',SV_ID) & SV_config_combo %in% c('--','++')])/nrow(bedpe_data) < 0.0002){
        print(paste0(
                "There are some bugs from old code that we will ignore here. ",
                "The number of deletions and duplications with SV_config_combo == '--' or '++' is ",
                nrow(bedpe_data[grepl('DUP|DEL',SV_ID) & SV_config_combo %in% c('--','++')])/nrow(bedpe_data)*100,
                "% of the total SVs."
        ))
        bedpe_data = bedpe_data[!(grepl('DUP|DEL',SV_ID) & SV_config_combo %in% c('--','++'))]
    }else{
        stop('too many bugs from old code')
    }
    bedpe_data[,svclass := ifelse(
        seqnames_breakend1 == seqnames_breakend2, 
        ifelse(
            cnt_type_breakend1 == cnt_type_breakend2, 
            ifelse(cnt_type_breakend1 == '+', "h2hINV", "t2tINV"), 
            ifelse(cnt_type_breakend1 == "+",  "DEL", "DUP")
        ),
        "TRA"
    )]
    bedpe_data$svclass = factor(bedpe_data$svclass, levels = c("DEL", "DUP", "h2hINV", "t2tINV", "TRA"))
    bedpe_data
}


if (!interactive()) {
    option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = "",
        help = "Input SvABA VCF File", metavar = "input"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = "",
        help = "Output BEDPE file path to write to", metavar = "output"
    )
    # make_option(c("-t", "--matchpen"),
    #   type = "numeric", default = 3,
    #   help = "match penalty for alignment", metavar = "match_pen"
    # )
    )

    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    vcf_path <- opt$input
    output_path <- opt$output
    print(paste0("Input VCF file: ", vcf_path))

    # load the VCF file -----------------------------------------------------------------------------------------------------------
    # vcf_path = '~/Documents/HMS/PhD/beroukhimlab/broad_mount/youyun/nti/data/tcga/dr_sm_vcf/d9c5493e-f969-4c04-a646-9a3134011021.broad-dRanger_snowman.20150918.somatic.sv.vcf'
    vcf = fread(cmd = paste("grep -v '^#'", vcf_path), sep = "\t")
    colnames11 <- c("seqnames", "start", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "SUPPORT_TYPE",'SUPPORT_TYPE_CLASS','SUPPORT_TYPE_VALUE')
    colnames(vcf) <- colnames11

    # Filter for something called by SvABA at least -----------------------------------------------------------------------------------------------------------
    vcf <- vcf[grepl("CALLER=S|CALLER=SD", INFO, ignore.case = TRUE)]
    vcf[
        ,both_end_pass_filter := .N == 2,.(gsub(':[1-2]$','',ID))
    ]
    vcf <- vcf[both_end_pass_filter == TRUE]

    # make sure chromosome format are correct (1-22, X, Y) -----------------------------------------------------------------------------------------------------------
    vcf$seqnames <- gsub("chr", "", vcf$seqnames)
    vcf$ALT <- gsub("chr", "", vcf$ALT)

    # Extract Information -----------------------------------------------------------------------------------------------------------
    # Sample
    vcf$Sample <- gsub(".*/|\\..*", "", vcf_path)
    # Extract the insertion sequence
    vcf$ins_seq = gsub("INSERTION=|FORSEQ=|SVINSSEQ=", "", str_extract(vcf$INFO, "INSERTION=[ACGT]*|FORSEQ=[ACGT]*|SVINSSEQ=[ACGT]*"))
    vcf$ins_len = nchar(vcf$ins_seq)
    # Extract the MH length
    vcf$mh_seq <- gsub("HOMSEQ=", "", str_extract(vcf$INFO, "HOMSEQ=[ACGT]*"))
    vcf$mh_len = nchar(vcf$mh_seq)
    # Extract the cnt type
    vcf$cnt_type <- unlist(lapply(grepl("^[AGCT]", vcf$ALT), function(x) {
        ifelse(x, "+", "-")
    }))
    # read support (N_ALT sometimes doesn't equal to N_ALT_RP + N_ALT_SR)
    vcf[,
        c('N_ALT','N_ALT_RP','N_ALT_SR') := list(
            as.numeric(strsplit(SUPPORT_TYPE_VALUE, ":")[[1]][1]),
            as.numeric(strsplit(SUPPORT_TYPE_VALUE, ":")[[1]][2]),
            as.numeric(strsplit(SUPPORT_TYPE_VALUE, ":")[[1]][3]) 
        ), .(ID)
    ]

    # Collapse into BEDPE format -----------------------------------------------------------------------------------------------------------
    bedpe = collapse_VCF_to_bedpe(vcf[,.(
        # breakend level features
        seqnames, start, ID, REF, ALT, QUAL, FILTER, INFO, cnt_type,
        # SV level features
        Sample, ins_seq, ins_len, mh_seq, mh_len, N_ALT, N_ALT_RP, N_ALT_SR
    )])
    # The columns resulting from the above processing are as follows:
    # SV_ID, Sample, SV_config_combo, ins_seq, ins_len, mh_seq, mh_len, N_ALT, N_ALT_RP, N_ALT_SR, 
    # ALT_breakend1, ALT_breakend2, FILTER_breakend1, FILTER_breakend2, ID_breakend1, ID_breakend2, 
    # INFO_breakend1, INFO_breakend2, QUAL_breakend1, QUAL_breakend2, REF_breakend1, REF_breakend2, 
    # breakend_num_breakend1, breakend_num_breakend2, cnt_type_breakend1, cnt_type_breakend2, 
    # seqnames_breakend1, seqnames_breakend2, start_breakend1, start_breakend2, svclass

    # We need to follow the standard bedpe format then put the other columns afterwards -----------------------------------------------------------------------------------------------------------

    ## standard bedpe format columns
    # seqnames_breakend1, start_breakend1, end_breakend1, seqnames_breakend2, start_breakend2, end_breakend2,
    # Sample, SV_ID, cnt_type_breakend1, cnt_type_breakend2
    
    ## Other important columns
    # ins_seq, ins_len, mh_seq, mh_len, N_ALT, N_ALT_RP, N_ALT_SR, svclass
    
    ## the rest of the columns I will keep
    # REF_breakend1, REF_breakend2, ALT_breakend1, ALT_breakend2, ID_breakend1, ID_breakend2
    
    ## the rest of the columns I will not keep
    # FILTER_breakend1, FILTER_breakend2, INFO_breakend1, INFO_breakend2, QUAL_breakend1, QUAL_breakend2,
    # breakend_num_breakend1, breakend_num_breakend2, SV_config_combo

    write.table(
        bedpe[,.(
            seqnames_breakend1, start_breakend1, end_breakend1 = start_breakend1,
            seqnames_breakend2, start_breakend2, end_breakend2 = start_breakend2,
            Sample, SV_ID, cnt_type_breakend1, cnt_type_breakend2,
            ins_seq, ins_len, mh_seq, mh_len, N_ALT, N_ALT_RP, N_ALT_SR, svclass,
            REF_breakend1, REF_breakend2, ALT_breakend1, ALT_breakend2, ID_breakend1, ID_breakend2
        )],
        output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
}
