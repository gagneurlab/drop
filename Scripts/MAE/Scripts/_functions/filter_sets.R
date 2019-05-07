#
# Author: mertes, baderda
#
# functions to create subsets for:
#	1. baseline filtering (quality)
#	2. exome filtering
#	3. protein affecting filtering
#	4. rare filtering
#	5. unique filtering
#	6. homozygous filtering
#
###############################################################################


filter_sets <- c(
	"filter_data",
	"filter_exome",
	"filter_prot_effect",
	"filter_rare",
	"filter_uniq"
)

VCF_CUTOFFS = list(
    qual = 90,
    mq   = 30,
    gq   = 90,
    dp   = 10,
    dp4  = 10
)

#'
#' returns the nth column of the dp4 column of the vcf file
#'
get_vcf_dp <- function(vcf, col){
	res <- sapply(info(vcf)$DP4, function(x) x[col])
	if(length(res) == 0){
		return(NA)
	}
	return(res)
}

#' filter_vcf_quality
#'
#' filter for the base line criteria
#' on data.table or vcf object
#'
filter_vcf_quality <- function(data, vcf_cutoffs_list=list(qual = 90, mq = 30, gq = 90, dp = 10, dp4 = 10)){
    if(dim(data)[1] == 0){
        return(data)
    }

	if(is.data.table(data)){
		return(subset(data,
					ref != 'N' &
					!is.na(gt) &
					(qual >= vcf_cutoffs_list[['qual']] | is.na(qual)) &
					(mq   >= vcf_cutoffs_list[['mq']]   | is.na(mq)) &
					(dp   >= vcf_cutoffs_list[['dp']]   | dp41 + dp42 + dp43 + dp44 >= vcf_cutoffs_list[['dp4']] |
						(is.na(dp) & is.na(dp41) & is.na(dp42) & is.na(dp43) & is.na(dp44) ))
		))
	}
	if (grepl("VCF", class(data))) {
		# calculate quality by sequence depth
		dpq      <- info(data)$DP    >= vcf_cutoffs_list[['dp']]
		dp4q     <- get_vcf_dp(data, 1) + get_vcf_dp(data, 2) + get_vcf_dp(data, 3) + get_vcf_dp(data, 4) >= vcf_cutoffs_list[['dp4']]
		dp_is_na <- is.na(info(data)$DP) & is.na(get_vcf_dp(data, 1)) & is.na(get_vcf_dp(data, 2)) & is.na(get_vcf_dp(data, 3)) & is.na(get_vcf_dp(data, 4))
		refn     <- rowData(data)$REF != 'N'
		qualq    <- rowData(data)$QUAL >= vcf_cutoffs_list[['qual']] | is.na(rowData(data)$QUAL)
		mqq      <- info(data)$MQ    >= vcf_cutoffs_list[['mq']]   | (is.null(info(data)$MQ) || is.na(info(data)$MQ))
		if(length(mqq) == 0){
			mqq <- TRUE
		}
		gt_is_na <- as.logical(apply(as.matrix(is.na(geno(data)$GT)), 1, any))
		gtq      <- as.logical(apply(as.matrix(geno(data)$GQ >= vcf_cutoffs_list[['gq']]), 1,
						FUN = function(x) any(x, na.rm = T)
		))

		# return filtered data
		return(data[
					refn &
					gtq  &
					qualq &
					ifelse(is.na(mqq), TRUE, mqq) &
					(dpq | ifelse(is.na(dp4q), FALSE, dp4q) | dp_is_na) &
					(!gt_is_na & gtq)
		])
	}
	stop("Class type '", class(data), "' is not supported!")
}

#'
#' filter for the exon close bases
#' this could be modified eventually, ignore missing values
#'
filter_exonic <- function(data, splice_distance = 5){
    if(dim(data)[1] == 0){
        return(data)
    }

    return(
        subset(data, grepl(".*splice.*", mstype, perl = TRUE)
                | exonic
                | sdistl < splice_distance
                | sdistr < splice_distance
                | (is.na(exonic) & is.na(sdistl) & is.na(sdistr))
        )
    )
}

#'
#' filter for the protein effecting mutations (altering the protein structure)
#'
filter_prot_effect <- function(data){
    if(dim(data)[1] == 0){
        return(data)
    }

	prot_affecting_ms_types <- c('ablation','del','frame-shift','incomplete','init','ins',
			'missense','splice','stop','stop_retain','unstart','unstop')

	if(is.data.table(data)){
    	return(
        	subset(data, (mstype == 'coding' & mtype == 'strucVar')
            	    | mstype %in% prot_affecting_ms_types
        	)
    	)
	}
	if (grepl("VCF", class(data))) {
		return(data[(
				#(mcols(data)$mstype == 'coding' & mcols(data)$mtype == 'strucVar')
				 mcols(data)$mstype %in% prot_affecting_ms_types
		)])
	}
	stop("Class type '", class(data), "' is not supported")
}

#'
#' filter for rare mutations based on database annotations
#'
#' by default use only ExAC database for filtering
#'
filter_rare <- function(data, col_maf = 'MAX_AF', maf_cutoff=0.001) {
    if(dim(data)[1] == 0){
        return(data)
    }

	# convert factor to double if needed
	if(sapply(data, class)[[col_maf]] != "numeric"){
		data[,c(col_maf):= list(as.double(as.character(get(col_maf))))]
	}
	
    return(data[get(col_maf) <= maf_cutoff | is.na(get(col_maf))])
}

#'
#' get the unique mutations in the data set
#' unique means the exact same mutation is seen only in one patient
#' this is based on the columns: chr, pos, ref and alt
#'
#' this function should use all data points (not filtered data)
#' since some variants could be filtered out because of low coverage or other QC steps
#'
filter_uniq <- function(data, data_raw = data) {
	sub_all  <- copy(filter_vcf_quality(data_raw, vcf_cutoffs_list = VCF_CUTOFFS)[,list(chr,pos,ref,alt)])
	setkey(sub_all,chr,pos,ref,alt)
	duplicated_vars <- unique(sub_all[duplicated(sub_all)])

	setkey(duplicated_vars,chr,pos,ref,alt)
	setkey(data,chr,pos,ref,alt)

	# uniq vars filtering
	return(data[!duplicated_vars])
}

#'
#' filters for compound heterozyougs mutations
#'
#' compound heteozygous means in this context:
#' 		1) the mutation is homozygous (there is no reference allele)
#' 		2) there are at least two heterozygous mutations in the same gene
#'
#' X chromosome: Is a single heterozygous mutation a compound?
#' Because of inactivation of the second chromosome.
#'
filter_compound_heterozygous <- function(data){
	if(is.data.table(data)){
		# 1) get all homozyougs variants
		homo_vars <- data[,!grepl("^0|0$", gt)]

		# 2) get compount heterozygous mutations
		# you can use a different approach:
		#		check for two mutation of a patient in the same gene
		compound_vars <- (
					duplicated(data[,list(sample,hgncid)]) |
					duplicated(data[,list(sample,hgncid)], fromLast=TRUE)
		)
	} else if(grepl("VCF", class(data))){
		# 1) get all homozyougs variants
		homo_vars <- as.vector(!grepl("^0|0$", geno(data)$GT))

		# 2) get compount heterozygous mutations
		compound_vars <- (
					duplicated(mcols(data)$SYMBOL) |
					duplicated(mcols(data)$SYMBOL, fromLast=TRUE)
		)
	} else {
		stop("Class type '", class(data), "' is not supported")
	}

	# return both variant lists
	return(
			data[homo_vars | compound_vars]
	)
}


filter_only_compound_heterozygous <- function(vdata){
    vdata <- vdata[grepl("^0|0$", gt)]
    compound_vars <- (
            duplicated(vdata[,list(sample,hgncid)]) |
            duplicated(vdata[,list(sample,hgncid)], fromLast=TRUE)
            )
    return(vdata[compound_vars])
}

#'
#' get the homozygous mutations
#'
#' means all genomic positions where there is no reference allele
#'
filter_homozygous <- function(data){
	return(data[!grepl("^0|0$", gt)]) # no 0 at start or end
}


filter_potential_biallelic <- function(vdata){
    unique(rbind(
        filter_only_compound_heterozygous(vdata),
        filter_homozygous(vdata)
    ))
}


get_qc_rare_protein_biallelic <- function(vdata){
    unique(minimize_factors_of_datatable(
        filter_potential_biallelic(
            filter_prot_effect(
                filter_rare(
                    filter_vcf_quality(
                        vdata,
                        vcf_cutoffs_list = VCF_CUTOFFS
                    ) [!is.na(hgncid)],
                    max_maf = MAX_MINOR_ALLELE_FREQ,
                    exacOnly = TRUE
                )[sample_freq < MAX_VARIANT_FREQ]
            )
        )
    ))
}



#'
#' removes the common variations within a database and within the given dataset
#'
#' it also remove sequencing errors since they should be common amoung the given dataset
#'
remove_common_vars <- function(data, max_maf=0.3, dummy=TRUE){
	data <- remove_seq_errors(data, MAX_MUT_FREQ)
	new_data <- filter_rare(data, max_maf, dummy)

	setkey(data,chr,pos)
	setkey(new_data,chr,pos)

	all_vars <- nrow(unique(data[,list(chr,pos)]))
	removed_vars <- all_vars - nrow(unique(new_data[,list(chr,pos)]))

	print_log(paste("removing ", removed_vars, " variations of ", all_vars,
					" which is ", removed_vars/all_vars*100, "%.", sep=""))

	return(new_data)
}

#'
#' a merge function for the annotation process during VEP annotation
#'
merge_vep_output_with_data <- function(vep_output, data){

	keys_data <- key(data)

	col_vep_output <- c("chr", "pos", "ref", "alt", "mstype", "mtype", "exonic",
			"tgmaf", "nhmaf", "sift1", "pph1", "rsid", "ensgid", "enstid",
			"hgncid", "ccds", "sever"
	)

	remove_col_from_data <- c(
			"mstype", "mtype", "exonic",
			"tgmaf", "nhmaf", "sift1", "pph1", "rsid", "ensgid", "enstid",
			"hgncid", "ccds"
	)

	setkey(vep_output, chr, pos, ref, alt)
	setkey(data, chr, pos, ref, alt)

	print_log("merging vep output with data")

	new_data <- vep_output[data][,list(
			sample=sample,
			chr=chr,
			pos=pos,
			ref=ref,
			alt=alt,
			gt=gt,
			gq=gq,
			hgncid=ifelse(is.na(hgncid), i.hgncid, hgncid),
			sift1=ifelse(is.na(sift1), i.sift1, sift1),
			pph1=ifelse(is.na(pph1), i.pph1, pph1),
			uniprot=uniprot,
			ucsc=ucsc,
			ccds=ifelse(is.na(ccds), i.ccds, ccds),
			refseq=refseq,
			qual=qual,
			dp=dp,
			dp41=dp41,
			dp42=dp42,
			dp43=dp43,
			dp44=dp44,
			mq=mq,
			sp=sp,
			pl=pl,
			ensgid=ifelse(is.na(ensgid), i.ensgid, ensgid),
			enstid=ifelse(is.na(enstid), i.enstid, enstid),
			exonic=ifelse(is.na(exonic), i.exonic, exonic),
			tunum=tunum,
			sdistl=i.sdistl,
			sdistr=i.sdistr,
			mtype=ifelse(is.na(mtype), i.mtype, mtype),
			mstype=ifelse(is.na(mstype), i.mstype, mstype),
			hmaf=hmaf,
			nhmaf=ifelse(is.na(nhmaf), i.nhmaf, nhmaf),
			tgmaf=ifelse(is.na(tgmaf), i.tgmaf, tgmaf),
			rsid=ifelse(is.na(rsid), i.rsid, rsid),
			rspm=rspm,
			rsg5a=rsg5a,
			rsg5=rsg5,
                        sever=sever
	)]

	setkeyv(new_data, keys_data)
	return(new_data)
}

