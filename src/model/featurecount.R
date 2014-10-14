## A list `railroadtracks_import` will be needed for this to run
## That list will contain the following elements:
## alignedreads_fn: file name for SAM/BAM file with alignments.
## annotation_fn: file name for GFF/GTF file with annotation corresponding to the reference
## ispairedend: a boolean to indicate if the data aligned are pair-ends
## strandspecific: 0 (no), 1 (yes), 2 (yes, reverse-stranded)
## gtf_featuretype: (default: "exon")
## gtf_attrtype: (default: "gene_id")
##
## results_fn: file name for results (CSV file)

require("Rsubread")

run <- local({
    ## Return R connection, deciding on whether it is gzip-compressed
    ## according to the file extension
    autoopen <- function(filename, open) {
        if (grepl("\\.gz$", filename)) {
                                        # gzip-compressed file
            conn <- gzfile(filename, open=open)
        } else {
            conn <- file(filename, open=open)
        }
        return(conn)
    }

    run <- function(railroadtracks_import) {
        p <- railroadtracks_import
        res <- Rsubread::featureCounts(files = p$alignedreads_fn,
                                       annot.ext = p$annotation_fn,
                                       isGTFAnnotationFile = TRUE,
                                       isPairedEnd = p$ispairedend,
                                       strandSpecific=p$strandspecific,
                                       GTF.featureType = p$gtf_featuretype,
                                       GTF.attrType = p$gtf_attrtype)
        counts <- res$counts
        #FIXME: should the tabular data follow a more defined format ?
        out_table <- as.data.frame(counts)
        ## Fix the column names (counts go the in column "count")
        out_table <- data.frame(ID=rownames(out_table),
                                count=out_table[[1]])
        out_conn <- autoopen(p$results_fn, "w")
        write.csv(out_table, file=out_conn, row.names=FALSE)
        close(out_conn)
    }

    run
})
