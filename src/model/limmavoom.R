## A list `railroadtracks_import` will be needed for this to run
## That list will contain the following elements:
## counttable_fn: file name for CSV file with counts.
##                The data table will only contain counts, an column names
##                will correspond to sample identifiers
## sampleinfo_fn: file name for CSV file with sample information, with one row per sample
##                The data table will only contain at least the following columns:
##                - sample_id: this must correspond to the column names in the file counttable_fn
##                - group: this will contain group information for 2-sample testing
##                - libsize: library size (FIXME: should this )
### FIXME: format for results ? should be common to all differential expression results
## results_fn: file name for results

require("limma")

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

    load_csv <- function(fn, row_names=NULL) {
        conn <- autoopen(fn, open="r")
        cds <- read.csv(conn, row.names=row_names)
        close(conn)
        return(cds)
    }

    run <- function(railroadtracks_import) {
        p <- railroadtracks_import
        dataf <- load_csv(p$counttable_fn, row_names=1)
        dataf_si <- load_csv(p$sampleinfo_fn, row_names="sample_id")
        l <- levels(dataf_si$group)
        # FIXME: test that only 2 levels, or should this be done in the Python wrapper ?
        x = paste(rev(l), sep="-")
        design <- with(dataf_si, model.matrix(~group))
        colnames(design) <- l
        #ct <- limma::makeContrasts(contrasts=x, levels=design)
        v <- limma::voom(dataf, design)
        fit <- lmFit(v, design)
        fit <- eBayes(fit)
        ##topTags(de)
        out_conn <- autoopen(p$diffexp_fn, "w")
        write.csv(topTable(fit, coef=2, number=Inf), file=out_conn)
        close(out_conn)
    }   
    run
})
