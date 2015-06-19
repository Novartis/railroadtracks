# Copyright 2014-2015 Novartis Institutes for Biomedical Research

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
## dispersion: numerical value [optional]
### FIXME: format for results ? should be common to all differential expression results
## results_fn: file name for results

require("edgeR")

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

    make_dgelist <- function(counttable_fn, sampleinfo_fn) {
        m <- load_csv(counttable_fn, row_names=1)
        m <- as.matrix(m)
        sampleinfo <- load_csv(sampleinfo_fn, row_names="sample_id")
        d <- edgeR::DGEList(counts=m, group=sampleinfo$group)#,
                            #lib.size=sampleinfo$libsize)
        return(d)
    }

    run <- function(railroadtracks_import) {
        p <- railroadtracks_import
        dgl <- make_dgelist(p$counttable_fn, p$sampleinfo_fn)
        ## 'dispersion' as a parameter
        d <- p$dispersion
        if (is.null(d)) {
            d <- 0.2 # FIXME value use in the R man page. Does having a default make sense ?
        }        
        de <- edgeR::exactTest(dgl, dispersion=d)
        res <- topTags(de, n=Inf)
        out_conn <- autoopen(p$diffexp_fn, "w")
        write.csv(res, file=out_conn)
        close(out_conn)
    }

    run
})
