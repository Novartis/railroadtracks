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
## maxround: numerical value [optional]
## qtrm: Qtrm in EBSeq::EBTest
## qtrmcut: QtrmCut in EBSeq::EBTest
### FIXME: format for results ? should be common to all differential expression results
## results_fn: file name for results

require("EBSeq")

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
        # data matrix
        m <- load_csv(p$counttable_fn, row_names=1)
        m <- as.matrix(m)
        #
        m_sizes = EBSeq::MedianNorm(m)
        sampleinfo <- load_csv(p$sampleinfo_fn, row_names="sample_id")

        maxround <- with(p, ifelse(is.null(maxround), 5, maxround))
        qtrm <- with(p, ifelse(is.null(qtrm), .75, qtrm))
        qtrmcut <- with(p, ifelse(is.null(qtrmcut), 10, qtrmcut))

        ebt <- EBSeq::EBTest(Data=m,
                             Conditions=sampleinfo$group,
                             sizeFactors=m_sizes,
                             maxround=maxround,
                             Qtrm=qtrm,
                             QtrmCut=qtrmcut)
        # matrix of PP (presumably "posterior probabilities")
        pp <- EBSeq::GetPPMat(ebt)
        # fold-change
        fc <- EBSeq::PostFC(ebt)
        # write all this in a result file
        res <- cbind(as.data.frame(pp),
                     fc[c("PostFC", "RealFC")],
                     ID=rownames(pp))
        out_conn <- autoopen(p$diffexp_fn, "w")
        write.csv(res, file=out_conn)
        close(out_conn)
    }

    run
})
