# Copyright (C) 2017 Johannes Helmuth & Ho-Ryun Chung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
setClassUnion("integerOrNULL", c("integer", "NULL"))

#' Container for configuration of read counting with bamsignals in normR
#'
#' This S4 class is a small wrapper for a configuration on obtaining counts
#' from bamfiles with \code{\link[bamsignals]{bamProfile}}. Herein, two
#' functions provide help for creating an instance of this class:
#' \code{\link{countConfigSingleEnd}} creates a configuration for single
#' end reads; and \code{\link{countConfigPairedEnd}} creates a configuration
#' for paired end reads.
#'
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#'
#' @slot type A \code{character} of value \code{paired.end} or
#' \code{single.end}.
#' @slot binsize An \code{integer} specifying the binsize in bp.
#' @slot mapq An \code{integer} specifying the minimal mapping quality for a
#' read to be counted.
#' @slot filteredFlag An \code{integer} to filter for in the SAMFLAG field.
#' For example, 1024 filters out marked duplicates (default). Refer to
#' \url{https://broadinstitute.github.io/picard/explain-flags.html} for
#' details.
#' @slot shift An \code{integer} specifing a shift of the read counting
#' position in 3'-direction. This can be handy in the analysis of chip-seq
#' data.
#' @slot midpoint Paired End data only: A \code{logical} indicating whether
#' fragment midpoints instead of 5'-ends should be counted.
#' @slot tlenFilter Paired End data only: An \code{integer} of length two
#' specifying the lower and upper length bound for a fragment to be considered.
#' The fragment length as estimated from alignment in paired end experiments
#' and written into the TLEN column.
#'
#' @return A \code{\link{NormRCountConfig}} with specified counting parameters
#' for \code{\link{normr}} methods (\code{\link{enrichR}}, \code{\link{diffR}},
#' \code{\link{regimeR}}
#'
#' @aliases NormRCountConfig BamsignalsConfig BamsignalsCountConfig
#' @seealso \link{normr} for functions that use this object.
#'
#' @example inst/examples/NormRCountConfig_example.R
#'
#' @include methods.R
#' @export
setClass("NormRCountConfig",
    representation = representation(type="character",
                                    binsize="integer",
                                    mapq="integer",
                                    filteredFlag="integer",
                                    shift="integer",
                                    midpoint="logical",
                                    tlenFilter="integerOrNULL"
                                    )
)
setValidity("NormRCountConfig",
    function(object) {
      if (!(object@type %in% c("single.end", "paired.end"))) {
        stop("invalid type")
      }
      if (object@binsize <= 0) stop("invalid binsize")
      if (object@mapq < 0 | object@mapq > 255) stop("invalid mapq")
      if (object@filteredFlag < -1 | object@filteredFlag > 4095) {
        stop("invalid filteredFlag")
      }
      if (object@shift < 0) stop("invalid shift")
      if (object@type == "paired.end") {
        if (length(object@tlenFilter) != 2) stop("invalid tlenFilter")
        if (object@tlenFilter[1] < 0) stop("invalid tlenFilter lower bound")
      }
      TRUE
    }
)

setGeneric("countConfigSingleEnd", function(...)
  standardGeneric("countConfigSingleEnd"))
#' @describeIn NormRCountConfig Setup single end count configuration
#'
#' @param binsize An \code{integer()} specifying the binsize in bp.
#' @param mapq An \code{integer()} specifying the minimal mapping quality for a
#' read to be counted.
#' @param filteredFlag An \code{integer()} to filter for in the SAMFLAG field.
#' For example, 1024 filters out marked duplicates (default). Refer to
#' \url{https://broadinstitute.github.io/picard/explain-flags.html} for
#' details.
#' @param shift An \code{integer()} specifing a shift of the read counting
#' position in 3'-direction. This can be handy in the analysis of chip-seq
#' data.
#'
#' @aliases countConfigSingleEnd configSingleEnd
#'
#' @export
setMethod("countConfigSingleEnd",
  definition=function(binsize=250L, mapq=20L, filteredFlag=1024L, shift=0L) {
    new("NormRCountConfig", type="single.end", binsize=as.integer(binsize),
        mapq=as.integer(mapq), filteredFlag=as.integer(filteredFlag),
        shift=as.integer(shift), tlenFilter=NULL)
})

setGeneric("countConfigPairedEnd", function(...)
  standardGeneric("countConfigPairedEnd"))
#' @describeIn NormRCountConfig Setup paired end count configuration
#'
#' @param midpoint Paired End data only: A \code{logical()} indicating whether
#' fragment midpoints instead of 5'-ends should be counted.
#' @param tlenFilter An \code{integer()} of length two specifying the lower and
#' upper length bound for a fragment to be considered. The fragment length as
#' estimated from alignment in paired end experiments and written into the TLEN
#' column.
#'
#' @aliases countConfigPairedEnd configPairedEnd
#'
#' @export
setMethod("countConfigPairedEnd",
  definition=function(binsize=250L, mapq=20L, filteredFlag=1024L, shift=0L,
                      midpoint=TRUE, tlenFilter=c(70L, 200L)) {
    new("NormRCountConfig", type="paired.end", binsize=as.integer(binsize),
        mapq=as.integer(mapq), filteredFlag=as.integer(filteredFlag),
        shift=as.integer(shift), midpoint=as.logical(midpoint),
        tlenFilter=as.integer(tlenFilter))
})

setGeneric("getFilter", function(x) standardGeneric("getFilter"))
#' @describeIn NormRCountConfig Get the filter compliant to
#' \code{\link[bamsignals]{bamProfile}}
#'
#' @aliases getFilter getBamsignalsFilter getNormRCountConfigFilter
#'
#' @export
setMethod("getFilter", "NormRCountConfig",
    function(x) {
      if (x@type == "single.end") return("ignore")
      if (x@type == "paired.end" & !x@midpoint) return("filter")
      if (x@type == "paired.end" & x@midpoint) return("midpoint")
    }
)

#' @describeIn NormRCountConfig Prints a given BamCounConfig
#'
#' @param x A \code{NormRCountConfig} object.
#'
#' @export
setMethod("print", "NormRCountConfig", function(x, ...) {
      cat("NormRCountConfig-class object\n\n",
          "Type:\t\t\t", x@type, "\n",
          "Binsize:\t\t", x@binsize, "bp\n",
          "Minimum MAPQ:\t\t", x@mapq, "\n",
          "Filtered SAMFLAG:\t", x@filteredFlag, "\n",
          "Shift of anchor:\t", x@shift, "bp \n")
      if (x@type == "paired.end") {
        cat("\nPaired End Options: \n",
            "\tMidpoint counting:\t", x@midpoint, "\n",
            "\tMinimal Fragmentlength:\t", x@tlenFilter[1], "bp\n",
            "\tMaximal Fragmentlength:\t", x@tlenFilter[2], "bp\n"
            )
      }
      invisible(x)
    }
)

#' @describeIn NormRCountConfig Shows a given BamCounConfig
#'
#' @param object A \code{NormRCountConfig} object.
#' @param ... optional arguments to be passed directly to the inherited
#' function without alteration and with the original names preserved.
#'
#' @export
setMethod("show", "NormRCountConfig", function(object) print(object))
