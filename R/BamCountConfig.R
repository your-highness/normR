# Copyright (C) 2016 Johannes Helmuth
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
#' Container for counting with bamsignals
#'
#' This S4 class is a small wrapper for a configuration on obtaining counts
#' from bamfiles with bamsignals::bamProfile().
#'
#' @field type A \code{character()} representing the type of bamfile
#' ("paired.end", "single.end").
#' @field chromosomes GenomicRanges object used to specify the chromosomes.
#' @field binsize An \code{integer()} giving the binsize to count in (in bp).
#' @field mapq discard reads with mapping quality strictly lower than this
#' parameter. The value 0 ensures that no read will be discarded, the value 254
#' that only reads with the highest possible mapping quality will be considered.
#' @field filteredFlag An \code{integer()} to filter for when counting reads (SAMFLAG).
#' For example, 1024 filters out marked duplicates (default), see
#' \url{https://broadinstitute.github.io/picard/explain-flags.html}
#' @field shift shift the read position by a user-defined number of basepairs.
#' This can be handy in the analysis of chip-seq data.
#' @field midpoint A \code{logical()} indicating whether fragment midpoints
#' instead of 5'-ends should be counted in paired end data.
#' @field tlenFilter A filter on fragment length as estimated from alignment
#' in paired end experiments (TLEN). If set to \code{c(min,max)} only reads are
#' considered where \code{min <= TLEN <= max}. If \code{paired.end=="ignore"},
#' this argument is set to \code{NULL} and no filtering is done. If
#' \code{paired.end!="ignore"}, this argument defaults to \code{c(0,1000)}.
#'
#' @aliases BamCountConfig
#' @aliases BamsignalsConfig
#' @aliases BamsignalsCountConfig
#'
#' @import GenomicRanges
#' @seealso \code{\link{normr-methods}} for the functions that require this
#' object.\cr
#' See \code{\link{bamsignals}}-package for counting in bam files.
#' @return return values are described in the Methods section.
setClassUnion("integerOrNULL", c("integer", "NULL"))

#' @export
setClass("BamCountConfig",
    representation = representation(type="character",
                                    chromosomes="GRanges",
                                    binsize="integer",
                                    mapq="integer",
                                    filteredFlag="integer",
                                    shift="integer",
                                    midpoint="logical",
                                    tlenFilter="integerOrNULL"
                                    )
)

setValidity("BamCountConfig",
    function(object) {
      if (!(object@type %in% c("single.end", "paired.end"))) {
        stop("invalid type")
      }
      if (is.null(object@chromosomes)) stop("invalid chromosomes")
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

setMethod("print", "BamCountConfig",
    function(x) {
      cat("BamCountConfig-class object\n\n",
          "Type:\t\t\t", x@type, "\n",
          "Binsize:\t\t", x@binsize, " bp\n",
          "Minimal MAPQ:\t\t", x@mapq, "\n",
          "Filtered SAMFLAG:\t", x@filteredFlag, "\n",
          "Shift of anchor:\t", x@shift, " bp \n")
      if (x@type == "paired.end") {
        cat("\nPaired End Options: \n",
            "\tMidpoint counting:\t", x@midpoint, "\n",
            "\tMinimal Fragmentlength:\t", x@tlenFilter[1], "\n",
            "\tMaximal Fragmentlength:\t", x@tlenFilter[2], "\n"
            )
      }
      invisible(x)
    }
)
setMethod("show", "BamCountConfig", function(object) print(object))

#' @export
setGeneric("countConfigPairedEnd", function(...)
  standardGeneric("countConfigPairedEnd"))
#' @describeIn BamCountConfig Setup a paired end counting configuration
#' @aliases countConfigPairedEnd
#' @export
setMethod("countConfigPairedEnd",
  definition=function(binsize=250L, mapq=20L, filteredFlag=1024L, shift=0L,
                      midpoint=T, tlenFilter=c(70L, 200L)) {
    new("BamCountConfig", type="paired.end", binsize=as.integer(binsize),
        mapq=as.integer(mapq), filteredFlag=as.integer(filteredFlag),
        shift=as.integer(shift), midpoint=as.logical(midpoint),
        tlenFilter=as.integer(tlenFilter))
})

#' @export
setGeneric("countConfigSingleEnd", function(...)
  standardGeneric("countConfigSingleEnd"))
#' @describeIn BamCountConfig Setup a single end counting configuration
#' @aliases countConfigSingleEnd
#' @export
setMethod("countConfigSingleEnd",
  definition=function(binsize=250L, mapq=20L, filteredFlag=1024L, shift=0L) {
    new("BamCountConfig", type="single.end", binsize=as.integer(binsize),
        mapq=as.integer(mapq), filteredFlag=as.integer(filteredFlag),
        shift=as.integer(shift), tlenFilter=NULL)
})

#' @export
setGeneric("getFilter", function(obj) standardGeneric("getFilter"))
#' @describeIn BamCountConfig Get the bamsignals relevant filter
#' @aliases getFilter
#' @export
setMethod("getFilter", "BamCountConfig",
    function(obj) {
      if (obj@type == "single.end") return("ignore")
      if (obj@type == "paired.end" & !obj@midpoint) return("filter")
      if (obj@type == "paired.end" & obj@midpoint) return("midpoint")
    }
)
