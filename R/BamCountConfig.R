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
#' from bamfiles with bamsignals.
#'
#' @slot type A \code{character()} representing the type of bamfile 
#' ("paired.end", "single.end")
#' @slot binsize An \code{integer()} giving the binsize to count in
#' @slot flag An \code{integer()} to filter for when counting reads (SAMFLAG).
#' 
#' @aliases BamCountConfig
#' @seealso \code{\link{normr-methods}} for the functions that require this 
#' object\cr
#' \code{\link{bamsignals}} for counting in bam files.
#' @return return values are described in the Methods section.
#' @export
setClass("BamCountConfig", 
    representation = representation(type="character",
                                    binsize="integer",
                                    mapq="integer",
                                    flag="character",
                                    shift="integer",
                                    midpoint="logical",
                                    tlenfilter="integer"
                                    )
)

setValidity("BamCountConfig",
    function(obj) {
      if (!(obj@type %in% c("single.end", "paired.end"))) stop("invalid type")
      if (obj@binsize <= 0) stop("invalid binsize")
      if (obj@mapq < 0 | obj@mapq > 255) stop("invalid mapq")
      if (!is.null(obj@flag)) stop("flag filtering not supported yet")
      if (obj@shift < 0) stop("invalid shift")
      if (obj@type == "paired.end") {
        if (length(obj@tlenfilter) != 2) stop("invalid tlenfilter")
        if (obj@tlenfilter[1] < 0) stop("invalid tlenfilter lower bound")
      }
      TRUE
    }
)

setMethod("print", "BamCountConfig",
    function(obj) {
      cat("BamCountConfig-class object\n\n",
          "Type:\t\t", obj@type, "\n",
          "Binsize:\t\t", obj@binsize, " bp\n",
          "Minimal MAPQ:\t", obj@mapq, "\n",
          "Filtered SAMFLAG:\t", obj@flag, "\n",
          "Shift of anchor:\t", obj@shift, " bp \n")
      if (obj@type == "paired.end") {
        cat("\nPaired End Options: \n",
            "\tMidpoint counting:\t", obj@midpoint, "\n",
            "\tMinimal Fragmentlength:\t", obj@tlenfilter[1], "\n",
            "\tMaximal Fragmentlength:\t", obj@tlenfilter[2], "\n"
            )
      }
      invisible(obj)
    }
)
setMethod("show", "BamCountConfig", function(obj) print(obj))

#' @export
setGeneric("countConfigPairedEnd", function(...) standardGeneric("countConfig"))
#' @describeIn BamCountConfig Setup a paired end counting configuration
#' @aliases countConfigPairedEnd
#' @export
setMethod("countConfigPairedEnd", 
    definition=function(binsize=250, mapq=20, flag=NULL, shift=0, midpoint=T,
                        tlenfilter=c(70, 200)) {
      new("BamCountConfig", type="paired.end", binsize=binsize, mapq=mapq, 
          flag=flag, shift=shift, midpoint=midpoint, tlenfilter=tlenfilter)
})

#' @export
setGeneric("countConfigSingleEnd", function(...) standardGeneric("countConfig"))
#' @describeIn BamCountConfig Setup a single end counting configuration
#' @aliases countConfigSingleEnd
#' @export
setMethod("countConfigSingleEnd", 
    definition=function(binsize=250, mapq=20, flag=NULL, shift=0) {
      new("BamCountConfig", type="single.end", binsize=binsize, mapq=mapq, 
          flag=flag, shift=shift)
})

#' @export
setGeneric("getFilter", function(obj) standardGeneric("getFilter"))
#' @describeIn BamCountConfig
#' @aliases getFilter
#' @export
setMethod("getFilter", "BamCountConfig",
    function(obj) {
      if (obj@type == "single.end") return("ignore")
      if (obj@type == "paired.end" & !obj@midpoint) return("paired.end")
      if (obj@type == "paired.end" & obj@midpoint) return("midpoint")
    }
)

