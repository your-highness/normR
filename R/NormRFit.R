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
#' Container for a NormR fit
#'
#' This S4 class wraps a NormR fit containing counts, run configuration and 
#' results of the fit. 
#'
#' Among other things the container provides an accessor method, 
#' that returns single signals as vectors and matrices, and the 
#' methods \code{as.list} and \code{alignSignals}, that convert the 
#' container to a list or an array/matrix respectively. A CountSignals
#' object is read-only, i.e. it cannot be modified.
#'
#' @param obj A NormRFit object.
#' @param i \code{integer()}-index for subsetting.
#' @param B The index of the background component.
#'
#' @slot type A \code{character()} representing the type of fit ("enrichR", 
#' "diffR", "regimeR")
#' @slot k An \code{integer()} giving the number of binomial mixture 
#' components to  befit to the data.
#' @slot B The index of the background component.
#' @slot counts A two-dimensional \code{list()} containing to \code{n}-length
#' \code{integer()}-vectors of count data used for the NormR fit.
#' @slot thetaStar A \code{numeric()} representing a naive background 
#' estimation, i.e. \code{sum(counts[[2]])/sum(counts[[1]]+counts[[2]])}
#' @slot theta A \code{k}-length \code{numeric()}-vector representing 
#' NormR-fitted parametrization of binomial mixture components.
#' @slot mixtures A \code{k}-length \code{numeric()}-vector representing 
#' NormR-fitted mixture proportions of binomial mixture
#' components.
#' @slot lnL A \code{numeric()}-vector representing the log-likelihood-trace of
#' the NormR model fit.
#' @slot eps A \code{numeric{}} representing preset threshold for NormR fit 
#' convergence.
#' @slot posteriors A \code{k}-dimensional \code{matrix()} containing 
#' posterior probabilities for every \code{n} 2-tupel in counts given each
#' \code{theta} and \code{mixtures}.
#' @slot enrichment A \code{n}-sized \code{numeric()}-vector of normalized 
#' enrichment over component \code{B} for each counts-tupel 
#' @slot p.vals A \code{n}-size \code{numeric()}-vector containing log-P-values 
#' every \code{n} 2-tupel in counts given \code{theta} of \code{B}.
#' @slot filteredT A \code{integer()}-vector giving indices of P-values to be
#' considered for FDR correction. These indices have been obtainted by filtering
#' P-values using the T method.
#' @slot q.vals A \code{n[filteredT]}-size \code{numeric()}-vector containing 
#' log-q-values for every \code{n[filteredT]} 2-tupel in counts given 
#' \code{theta} of \code{B}. These are P-values corrected for multiple testing 
#' using Storey's method.
#' @aliases NormRFit
#' @seealso \code{\link{normr-methods}} for the functions that produce 
#' this object
#' @return return values are described in the Methods section.
#' @example inst/examples/class_example.R
#' @export
setClass("NormRFit", 
    representation = representation(type = "character",
                                    k = "integer",
                                    B = "integer",
                                    counts = "list", 
                                    thetastar = "numeric",
                                    theta = "numeric",
                                    mixtures = "numeric",
                                    lnL = "numeric",
                                    eps = "numeric",
                                    posteriors = "matrix",
                                    enrichment = "numeric",
                                    p.vals = "numeric",
                                    filteredT = "numeric",
                                    q.vals = "numeric" )
)

setValidity("NormRFit",
    function(obj) {
      if (!(obj@type %in% c("enrichR", "diffR", "regimeR"))) {
          return("invalid type slot")
      }
      if (length(obj@k) == 0 | obj@k <= 0) return("invalid k slot")
      if (obj@B <= 0 || obj@B > k) return("invalid B slot")
      if (length(obj@counts) != 2 || 
          length(obj@counts[[1]]) != length(obj@counts[[2]])) {
        return("invalid counts slot")
      }
      if (length(obj@thetastar) != 1 || obj@thetastar < 0 || 
          obj@thetastar > 1) {
        return("invalid thetastar slot")
      }
      if (length(obj@theta) != obj@k) return("invald theta slot")
      if (length(obj@mixtures) != obj@k) return("invald mixtures slot")
      n <- length(obj@counts[[1]])
      if (NCOL(obj@posteriors) != obj$k || 
          NROW(obj@posteriors) != n) {
        return("invald posterios slot")
      }
      if (length(obj@enrichment) != n) return("invaled enrichment slot")
      if (length(obj@p.vals) != n) return("invaled p.vals slot")
      if (max(obj@filteredT) > n) return("invaled filteredT slot")
      if (length(obj@q.vals) != length(obj@filteredT)) {
        return("invaled q.vals slot")
      }
      TRUE
    }
)

#' @describeIn NormRFit Number of regions analyzed.
#' @aliases length
#' @export
setMethod("length", "NormRFit", function(obj) length(obj@counts[[1]]))



