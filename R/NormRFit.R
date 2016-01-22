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
#' Container for NormR fits
#'
#' This s4 class is a tiny wrapper around a normal list (stored in the
#' \code{signals} slot) and it is the output of the methods in the
#' bamsignals package. 
#' Among other things the container provides an accessor method, 
#' that returns single signals as vectors and matrices, and the 
#' methods \code{as.list} and \code{alignSignals}, that convert the 
#' container to a list or an array/matrix respectively. A CountSignals
#' object is read-only, i.e. it cannot be modified.
#'
#' @param x A CountSignals object
#' @param i Index for subsetting. It can be a single index as well as
#' a vector of indices.
#' @param drop In case \code{i} is a vector of length 1, after subsetting, 
#' collapse the CountSignal object to a single signal or not.
#' @slot ss A single boolean value indicating whether all
#'     signals are strand-specific or not
#' @slot signals A list of integer vectors (if \code{ss==TRUE}) or of integer
#'     matrices, representing each signal
#' @aliases CountSignals
#' @seealso \code{\link{bamsignals-methods}} for the functions that produce 
#' this object
#' @return return values are described in the Methods section.
#' @example inst/examples/class_example.R
#' @export
setClass( "CountSignals", 
    representation = representation(signals = "list", ss = "logical" )
)

