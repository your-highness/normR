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

#getLauncher <- function(dest="normR"){
#    RscriptPath <- file.path(Sys.getenv("R_HOME"), "bin", "Rscript")
#    shebang <- NULL
#    if (!file.exists(RscriptPath)) {
#        warning("Rscript executable not found at the expected location")
#    } else shebang <- paste0("#!", RscriptPath)
#    
#    normrPath <- path.package("normr")
#    #if the package was loaded with devtools::load_all it is not properly
#    #installed and we need a slightly different launcher
#    if (tryCatch("normr" %in% devtools::dev_packages(), error=function(e) FALSE)){
#        devtoolsPath <- find.package("devtools")
#        loadDevtools <- paste0("library(devtools, lib.loc=\"", dirname(devtoolsPath), "\")")
#        loadnormr <- paste0("devtools::load_all(\"", normrPath, "\", quiet=TRUE)")
#        loadLibs <- paste(sep="\n", loadDevtools, loadnormr)
#    } else {
#        loadLibs <- paste0("library(normr, lib.loc=\"",
#                           dirname(path.package("normr")), "\")")
#    }
#    
#    writeLines( c(shebang, "cat('loading normr\\n')", loadLibs,
#    "normr:::CLI(args=commandArgs(trailingOnly=TRUE), normr:::getProg())"),
#    dest)
#}
#
##the CLI subprograms are:
##segmentCLI, in the file segment.R
##reportCLI, in the file report.R
##getcountsCLI, in the file getcounts.R
##normalizecountsCLI, in the file normalizecounts.R
#getCLIsubprograms <- function(){list(
#    getcounts=list(desc="Produce a counts matrix from several bam files", 
#    fun=getcountsCLI, cliargs=getGetcountsOptions),
#    normalizecounts=list(desc="Normalize several count matrices", 
#    fun=normalizecountsCLI, cliargs=getNormalizeCountsOptions),
#    segment=list(desc="Produce a segmentation and a report",
#    fun=segmentCLI, cliargs=getSegmentOptions),
#    report=list(desc="Produce a report for a given segmentation", 
#    fun=reportCLI, cliargs=getReportOptions))}
#    
#CLI <- function(args, prog){
#    CLIsubprograms <- getCLIsubprograms()
#    if (args[1] %in% names(CLIsubprograms)){
#        CLIsubprograms[[args[1]]]$fun(args[-1], paste(prog, args[1]))
#    } else {
#        print_CLI(prog, CLIsubprograms)
#        quit(status=1)
#    }
#}
#
