### Create an instance of this class (see below for helper functions)
# 250bp bins; single end reads; MAPQ>=10; no duplicates
countConfig <- new("NormRCountConfig", 
                   type = "single.end", binsize = 250L, mapq = 10L, 
                   filteredFlag = 1024L, shift = 0L, midpoint = FALSE, 
                   tlenFilter = NULL)

### Counting configuration for Single End alignment files
# 250bp bins (default); only reads with MAPQ>=20; move counting position 100bp
countConfigurationSE <- countConfigSingleEnd(mapq = 20, shift = 100)
countConfigurationSE

### Counting configuration for Paired End alignment files
# 250bp bins; count center of fragments; only fragments with 70bp<=length<=200
countConfigurationPE <- countConfigPairedEnd(midpoint = TRUE,
                                             tlenFilter = c(70, 200))
countConfigurationPE
