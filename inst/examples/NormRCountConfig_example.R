### Counting configuration for Single End alignment files
# 250bp bins (default); only reads with MAPQ>=20; move counting position 100bp
countConfigurationSE <- countConfigSingleEnd(mapq = 20, shift = 100)
countConfigurationSE

### Counting configuration for Paired End alignment files
# 250bp bins; count center of fragments; only fragments with 70bp<=length<=200
countConfigurationPE <- countConfigPairedEnd(midpoint = T, 
                                             tlenFilter = c(70, 200))
countConfigurationPE
