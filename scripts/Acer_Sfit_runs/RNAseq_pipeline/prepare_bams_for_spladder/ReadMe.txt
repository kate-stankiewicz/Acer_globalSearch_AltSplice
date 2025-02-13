## steps for preparing STAR bam output for SplAdder

# first sort and index the original bams 

# next use one of the sorted, indexed, original bams to get the chr list, should be the same across all bams so only need to use one)
# check that this matches the chrName.txt file output by STAR genomeGenerate (it does match, so can simply use this)

# next split the bams by host and sym

# then sort and index the split bams

# run spladder on each separately


## should delete intermediate bam files to save space
