#!/bin/bash

# Clean up
rm 1-scripts/ase/jobs/PURE-MAP*

# launch scripts for Colosse
for file in $(ls /home1/scratch/jleluyer/hybrids/3-trimmed/*paired.fastq.gz |grep -E "/(AA|GG)"|perl -pe 's/_R[12].paired.f(ast)?q.gz//'|sort -u)

# fo testing
#for file in $(ls 3-trimmed/*paired.fastq.gz|perl -pe 's/_R[12].paired.f(ast)?q.gz//'|sort -u)
do
	base=$(basename "$file")
	toEval="cat 1-scripts/ase/3-star_mapping_pures.pbs | sed 's/__BASE__/$base/g'"; eval $toEval > 1-scripts/ase/jobs/PURE-MAP_"$base".pbs
done


#change jobs header

#Submit jobs
for i in $(ls 1-scripts/ase/jobs/PURE-MAP*pbs); do qsub $i; done
