#!/bin/bash
set -euo pipefail

# script calculates gaps and coverage using original ROI files. This has been requested for
# truSightCancer as the method currently used assumes all exons in a given gene are included
# in the panel - which is often not the case.

# USE: bash inside run folder - will iterate over all samples. Depends on R script
# calculateTargetCoverage.R (located in same folder).

version="2.5.3"

RUN_DIR=$PWD

for i in */;
do

    echo $RUN_DIR/$i
    cd $RUN_DIR/$i

    #load sample & pipeline variables
    . ./*.variables

    . /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

    echo "making PASS file"

    #Make PASS BED
    /share/apps/htslib-distros/htslib-1.4.1/tabix -R /data/results/$seqId/$panel/IlluminaTruSightCancer_CustomROI_b37.bed \
    "$seqId"_"$sampleId"_DepthOfCoverage.gz | \
    awk -v minimumCoverage="$minimumCoverage" '$3 >= minimumCoverage { print $1"\t"$2-1"\t"$2 }' | \
    sort -k1,1V -k2,2n -k3,3n | \
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge > "$seqId"_"$sampleId"_customPASS.bed

    echo "making Gaps file"

    #Make GAP BED
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools subtract \
    -a /data/results/$seqId/$panel/IlluminaTruSightCancer_CustomROI_b37.bed \
    -b "$seqId"_"$sampleId"_customPASS.bed | \
    sort -k1,1V -k2,2n -k3,3n \
    > "$seqId"_"$sampleId"_customGaps.bed

    echo "making clincoverage file"

    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools coverage \
    -a /data/results/$seqId/$panel/IlluminaTruSightCancer_CustomROI_b37.bed \
    -b "$seqId"_"$sampleId"_customPASS.bed | \
    tee "$seqId"_"$sampleId"_customClinicalCoverageTargetMetrics.txt | \
    awk '{pass[$4]+=$6; len[$4]+=$7} END { for(i in pass) printf "%s\t %.2f%\n", i, (pass[i]/len[i]) * 100 }' | \
    sort -k1,1 > "$seqId"_"$sampleId"_customClinicalCoverageGeneCoverage.txt

    echo "start R"

    /share/apps/R-distros/R-3.3.1/bin/Rscript --vanilla `dirname $0`/calculateTargetCoverage.R $seqId $sampleId

    rm "$seqId"_"$sampleId"_customClinicalCoverageGeneCoverage.txt
    rm "$seqId"_"$sampleId"_customClinicalCoverageTargetMetrics.txt

done
