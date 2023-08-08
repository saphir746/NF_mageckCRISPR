#!/bin/sh

module load Nextflow/22.04.0
module load Singularity/3.6.4


#/nemo/stp/babs/scratch/schneid/PM18257_CRISPR_2/
PROJECT_dir=/nemo/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/
WORK_DIR=${PWD}"/werk"
POOLPATH=${PROJECT_dir}"results_reSeqed_bis/pools"
DESIGNSPATH=${PROJECT_dir}"input/design_mats_bis"
COUNTSFILE=${PROJECT_dir}"results_reSeqed_bis/counts/Brie_modified.count.txt"
SamSheet=${PROJECT_dir}"input/sample_sheets/Sample_sheet_PM18257_runs_combined_June2023.csv"

export NXF_SINGULARITY_CACHEDIR=${PROJECT_dir}/images/cachedir/
##

nextflow run mageck_mle.nf \
              		--Pool ${POOLPATH} \
			--Designs ${DESIGNSPATH} \
			--COUNTS ${COUNTSFILE} \
                        --outDIR ${PROJECT_dir}"results_reSeqed_bis/" \
                        -params-file params.yml \
                        -work-dir $WORK_DIR \
                        -resume

