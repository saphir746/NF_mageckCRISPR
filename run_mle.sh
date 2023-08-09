#!/bin/sh

module load Nextflow/22.04.0
module load Singularity/3.6.4


#/nemo/stp/babs/scratch/schneid/PM18257_CRISPR_2/
PROJECT_dir=/nemo/stp/babs/working/schneid/projects/vousdenk/julia.weber/CRISPR_screen_for_metabolic_modulators_of_cellular_response_to_serine_glycine_starvation-jw392/
WORK_DIR=${PWD}"/work"
POOLPATH=${PROJECT_dir}"results_count/pools"
DESIGNSPATH=${PROJECT_dir}"input/design_mat"
COUNTSFILE=${PROJECT_dir}"results_count/counts/Metabolic.count.txt"
SamSheet=${PROJECT_dir}"input/sample_sheets/Sample_sheet_PM22319_runs_combined.csv "

export NXF_SINGULARITY_CACHEDIR=${PROJECT_dir}/images/cachedir/
##

nextflow run mageck_mle.nf \
              		--Pool ${POOLPATH} \
			--Designs ${DESIGNSPATH} \
			--COUNTS ${COUNTSFILE} \
                        --outDIR ${PROJECT_dir}"results_mle/" \
                        -params-file params.yml \
                        -work-dir $WORK_DIR \
                        -resume

