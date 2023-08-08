#!/bin/sh

module load Nextflow/20.12.0-edge
module load Singularity/3.6.4

WORK_DIR=/flask/scratch/babs/schneid/PM22319_CRISPR/
PROJECT_dir=/nemo/stp/babs/working/schneid/projects/vousdenk/julia.weber/CRISPR_screen_for_metabolic_modulators_of_cellular_response_to_serine_glycine_starvation-jw392/
POOOOLPATH=${PROJECT_dir}"input/sgRNA_lib"

export NXF_SINGULARITY_CACHEDIR=${PROJECT_dir}/images/cachedir/
##

nextflow run mageck_count.nf \
			--Samplesheet ${PROJECT_dir}"input/sample_sheets/Sample_sheet_PM22319_runs_combined.csv" \
			--PoolPath ${POOOOLPATH} \
			--BAM_ALIGN FALSE \
			--outDIR ${PROJECT_dir}"results_count/" \
                        -params-file params.yml \
                        -work-dir $WORK_DIR \
                        -resume

