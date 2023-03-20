#!/bin/sh

module load Nextflow/22.04.0
module load Singularity/3.6.4

WORK_DIR=/camp/stp/babs/scratch/schneid/PM18257_CRISPR_2/
PROJECT_dir=/camp/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/
SAMSHEET=${PROJECT_dir}"input/sample_sheets/Sample_sheet_PM18257_runs_combined_2023.csv"
POOOOLPATH=${PROJECT_dir}"input/sgRNA_lib_mod"
DESIGNSPATH=${PROJECT_dir}"input/design_mats"

export NXF_SINGULARITY_CACHEDIR=${PROJECT_dir}/images/cachedir/
##

nextflow run mageck_count.nf \
			--Samplesheet ${SAMSHEET} \
			--PoolPath ${POOOOLPATH} \
			--designsMLE ${DESIGNSPATH} \
			--BAM_ALIGN FALSE \
			--outDIR ${PROJECT_dir}"results_mod/" \
                        -params-file params.yml \
                        -work-dir $WORK_DIR \
                        -resume

