#!/bin/sh

module load Nextflow/20.12.0-edge
module load Singularity/3.6.4

WORK_DIR=/camp/stp/babs/scratch/schneid/PM18257_CRISPR_2/
PROJECT_dir=/camp/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/
SAMSHEET=${PROJECT_dir}"input/sample_sheets/Sample_sheet_PM18257_runs_combined.csv"
POOOOLPATH=${PROJECT_dir}"input/sgRNA_lib"
DESIGNSPATH=${PROJECT_dir}"input/design_mats"


nextflow run mageck_count.nf \
			--Samplesheet ${SAMSHEET} \
			--PoolPath ${POOOOLPATH} \
			--designsMLE ${DESIGNSPATH}  \
			--outDIR ${PROJECT_dir}"results/" \
                        -params-file params.yml \
                        -work-dir $WORK_DIR \
                        -resume

