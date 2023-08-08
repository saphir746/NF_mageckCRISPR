#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

////////////////////////////////////////////////////////////////////////////////
////// SCRIPTS ////////////////////////////////////////////////////////////////

//filter_bam_script = Channel.fromPath("scripts/filter_bam.cpp")

include { BAM_ALIGN_COUNT } from './workflows/align_bams.nf'
include { mageck_count_fastq } from './modules/mageck.nf'

//// PROCESSES ////////////////////////////////////////////////////////////////

publish_mode = "symlink"
publish_overwrite = true

process merge_fastq {

        executor "local"
        tag { sample_name }

        publishDir Paths.get( params.outDIR , "samples" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                val DATA

        output:
                tuple val(DATA), path("${Sample_name}.fastq.gz")

        script:
                fastQ=DATA["fastqs"]
                sample_name=DATA['Sample_ID']
                """
                zcat $fastQ | gzip -c > "${sample_name}.fastq.gz"
                """
}

process unzip {

        executor "local"

        tag { sample_name }

        publishDir Paths.get( params.outDIR , "samples" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(DATA), path(fastQgz)

        output:
                tuple val(DATA), path("${Sample_name}.fastq")

        script:
                sample_name=DATA['Sample_ID']
                """
                zcat -c $fastQgz > ${sample_name}.fastq
                """
}

process reformat_pool {

        executor "local"

        container "docker://bahnk/sequencing:v2"

        tag { pool }

        publishDir Paths.get( params.outDIR , "pools" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(pool), val(csv)

        output:
                tuple val(pool), path("${pool}.fa"), emit: fasta
                tuple val(pool), path("${pool}.revcomp.fa"), emit: fasta_revcomp
                tuple val(pool), path("${pool}_library.txt"), emit: library
                tuple val(pool), path("${pool}_sequence.txt"), emit: sequence
                tuple val(pool), path("${pool}_control.txt"), emit: control
                tuple val(pool), path("${pool}_targetsOnly.txt"), emit: targets

        script:
                """
                cat "${csv}" \
                        | sed '1d' \
                        | awk -F "," '{ printf ">%s\\n%s\\n", \$1, \$2 }' \
                        > "${pool}.fa"

                seqtk seq -r "${pool}.fa" > "${pool}.revcomp.fa"

                cat "${csv}" \
                        | sed '1d' \
                        | awk -F "," '{ printf "%s\\t%s\\t%s\\n", \$1, \$2 , \$3 }' \
                        > "${pool}_library.txt"

                cat "${csv}" \
                        | sed '1d' \
                        | awk -F "," '{ print \$3 }' \
                        | sort \
                        > "${pool}_sequence.txt"

                cat "${csv}" \
                        | grep -i NoTarget \
                        | awk -F "," '{ print \$1 }' \
                        | sort \
                        | uniq \
                        > "${pool}_control.txt"

               cat "${csv}" \
                        | grep -v NoTarget \
                        | awk -F "," '{ printf "%s\\t%s\\t%s\\n", \$1, \$2 , \$3 }' \
                        > "${pool}_targetsOnly.txt"
                """
}

///////////////////////////////////////////////////////////////////////////////
//// MAIN WORKFLOW ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

workflow {
 
         Channel
                .fromPath(params.Samplesheet)
                .splitCsv(header: true)
                .set{ SAMPLES }
                
         Channel
                .fromPath(Paths.get(params.PoolPath, "*_library.csv"))
                .map{[
                        it.toString().replaceAll("(.*)/([a-zA-Z0-9_]+)_library.csv", "\$2"),
                        it
                        ]}
                .set{ POOL }
        
       ////////////////////////////////////////////////////////////////////////////
        merge_fastq(SAMPLES)
        unzip(merge_fastq.out)
        
        reformat_pool(POOL)
        unzip
                .out
                .map{ [ it[0]["Sample_name"] , it[1] ] }
                .combine(POOL)
                .map{ [ it[2], it[0], it[1] ] }
                .groupTuple()
                .set{ FASTQ }


        reformat_pool
                .out
                .library
                .combine( FASTQ , by: 0 )
                .set{ TO_COUNTS }

        mageck_count_fastq(TO_COUNTS)
        ////////////////////////////////////////////////////////////////////////////
        if (params.BAM_ALIGN) {
          
            BAM_ALIGN_COUNT(
              reformat_pool.out.fasta
              unzip.out
              )
        }
}


