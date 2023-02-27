#!/usr/en/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

////////////////////// alignments ///////////////////////////////////


process bowtie2_index {

        executor "local"
        container "/camp/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/input/assests/sequencing/sequencing.sif"

        tag { pool }

        publishDir Paths.get( "results" , "pools" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(pool), path(fasta)

        output:
                tuple val(pool), path("pool${pool}*.bt2")

        script:

                name = "pool${pool}"

                """
                bowtie2-build $fasta $name
                """
}


process trim {

        executor "slurm"
        container "/camp/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/input/assests/sequencing/sequencing.sif"

        cpus 6
        maxForks 3

        tag { name }

        publishDir Paths.get( "results" , "samples" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(DATA), path(fastq)

        output:
                tuple val(DATA), path("${name}.trimmed.fastq")

        script:

                name = DATA["Sample_ID"]

                """
                cutadapt -g TTGTGGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTA $fastq > "${name}.trimmed.fastq"
                """
}

//                cutadapt -g TTGTGGAAAGGACGAAACACCG $fastq > tmp.fastq
//                cutadapt -a GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACT tmp.fastq > "${name}.trimmed.fastq"




process bowtie2 {

        executor "slurm"
        container "/camp/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/input/assests/sequencing/sequencing.sif"
        cpus 12
        time "01:00:00"

        tag { name }

        publishDir Paths.get( "results" , "bowtie2" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(metadata), path(fastq), val(lib_pool), path(index)

        output:
                tuple val(metadata), val(lib_pool), path("${name}.bam")

        script:

                sample = metadata["Sample_ID"]
                name = "${sample}_${lib_pool}"

                """
                bowtie2 \
                        --threads $task.cpus \
                        -x "pool${lib_pool}" \
                        -U $fastq \
                        --norc \
                        | samtools view -h -F 4 -x XS - \
                        | samtools view -h -bS - > "${name}.bam"
                """
}

process filter_bam {

        executor "local"
        container "/camp/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/input/assests/sequencing/sequencing.sif"

        cpus 2
        time "01:00:00"

        tag { sample }

        publishDir Paths.get( "results" , "bowtie2" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(metadata), val(lib_pool), path(bam), path(script)

        output:
                tuple val(metadata), val(lib_pool), path("${sample}.filtered.bam")

        script:

                sample = metadata["Sample_ID"]

                """
                samtools view -F 4 -x XS $bam | samtools view -Sb - > "${sample}.filtered.bam"
                """
}

process count_bams {

        executor "slurm"
        cpus 2
        time "1h"
        memory "5G"
        //executor "local"
        container "/camp/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/input/assests/mageck/mageck.sif"

        tag { pool }

        publishDir Paths.get( "results" , "counts_bowtie2" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(pool), path(library), val(samples), path(bams)

        output:
                tuple val(pool), path("${pool}_aligned.count.txt"), emit: counts
                tuple val(pool), path("${pool}_aligned.count_normalized.txt"), emit: norm_counts
                tuple val(pool), path("${pool}_aligned.countsummary.txt"), emit: summary
        script:

                labels = samples.join(",")
                prefix = "${pool}_aligned"

                """
                mageck count \
                        --list-seq $library \
                        --fastq $bams \
                        --output-prefix $prefix \
                        --sample-label $labels
                """
}

////////////////////////////////////////////////////////////////////

workflow BAM_ALIGN_COUNT {
   take:
		fasta_pools
		unzipped
   main:
       
       // bowtie2_index( reformat_pool.out.fasta )
       bowtie2_index(fasta_pools)
	// trim(unzip.out)
       trim(unzipped)
	trim
                .out
                .combine(bowtie2_index.out)
                .set{ TO_ALIGN }

        bowtie2(TO_ALIGN)

//      filter_bam( bowtie2.out.combine(filter_bam_script) )

        bowtie2
                .out
                .map{ [ it[1], it[0]["Sample_name"], it[2] ] }
                .groupTuple()
                .combine(fasta_pools.library)
                .map{ [ it[0], it[4], it[1], it[2] ] }
                .set{ BAMS_TO_COUNT }

        count_bams(BAMS_TO_COUNT)

        fasta_pools
          .control
          .combine( count_bams.out.counts , by: 0 )
          .set{ MLE_INPUT_BAMS }

    emit:
	  MLE_INPUT_BAMS

}
