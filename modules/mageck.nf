import java.nio.file.Paths

publish_mode = "symlink"
publish_overwrite = true

process mageck_count_fastq {

        executor "slurm"
        cpus 2
        time "1h"
        memory "5G"

        container "docker://davidliwei/mageck:latest"

        tag { pool }

        publishDir Paths.get( params.outDIR , "counts" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(pool), path(library), val(samples), path(fastqs)

        output:
                tuple val(pool), path("${pool}.count.txt"), emit: counts
                tuple val(pool), path("${pool}.count_normalized.txt"), emit: norm_counts
                tuple val(pool), path("${pool}.countsummary.txt"), emit: summary

        script:
                samples_all = samples.join(",")
                """
                mageck count \
                        --list-seq $library \
                        --fastq $fastqs \
                        --output-prefix $pool \
                        --sample-label $samples_all
                """
}


process mageck_test {

        executor "slurm"
        cpus 4
        time "3h"
        memory "10G"

        container "docker://davidliwei/mageck:latest"

        tag { postfix }

        publishDir Paths.get( params.outDIR , "univariate" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
               tuple val(pool), path(counts)

        output:
                path("${postfix}.gene_summary.txt"), emit: gene_summary
                path("${postfix}.sgrna_summary.txt"), emit: sgrna_summary
                path("${postfix}.log"), emit: log

        script:
                postfix="EV_vs_FRC"
                """
                mageck test \
                        -k $counts \
                        -t A6 \
                        -c A7,A3 \
                        -n $postfix
                """
}

process mageck_mle {

        executor "slurm"
        cpus 4
        time "3h"
        memory "10G"

        container "docker://davidliwei/mageck:latest"

        tag { postfix }

        publishDir Paths.get( params.outDIR , "mle" ),
                mode: publish_mode,
                overwrite: publish_overwrite

        input:
                tuple val(pool), val(postfix), path(design), path(control), path(counts)

        output:
                path("${postfix}.gene_summary.txt"), emit: gene_summary
                path("${postfix}.sgrna_summary.txt"), emit: sgrna_summary
                path("${postfix}.log"), emit: log

        script:

                """
                mageck mle \
                        --design-matrix $design \
                        --count-table $counts \
                        --norm-method control \
                        --update-efficiency \
                        --control-sgrna $control \
                        --output-prefix $postfix
                """
}

