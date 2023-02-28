#!/usr/ein/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

///////////////////////////////////////////////////////////////////////////////
/////////////////// Heavily modified  ///////////////////////////////// //////
/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////// SCRIPTS ////////////////////////////////////////////////////////////////

//filter_bam_script = Channel.fromPath("scripts/filter_bam.cpp")

include { BAM_ALIGN_COUNT } from './workflows/align_bams.nf'
include { mageck_count_fastq; mageck_test; mageck_mle } from './modules/mageck.nf'

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
		tuple val(DATA), path("${sample_name}.fastq.gz")

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
		tuple val(DATA), path("${sample_name}.fastq")
	
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
			| grep -i nonTarget \
                        | awk -F "," '{ print \$1 }' \
                        | sort \
                        | uniq \
                        > "${pool}_control.txt"

               cat "${csv}" \
                        | grep -v nonTarget \
                        | awk -F "," '{ printf "%s\\t%s\\t%s\\n", \$1, \$2 , \$3 }' \
                        > "${pool}_targetsOnly.txt"
		"""
}

////////////////////////////////////////////////////////////////////

def listSpread(item) {
   def array = []
   item[1].each{ array.add( [ item[0] , it , item[2] , item[3] ] ) }
   return(array)
}

///////////////////////////////////////////////////////////////////////////////
//// MAIN WORKFLOW ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

workflow {

	Channel
        	.fromPath(params.Samplesheet)
        	.splitCsv(header: true)
        	.map{
                it << [
                        "fastqs": sprintf("%s %s",
                                it["fastq_1"], it["fastq_2"])
                	]}
        	.set{ SAMPLES }

	Channel
        	.fromPath(Paths.get(params.PoolPath, "*_library.csv"))
        	.map{[
                	it.toString().replaceAll("(.*)/sgRNA_lib/(.*)_library.csv", "\$2"),
                	it
        		]}
        	.set{ POOL }


	Channel
        	.fromPath(Paths.get(params.designsMLE, "design_*.txt"))
        	.map{[
                	'Brie',
                	it.toString().replaceAll("(.*)/design_mats/design_(.*).txt", "\$2"),
                	it
        		]}
        	.set{ Design }
        
	////////////////////////////////////////////////////////////////////////////
	merge_fastq(SAMPLES)
	unzip(merge_fastq.out)	
	
	reformat_pool(POOL)
        unzip
		.out
		.map{ [ "Brie" , it[0]["Sample_name"] , it[1] ] }
		.groupTuple()
		.set{ FASTQ }

       	reformat_pool
		.out
		.library
		.join(FASTQ)
		.set{ TO_COUNTS }
	
	mageck_count_fastq(TO_COUNTS)

//        reformat_pool
//		.out
//		.control
//		.combine( count_fastq.out.counts , by: 0 )
//		.set{ MLE_INPUT }
////        MLE_INPUT.view()
//
//        //////////////// ALIGN TRIM and BAM COUNT /////////////////////////////////
//	if( params.BAM_ALIGN )
//		BAM_ALIGN_COUNT(reformat_pool.out,
//			        unzip.out)
//			.set{ MLE_INPUT }
//        else
//		mageck_count_fastq(TO_COUNTS)
//
//                reformat_pool
//                	.out
//                	.control
//                	.combine( count_fastq.out.counts , by: 0 )
//                	.set{ MLE_INPUT }
//	/////////////////////////////////////////////////////////////////////////// 
//
//	Design
//		.combine( MLE_INPUT, by: 0 )
//		.set{ MLE_INPUT_2 }
//	
//	mageck_test(count_fastq.out.counts)
//	mageck_mle(MLE_INPUT_2)	
}
