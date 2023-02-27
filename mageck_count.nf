#!/usr/ein/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

///////////////////////////////////////////////////////////////////////////////
//// Stolen from Nourdine Bah (nourdine.bah@crick.ac.uk) and modified   //////
/////////////////////////////////////////////////////////////////////////////

MD_ANACONDA = "Anaconda2/2019.03"
CONDA_ENV = '/camp/stp/babs/working/schneid/conda/envs/RegenieGWA'
R_ENV = '/camp/stp/babs/working/software/anaconda/envs/R-3.6.1-rstudio-BABS'


PY_DIR = Paths.get(workflow.projectDir.toString(),"scripts").toString()
OUT_DIR = Paths.get( workflow.projectDir.toString() , "results", "counts" ).toString()

////////////////////////////////////////////////////////////////////////////////
////// SCRIPTS ////////////////////////////////////////////////////////////////

filter_bam_script = Channel.fromPath("scripts/filter_bam.cpp")

include { BAM_ALIGN_COUNT } from './workflows/align_bams.nf'
include { mageck_count_fastq; mageck_test; mageck_mle } from './modules/mageck.nf'

///////////////////////////////////////////////////////////////////////////////
//// SAMPLES //////////////////////////////////////////////////////////////////

Channel
	.fromPath(params.Samplesheet)
	.splitCsv(header: true)
	.map{
                it << [
                        "fastqs": sprintf("%s %s",
                                it["fastq_1"], it["fastq_2"])
                ]
        }
	.set{ SAMPLES }

///////////////////////////////////////////////////////////////////////////////
//// POOLS ////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Channel
	.fromPath("input/sgRNA_lib/*_library.csv")
	.map{[
		it.toString().replaceAll("(.*)/sgRNA_lib/(.*)_library.csv", "\$2"),
		it
	]}
	.set{ POOL }


///////////////////////////////////////////////////////////////////////////////
//// DESIGN MATRICES //////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


Channel
        .fromPath("input/design_mats/design_*.txt")
        .map{[
                'Brie',
		it.toString().replaceAll("(.*)/design_mats/design_(.*).txt", "\$2"),
                it
        ]}
        .set{ Design }


///////////////////////////////////////////////////////////////////////////////
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

	container "/camp/stp/babs/working/schneid/projects/tybulewiczv/edina.schweighoffer/CRISPR_screen_Mice_Bcells/input/assests/sequencing/sequencing.sif"

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

//
workflow {
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

        reformat_pool
		.out
		.control
		.combine( count_fastq.out.counts , by: 0 )
		.set{ MLE_INPUT }
//        MLE_INPUT.view()

        //////////////// ALIGN TRIM and BAM COUNT /////////////////////////////////
	BAM_ALIGN_COUNT(reformat_pool.out,
		        unzip.out)
	/////////////////////////////////////////////////////////////////////////// 

	Design
		.combine( MLE_INPUT, by: 0 )
		.set{ MLE_INPUT_2 }
  
//	MLE_INPUT_2.flatMap{ listSpread(it) }
//		.set{ MLE_INPUT_3 }
	
	mageck_test(count_fastq.out.counts)
	mageck_mle(MLE_INPUT_2)	
}
