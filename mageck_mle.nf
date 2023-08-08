#!/usr/ein/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

///////////////////////////////////////////////////////////////////////////////
/////////////////// Heavily modified  ///////////////////////////////// //////
/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////// SCRIPTS ////////////////////////////////////////////////////////////////

include { mageck_count_fastq; mageck_test; mageck_mle } from './modules/mageck.nf'

//// PROCESSES ////////////////////////////////////////////////////////////////

publish_mode = "symlink"
publish_overwrite = true


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
	        .fromPath(Paths.get(params.Designs,"*.txt"))
	        .map{[
	                it.toString().replaceAll("(.*)/design_(.*).txt", "\$2"),
	                it
	        ]}
	        .set{ Design }

	Channel
	        .fromPath(Paths.get(params.Pool,"*_control.txt"))
//	        .map{[
//	                it.toString().replaceAll("(.*)/Brie_(.*).txt", "\$2"),
//	                it             
//	        ]}
	        .set{ reformat_pool }
    
	Channel
	        .fromPath(params.COUNTS)
	  	.set{ count_fastqs }

	reformat_pool
//	        .control
	        .combine( count_fastqs )
	        .set{ MLE_INPUT }
	
	Design
	      .combine( MLE_INPUT )
	      .set{ MLE_INPUT_2 }
	MLE_INPUT_2.view()	  
//	MLE_INPUT_2.flatMap{ listSpread(it) }
//	    .set{ MLE_INPUT_3 }
	      
	mageck_mle(MLE_INPUT_2)
}
