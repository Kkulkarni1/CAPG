#include <stdlib.h>

#include "simulate_options.h"
#include "error.h"
#include "initialize.h"
#include "simulate.h"


int make_simulate_options(simulate_options **simo_in, initialize_options *inio, model_options *modo){
	
	simulate_options *so;
	*simo_in = malloc(sizeof **simo_in);
	
	if (!*simo_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"simulate_options");

	so = *simo_in;
	so->inio = inio;
	so->modo = modo;
	
	so->fastq_outfile = NULL;
	so->true_K = 0;
	so->cluster_spread = 0;
	so->haplotype_spread = 0;
	so->hap_control_file = NULL;
	so->ancestor = NULL;
	so->background_model = 0;
	
	so->simulation_outfile = NULL;
	so->simulation_infile = NULL;
	so->quality_change = 0;
	so->sim_n_reads = 0;
	so->sim_read_length = 0;
	//	op->simulation_method = SIMULATE_FROM_MODEL_WITH_DATA;	/* [KSD, TMP] */
	so->simulation_initialization = INIT_PARTITION;
	so->simulation_control = CONTROL_NOTHING; //CONTROL_ABUNDANCE, EXTRA_HAPLOTYPE
	
	so->simulation_random = MIXING | NUCLEOTIDES | HAPLOTYPES;
	// [KSD, TMP]		| QUALITIES;	/* all of it */
	so->simulation_from_data = QUALITIES;
	// [KSD, TMP] | NUCLEOTIDES | PARTITION | HAPLOTYPES;
	so->sim_pi_alpha = 100.0;	/* uniform */
	so->sim_gamma_alpha = 100.0;	/* uniform */

	return NO_ERROR;
}/* make_simulate_options */
