/*
 * optimization_algorithms.h
 *
 *  Created on: Mar 23, 2018
 *      Author: paul
 */

#ifndef OPTIMIZATION_ALGORITHMS_H_
#define OPTIMIZATION_ALGORITHMS_H_

/*
 * optimization_algorithms.c
 *
 *  Created on: Mar 23, 2018
 *      Author: paul
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include <getopt.h>
#include <float.h>
#include "cgp.h"
#include "mutation.h"
#include "optimization_algorithms.h"
#include "statistical_methods.h"

// #include <values.h>

// algorithm selection
enum algorithm_type_def {
	LAHC, SIMULATED_ANNEALING, RANDOM_WALK, PLUS_ES, HC, VND, EA, SSEA
};

extern enum algorithm_type_def algorithm_type;

void generate_new_pop_es(int** chromosomes, int* best_chromosome);
double EA_obsolete(int *gen_of_best, int *generations_to_termination,
		int* num_nodes_active, int run, char *prog);
double plus_es(int *evals_to_perfect, int *num_nodes_active,
		long *nodes_processed);
double SA(int *gen_of_best, int *num_nodes_active, long *nodes_processed);
double run_optimizer(int num_runs_total);
double lahc(int *evals_to_perfect, int *num_nodes_active,
		long *nodes_processed);
double hc(int *evals_to_perfect, int *num_nodes_active, long *nodes_processed);

#endif /* OPTIMIZATION_ALGORITHMS_H_ */
