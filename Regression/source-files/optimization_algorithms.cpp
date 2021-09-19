/*
 * optimization_algorithms.c
 *
 *  Created on: Mar 23, 2018
 *      Author: paul
 */

#include "optimization_algorithms.h"

#include <iomanip>
#include <iostream>

#include "optimizationlog.h"

enum algorithm_type_def algorithm_type;

//struct sa_parameters sa_parameters_global;

// *************************************************************************************************************************************************
// Miller's implementation --- Miller's implementation --- Miller's implementation --- Miller's implementation
// *************************************************************************************************************************************************

/* (1+lambda evolutionary strategy where lambda = population size -1 */
void generate_new_pop_es(int **chromosomes, int *best_chromosome) {
	int j, k;
	int num_mutant;
	optimization_log l;

	num_mutant = get_num_mutant(num_genes, per_cent_mutate);

	/* copy best_chromosome into last member of chromosome array */
	for (j = 0; j < num_genes; j++)
		chromosomes[population_size - 1][j] = best_chromosome[j];

	/* generate new population by mutating all but last */

	for (k = 0; k < population_size - 1; k++) {
		for (j = 0; j < num_genes; j++) /* copy best chromosome */
			chromosomes[k][j] = best_chromosome[j];

		/* mutate the chromosome */
		//mutate_a_chromosome(chromosomes[k], num_mutant);
		mutate_active_chromosome(chromosomes[k], num_mutant, -1, &l);
	}
}

///* Do a run of the EA, Millers original implementation */
//// this implementation should not be used anymore, please use the general implementation of mu plus/comma lambda EA optimization scheme
//double EA_obsolete(int *gen_of_best, int *generations_to_termination, int* num_nodes_active, int run, char *prog) {
//	int gen, best_gen = 0;
//	int generations = 0;
//	int nodes_to_process[MAX_NUM_NODES];
//	int** chromosomes;
//	int* best_chromosome;
//	double previous_best_fit = 99999999999.0;
//	double best_fit = previous_best_fit;
//
//	chromosomes = create_2d_space(num_genes, population_size);
//	best_chromosome = create_1d_space(num_genes);
//
//	initialise(chromosomes);
//
//	for (gen = 1; gen <= num_generations; gen++) {
//		/*write_generation_to_screen(gen);*/
//		best_fit = get_best_chromosome(chromosomes, best_chromosome, previous_best_fit, gen);
//
//		check_if_improvement(best_fit, &previous_best_fit, &best_gen, gen, prog, best_chromosome);
//
//		/* jump out of run if maximum fitness acheived */
//		if ((best_fit <= perfect) && (shrink_phenotype == 0)) {
//			generations = gen;
//			break;
//		} else {
//			/* create a new population */
//			generate_new_pop_es(chromosomes, best_chromosome);
//		}
//
//	}
//
//	if (best_fit > perfect) {
//		generations = num_generations;
//	}
//
//	*num_nodes_active = get_nodes_to_process(best_chromosome, nodes_to_process);
//
//	/*write_result_of_EA_to_file(run, best_gen, best_fit, best_chromosome);*/
//
//	*gen_of_best = best_gen;
//	*generations_to_termination = generations;
//
//	/* write the raw best chromosome to cgp.chr */
//	/*fprint_a_raw_chromosome(best_chromosome,"cgp.chr",0);*/
//
//	/*printf("%lf", best_fit);*/
//
//	free_array2d(population_size, chromosomes);
//	free(best_chromosome);
//
//	return best_fit;
//}

// *************************************************************************************************************************************************

// Proper mu + lambda Local Search implementation -- Proper mu + lambda Local Search implementation -- Proper mu + lambda Local Search implementation

// *************************************************************************************************************************************************

/* mu+lambda ES implementation */
double plus_es(long long *evals_to_perfect, int *num_nodes_active, long *nodes_processed) {

	// reset variables
	*evals_to_perfect = INT_MAX;
	*num_nodes_active = -1;
	*nodes_processed = 0;

	int lambda = population_size;

	int num_mutant = get_num_mutant(num_genes, per_cent_mutate);
	printf("num mutated genes: %d\n", num_mutant);

	if (lambda < 1 || mu < 1)
		return INT_MAX;

	if (num_mutant < 1)
		return INT_MAX;

	optimization_log l;
	int nodes_to_process[MAX_NUM_NODES];
	int **population;
	int *population_tmp[mu + lambda];
	double my_fitness[mu + lambda];
	double my_fitness_tmp[mu + lambda];
	population = create_2d_space(num_genes, mu + lambda);
	int *best;
	best = create_1d_space(num_genes);
	double best_fitness = -FLT_MAX;
	int i, j;
	int number_of_evaluations = 1;

	// initialize parents
	for (i = 0; i < mu; i++)
		generate_a_random_chromosome(population[i]);

	// evaluate fitness of parent individuals
	for (i = 0; i < mu; i++) {
		my_fitness[i] = fitness(population[i], nodes_processed);
		number_of_evaluations++;
	}

	do {
		// create children
		for (i = mu; i < mu + lambda; i++) {
			int random_individual = newrand(mu);
			memcpy(population[i], population[random_individual], num_genes * sizeof(int));
			//mutate_a_chromosome(population[i], num_mutant);
			//mutate_active_chromosome(population[i], num_mutant, number_of_evaluations, &l);
			mutate(population[i], number_of_evaluations, &l);
//			my_fitness[i] = fitness2(population[i], nodes_processed, test_data_set);
			my_fitness[i] = fitness(population[i], nodes_processed);
			number_of_evaluations++;
		}

		// select new parents, respect neutrality i.e. parent propagate to next population only if strictly better than children
		// population contains parents (0..(mu-1)) and then children (mu..(mu+lambda-1))
		// strategy: start from mu+lambda-1 to 0, find FIRST best individual
		// copy its pointer and fitness to population_tmp and fitness_tmp
		// remove pointer to this BEST FIRST individual from population
		// at the end: population_tmp and fitness_tmp contain from 0 to mu+lambda-1 pointers to individuals
		// with rising costs; children come first then parents that are on par
		for (i = 0; i < mu + lambda; i++) {

			// find best individual, children are favored if on par or better than parents
			int tmp = -1;
			for (j = mu + lambda - 1; j >= 0; j--) {

				// ignore individuals already copied to new population[]
				if (population[j] == NULL)
					continue;

				// remember the index of last not NULL element in population
				if (tmp < 0)
					tmp = j;
				// else if costs of an individual with lower index is
				else if (my_fitness[j] > my_fitness[tmp])
					tmp = j;
			}
			population_tmp[i] = population[tmp];
			my_fitness_tmp[i] = my_fitness[tmp];
			population[tmp] = NULL;
		}

		memcpy(population, population_tmp, sizeof(int*) * (mu + lambda));
		memcpy(my_fitness, my_fitness_tmp, sizeof(double) * (mu + lambda));

		// best individual ever found?
		if (my_fitness[0] > best_fitness) {
			memcpy(best, population[0], num_genes * sizeof(int));
			best_fitness = my_fitness[0];
			if (progress_report)
				write_progress_info_to_screen(number_of_evaluations, *nodes_processed, best_fitness);
		}

	} while (best_fitness < perfect && number_of_evaluations < max_number_evaluations);

	if (best_fitness >= perfect)
		*evals_to_perfect = number_of_evaluations;

	*num_nodes_active = get_nodes_to_process(best, nodes_to_process);

	double best_test_fitness = fitness2(best, nodes_processed, test_data_set);
	printf("\ntraining fitness %f  test fitness %f\n", best_fitness, best_test_fitness);

	free_array2d(mu + lambda, population);
	free(best);
	return best_test_fitness;
}

int tournament_selection_best(int **population, double *fitness, int pop_size, int t) {

	int candidate, winner = newrand(pop_size);
	while (t-- > 0) {
		candidate = newrand(pop_size);
		if (fitness[candidate] > fitness[winner])
			winner = candidate;
	}
	return winner;
}

int tournament_selection_worst(int **population, double *fitness, int pop_size, int t) {

	int candidate, winner = newrand(pop_size);
	while (t-- > 0) {
		candidate = newrand(pop_size);
		if (fitness[candidate] < fitness[winner])
			winner = candidate;
	}
	return winner;
}

/**
 * regular population-based EA with tournament selection scheme
 */
double ea(long long *evals_to_perfect, int *num_nodes_active, long *nodes_processed) {

	// reset variables
	*evals_to_perfect = INT_MAX;
	*num_nodes_active = -1;
	*nodes_processed = 0;

	if (mu < 1)
		return INT_MAX;

	optimization_log *l = NULL;
	//l = new (optimization_log);

	int nodes_to_process[MAX_NUM_NODES];
	int **population;
	double my_fitness[2 * mu];
	double best_float;

	population = create_2d_space(num_genes, 2 * mu);
	int i, j, best, number_of_evaluations = 1;

	// initialize parents
	for (i = 0; i < mu; i++)
		generate_a_random_chromosome(population[i]);

	// evaluate fitness of parent individuals
	for (i = 0; i < mu; i++) {
		my_fitness[i] = fitness(population[i], nodes_processed);
		number_of_evaluations++;
	}

	double old_best_fitness = my_fitness[0];
	do {
		// elitism, copy best to new population
		memcpy(population[mu], population[0], num_genes * sizeof(int));
		my_fitness[mu] = my_fitness[0];
		best = mu;
		best_float = my_fitness[0];

		// create children
		for (i = mu + 1; i < 2 * mu; i++) {
			int random_individual = tournament_selection_best(population, my_fitness, mu, 10);
			memcpy(population[i], population[random_individual], num_genes * sizeof(int));
			mutate(population[i], number_of_evaluations, l);
			my_fitness[i] = fitness(population[i], nodes_processed);
			number_of_evaluations++;

			if (my_fitness[i] > my_fitness[best])
				best = i;
		}

		if (my_fitness[mu] != my_fitness[best])
			//if (progress_report)
			write_progress_info_to_screen(number_of_evaluations, *nodes_processed, my_fitness[best]);

		// copy best to position zero
		memcpy(population[0], population[best], num_genes * sizeof(int));
		my_fitness[0] = my_fitness[best];

		// copy rest of the new population and replace the old population
		i = 1;
		j = mu;
		while (j < 2 * mu) {
			if (j != best) {
				memcpy(population[i], population[j], num_genes * sizeof(int));
				my_fitness[i] = my_fitness[j];
				i++;
			}
			j++;
		}

		if (best_float > my_fitness[0]) {
			printf("best fitness reduced from %f to %f\n", best_float, my_fitness[0]);
			exit(1);
		}
	} while (my_fitness[0] < perfect && number_of_evaluations < max_number_evaluations);

	if (my_fitness[0] == perfect)
		*evals_to_perfect = number_of_evaluations;

	if (l) {
		l->print();
		delete l;
	}

	double best_test_fitness = fitness2(population[0], nodes_processed, test_data_set);
	printf("\ntraining fitness %f  test fitness %f\n", my_fitness[0], best_test_fitness);

	*num_nodes_active = get_nodes_to_process(population[0], nodes_to_process);
	free_array2d(2 * mu, population);
	printf("\nLeaving EA\n");
	return best_test_fitness;
}

/**
 * Hill Climber
 */
double hc(long long *evals_to_perfect, int *num_nodes_active, long *nodes_processed) {

	int nodes_to_process[MAX_NUM_NODES];

	optimization_log l;
	int *parent = create_1d_space(num_genes);
	generate_a_random_chromosome(parent);
	*nodes_processed = 0;
	double costs = fitness(parent, nodes_processed);

	int *child = create_1d_space(num_genes);
	double child_costs;

	int I = 0;
	while (costs < perfect && I < max_number_evaluations) {

		memcpy(child, parent, num_genes * sizeof(int));
		mutate_active_chromosome(child, 1, I, &l);
		child_costs = fitness(child, nodes_processed);

		// set fitness change type in the log
		if (child_costs < costs) // fitness worsening
			*l.getLastImpact() = -1;
		else if (child_costs > costs) // functional quality improving
			*l.getLastImpact() = 1;
		else
			*l.getLastImpact() = 0; // functional quality same as parent's

		if (child_costs >= costs) {
			memcpy(parent, child, num_genes * sizeof(int));
			if (child_costs > costs && progress_report)
				write_progress_info_to_screen(I, *nodes_processed, child_costs);
			costs = child_costs;
		}

		I++;
	}

	*evals_to_perfect = (costs == perfect) ? I : INT_MAX;
	*num_nodes_active = get_nodes_to_process(parent, nodes_to_process);

	free(parent);
	free(child);

	l.print();
	return costs;

}


// *************************************************************************************************************************************************

// General optimizer calling specific optimizers num_runs_total amount of times

// *************************************************************************************************************************************************

/* do multiple runs of EA and write out results */
double run_optimizer(int num_runs_total) {

	int run;
	int active_nodes_in_run[num_runs_total];
	double fitnesses[num_runs_total];
	long long fitness_evals_to_perfect[num_runs_total];
	long nodes_processed[num_runs_total];

	for (run = 0; run < num_runs_total; run++) {
		printf("\n\nRUN %d\n", run);

		if (algorithm_type == PLUS_ES)
			fitnesses[run] = plus_es(&fitness_evals_to_perfect[run], &active_nodes_in_run[run], &nodes_processed[run]);
		else if (algorithm_type == HC)
			fitnesses[run] = hc(&fitness_evals_to_perfect[run], &active_nodes_in_run[run], &nodes_processed[run]);
		else if (algorithm_type == EA)
			fitnesses[run] = ea(&fitness_evals_to_perfect[run], &active_nodes_in_run[run], &nodes_processed[run]);
		else
			printf("no algorithm found: %d\n", algorithm_type);

	}

	printf("\n\nfitness values:   ");
	for (run = 0; run < num_runs_total; run++)
		printf("%f ", fitnesses[run]);
	printf("\n\nevals to perfect: ");
	for (run = 0; run < num_runs_total; run++)
		printf("%lld ", fitness_evals_to_perfect[run]);
	printf("\n\nnodes processed:   ");
	for (run = 0; run < num_runs_total; run++)
		printf("%ld ", nodes_processed[run]);
	printf("\n");

	report_final_results(stdout, num_runs_total, fitnesses, fitness_evals_to_perfect, active_nodes_in_run, nodes_processed);
	cout << endl;

	return 0;

	double ce_fitness;
	long long stop_at;
	ce_fitness = calc_CE(fitness_evals_to_perfect, num_runs_total, max_number_evaluations, &stop_at);
	cout << setw(41) << "computational effort [no. fit. evals.] " << ce_fitness << endl;
	cout << setw(41) << "restart at [no. fitness evals.] " << stop_at << endl << endl;

	double ce_nodes;
	long stop_at_long;
	ce_nodes = calc_CE(nodes_processed, num_runs_total, max(nodes_processed, num_runs_total), &stop_at_long);
	cout << setw(41) << "computational effort [no. node evals.] " << ce_nodes << endl;
	cout << setw(41) << "restart at [no. node evals.] " << stop_at_long << endl;

	double q1, retval, q3;
	long long q1_int, retval_int, q3_int;
	long q1_long, retval_long, q3_long;

	switch (fitness_type) {

	case MEDIAN_FITNESS:
		quartiles(fitnesses, num_runs_total, &q1, &retval, &q3);
		break;
	case CE_FITNESS:
		retval = ce_fitness;
		break;
	case CE_NODES_PROCESSED:
		retval = ce_nodes;
		break;
	case MEDIAN_EVALUATIONS:
		quartiles(fitness_evals_to_perfect, num_runs_total, &q1_int, &retval_int, &q3_int);
		retval = (double) retval_int;
		break;
	case NODES_PROCESSED:
		quartiles(nodes_processed, num_runs_total, &q1_long, &retval_long, &q3_long);
		retval = (double) retval_long;
		break;
	default:
		retval = 0.0;
	}
	return retval;
}

