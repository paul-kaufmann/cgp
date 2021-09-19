/*
 * statistical_methods.c
 *
 *  Created on: Mar 23, 2018
 *      Author: paul
 */

#include "statistical_methods.h"
#include "cgp.h"
#include "mutation.h"
#include "optimization_algorithms.h"

#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

enum fitness_type_def fitness_type;

double get_standard_deviation(int num_items, double average, double *items) {
	int i;
	double temp, st_dev;

	st_dev = 0.0;
	for (i = 0; i < num_items; i++) {
		temp = (items[i] - average);
		temp = temp * temp;
		st_dev = st_dev + temp;
	}

	st_dev = st_dev / (double) num_items;
	st_dev = sqrt(st_dev);

	return st_dev;
}

double get_standard_deviation_int(int num_items, double average, int *items) {
	int i;
	double temp, st_dev;

	st_dev = 0.0;
	for (i = 0; i < num_items; i++) {
		temp = (items[i] - average);
		temp = temp * temp;
		st_dev = st_dev + temp;
		if (items[i] == INT_MAX || items[i] == INT_MIN)
			return 0;
	}

	st_dev = st_dev / (double) num_items;
	st_dev = sqrt(st_dev);

	return st_dev;
}

void report_final_results(FILE* outfile, int num_runs_total, double *fitnesses, long long *fitness_evals_to_perfect, int *active_nodes_in_run, long *nodes_processed) {

	printf("\n");

	print_avg_quartiles(outfile, "absolute fitness/cost                    ", fitnesses, num_runs_total);

	double rel_fitness[num_runs_total];
	for (int i = 0; i < num_runs_total; i++)
		rel_fitness[i] = min(1.0, fitnesses[i] / perfect);
	print_avg_quartiles(outfile, "relative fitness/cost                    ", rel_fitness, num_runs_total);

	print_avg_quartiles(outfile, "evaluations until perfect solution found ", fitness_evals_to_perfect, num_runs_total);
	print_avg_quartiles(outfile, "nodes processed                          ", nodes_processed, num_runs_total);
	print_avg_quartiles(outfile, "active nodes                             ", active_nodes_in_run, num_runs_total);
	cout << setw(41) << "lowest obj. function" << setw(14) << min(fitnesses, num_runs_total) << endl;
	cout << setw(41) << "highest obj. function" << setw(14) << max(fitnesses, num_runs_total) << endl;
}

