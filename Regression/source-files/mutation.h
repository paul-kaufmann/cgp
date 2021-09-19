/*
 * mutation.h
 *
 *  Created on: Mar 23, 2018
 *      Author: paul
 */

#ifndef MUTATION_H_
#define MUTATION_H_

#include "cgp.h"
#include "optimizationlog.h"

enum mutation_type_def {
	FUNCTION_MUTATION, CONNECTION_MUTATION, OUTPUT_MUTATION
};

extern const char *mutation_type_def2str[];

int get_num_mutant(int num_genes, double per_cent_mutate);
//int mutate_a_gene(int *chromosome, int *active_nodes, int *mutation_type);
//void mutate_a_chromosome(int *chromosome, int num_mutations);
//void get_used_nodes(int* chromosome, int* node_used);
void mutate_active_chromosome(int *chromosome, int num_mutations, int evaluation_no, class optimization_log *log);
//void randomize_inactive_nodes(int *chromosome, int *active_nodes);

void mutate(int *chromosome, int evaluation_no, optimization_log *log);

#endif /* MUTATION_H_ */
