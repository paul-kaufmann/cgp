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

extern const char* mutation_type_def2str[];

int get_num_mutant(const int num_genes, double per_cent_mutate);
//int mutate_a_gene(int* chromosome, int *active_nodes, int* mutation_type);
//void mutate_a_chromosome(int *chromosome, int num_mutations);
void get_used_nodes(const int* chromosome, int* node_used);
void convert_active_nodes_to_genes(const int *active_nodes, int* active_genes);
void get_used_genes(const int* chromosome, int* active_nodes, int* active_genes);
void get_used_genes(const int* chromosome, int* active_genes);
void mutate_active_genes(int* chromosome, int num_mutations, long long evaluation_no, optimization_log *log);
void randomize_inactive_nodes(int* chromosome, const int *active_nodes);

void mutate(int* chromosome, int evaluation_no, optimization_log *log);

#endif /* MUTATION_H_ */
