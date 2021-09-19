/* cgp-functions.c
 Julian F. Miller (c) 2009
 version 1.1 of first public release on 20-July-2009
 Dept. of Electronics, University of York, UK

 IMPORTANT: For Boolean unsigned problems
 program outputs are arranged most significant on the left
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "cgp.h"
#include "mutation.h"
#include "optimization_algorithms.h"
#include "statistical_methods.h"
#include <iostream>
#include <fstream>
#include <sstream>

#ifndef FLT_MAX
#define FLT_MAX __FLT_MAX__
#endif				/*  */

/* write out parameter values in results file */
void write_cgp_info(FILE *fp, char *command, char *parameterfile, char *datafile) {
	int i;
	fprintf(fp, "The program is        		%s\n", command);
	fprintf(fp, "The parameter file is 		%s\n", parameterfile);
	fprintf(fp, "The data file is      		%s\n", datafile);
	fprintf(fp, "The algorithm is			%d\n", algorithm_type);
	fprintf(fp, "The fitness type is		%d\n", fitness_type);
	fprintf(fp, "population_size is    		%d\n", population_size);
	fprintf(fp, "mutation rate is      		%6.2lf\n", per_cent_mutate);
	fprintf(fp, "explicit neutrality is		%d\n", explicit_neutrality);
	fprintf(fp, "start temperature is  		%lf\n", sa_parameters_global.tstart);
	fprintf(fp, "stop temperature is   		%lf\n", sa_parameters_global.tstop);
	fprintf(fp, "cooling scheme is     		%d\n", sa_parameters_global.cooling_scheme);
	fprintf(fp, "steps in temp level   		%d\n", sa_parameters_global.steps_in_temp_level);
	fprintf(fp, "rejection factor is   		%lf\n", sa_parameters_global.rejection_factor);
	fprintf(fp, "max_number_evaluations is	%lld\n", max_number_evaluations);
	fprintf(fp, "num_runs is          		%d\n", num_runs_total);
	fprintf(fp, "num_rows is          		%d\n", num_rows);
	fprintf(fp, "num_cols is         		%d\n", num_cols);
	fprintf(fp, "levels_back is       		%d\n", levels_back);
	fprintf(fp, "progress report is   		%d\n", progress_report);
	fprintf(fp, "report interval is   		%d\n", report_interval);
	fprintf(fp, "global_seed is       		%u\n", global_seed);
	fprintf(fp, "save_best_chrom is   		%d\n", save_best_chrom);
	fprintf(fp, "run_from_chrom is    		%d\n", run_from_chrom);
	fprintf(fp, "shrink_phenotype is  		%d\n", shrink_phenotype);
	fprintf(fp, "num_finctions is     		%d\n", num_functions);
	fprintf(fp, "num_genes is        		%d\n", num_genes);
	fprintf(fp, "num_nodes is         		%d\n", num_nodes);
	fprintf(fp, "num_genes_per_node is		%d\n", num_genes_per_node);
	fprintf(fp, "num_inputs is				%d\n", num_inputs);
	fprintf(fp, "num_outputs is				%d\n", num_outputs);
	fprintf(fp, "num_tests is				%d\n", num_tests);
	fprintf(fp, "perfect is					%lf\n", perfect);
	fprintf(fp, "bit_width is 				%d\n", bit_width);
	fprintf(fp, "mutate_active_node_func    %d\n", mutate_active_node_func);
	fprintf(fp, "mutate_inactive_node_func  %d\n", mutate_inactive_node_func);
	fprintf(fp, "mutate_active_wires        %d\n", mutate_active_wires);
	fprintf(fp, "mutate_inactive_wires      %d\n", mutate_inactive_wires);

	for (i = 0; i < MAX_NUM_FUNCTIONS; i++) {
		fprintf(fp, "%d %s\n", number[i], node_types[i]);
	}
	fprintf(fp, "\nHere are the Results\n");
}

/*  returns a random integer between 0 and range-1 */
int newrand(int range) {
	int temp;
	temp = rand() % range;
	return (temp);
}

/* read an item of data from problem
 specification file. Format of data
 depends on defined data type
 */
data_type myfscanf(FILE *fp) {
	data_type datum_read;

#ifdef DATA_IS_UNSIGNED_INT
	if (fscanf(fp, "%lu", &datum_read) == EOF) {
		printf("reading data from problem specification file failed\n");
		exit(1);
	}
#endif				/*  */
#ifdef DATA_IS_INT
	if (fscanf(fp, "%i", &datum_read) == EOF) {
		printf("reading data from problem specification file failed\n");
		exit(1);
	}
#endif				/*  */
#ifdef DATA_IS_DOUBLE
	if (fscanf(fp, "%lf", &datum_read) == EOF) {
		printf("reading data from problem specification file failed\n");
		exit(1);
	}
#endif				/*  */
	return datum_read;
}

/* reads input and output data from a file  (e.g. compressed truth table .plu)
 defines number of inputs to be num_inputs+num_constant_inputs
 */
void read_data(const char datafile[MAX_NUM_LETTERS]) {
	int i, j;
	char dummy[MAX_NUM_LETTERS];
	FILE *fp;
	fp = fopen(datafile, "r");
	if (!fp) {
		puts("ERROR. Missing input data file (e.g. .plu (compressed Boolean) .dat");
		exit(1);
	} else {
		if (fscanf(fp, "%s %d", dummy, &num_inputs) == EOF) {
			puts("Error, can't read number of inputs");
			exit(1);
		}
		if (fscanf(fp, "%s %d", dummy, &num_outputs) == EOF) {
			puts("Error, can't read number of outputs");
			exit(1);
		}
		if (fscanf(fp, "%s %d", dummy, &num_tests) == EOF) {
			puts("Error, can't read number of tests");
			exit(1);
		}
		if (num_tests >= MAX_NUM_DATA) {
			printf("\nERROR. Too many test cases (in datafile)\n");
			exit(0);
		}

		for (i = 0; i < num_tests; i++) {

			for (j = 0; j < num_inputs; j++)
				data_inputs[i][j] = myfscanf(fp);

			for (j = 0; j < num_outputs; j++)
				data_outputs[i][j] = myfscanf(fp);
		}

		fclose(fp);
	}
	num_genes = num_genes_per_node * num_nodes + num_outputs;

#ifdef DATA_IS_UNSIGNED_INT
	num_constant_inputs = 0;

#endif				/*  */
#ifdef DATA_IS_INT
	num_constant_inputs = MAX_NUM_CONSTANT_INPUTS;

#endif				/*  */
#ifdef DATA_IS_DOUBLE
	num_constant_inputs = MAX_NUM_CONSTANT_INPUTS;

#endif				/*  */
	num_inputs = num_inputs + num_constant_inputs;

	printf("ni=%d no=%d ntest=%d ngenes=%d ngenespernode=%d\n", num_inputs, num_outputs, num_tests, num_genes, num_genes_per_node);
}

/**
 * reads a data file and save the parameter in a private class
 */
data_set* read_data2(const char datafile[MAX_NUM_LETTERS]) {
	int i, j;
	char dummy[MAX_NUM_LETTERS];
	data_set *d = new (data_set);

	FILE *fp;
	fp = fopen(datafile, "r");
	if (!fp) {
		puts("ERROR. Missing input data file (e.g. .plu (compressed Boolean) .dat");
		exit(1);
	}
	if (fscanf(fp, "%s %d", dummy, &(d->no_inputs)) == EOF) {
		puts("Error, can't read number of inputs");
		exit(1);
	}
	if (fscanf(fp, "%s %d", dummy, &(d->no_outputs)) == EOF) {
		puts("Error, can't read number of outputs");
		exit(1);
	}
	if (fscanf(fp, "%s %d", dummy, &(d->no_tests)) == EOF) {
		puts("Error, can't read number of tests");
		exit(1);
	}
	if (d->no_tests >= MAX_NUM_DATA) {
		printf("\nERROR. Too many test cases (in datafile)\n");
		exit(0);
	}

	d->din = (data_type*) malloc(sizeof(data_type) * (d->no_tests * d->no_inputs));
	d->dout = (data_type*) malloc(sizeof(data_type) * (d->no_tests * d->no_outputs));

	for (i = 0; i < d->no_tests; i++) {
		for (j = 0; j < d->no_inputs; j++)
			d->din[i * d->no_inputs + j] = myfscanf(fp);

		for (j = 0; j < d->no_outputs; j++)
			d->dout[i * d->no_outputs + j] = myfscanf(fp);
	}

	fclose(fp);

	printf("Data file %s loaded: inputs=%d outputs=%d tests=%d \n", datafile, d->no_inputs, d->no_outputs, d->no_tests);
	return d;
}

/* calculates a perfect score (if any) */
void define_perfect(void) {

#ifdef DATA_IS_UNSIGNED_INT	/* Boolean case: using compressed truth tables */
	perfect = pow2(num_inputs) * num_outputs;
	if (num_inputs == 2)
		bit_width = 4;

	else if (num_inputs == 3)
		bit_width = 8;

	else if (num_inputs == 4)
		bit_width = 16;

	else
		bit_width = 32;

#endif				/*  */
#ifdef DATA_IS_INT
	perfect = num_tests;

#endif				/*  */
#ifdef DATA_IS_DOUBLE
	perfect = -0.01; 	// Miller style
	perfect = -0.000001;
#endif				/*  */
}

/* prints a chromosome to a file
 when append is 1, the function appends the information to the file
 when append is 0, the function creates a new file
 */
void fprint_a_chromosome(int *chromosome, const char *name, int append) {
	int i, node_label;
	int write_bracket = 1;
	FILE *fp;

	if (append)
		fp = fopen(name, "a");
	else
		fp = fopen(name, "w");
	node_label = num_inputs - 1;
	for (i = 0; i < num_nodes * num_genes_per_node; i++) {
		if ((i + 1) % num_genes_per_node == 0) {
			node_label++;
			fprintf(fp, "[%d]:%d)\t", chromosome[i], node_label);
			write_bracket = 1;
		} else {
			if (write_bracket)
				fprintf(fp, "(");

			fprintf(fp, "%d,", chromosome[i]);
			write_bracket = 0;
		}
	}
	fprintf(fp, "\t\t");
	for (i = 0; i < num_outputs; i++)
		fprintf(fp, " %d", chromosome[num_nodes * num_genes_per_node + i]);
	fprintf(fp, "\n\n");
	fclose(fp);
}

/* prints a chromosome to the screen */
void print_a_chromosome(int *chromosome) {
	int i;
	for (i = 0; i < num_nodes * num_genes_per_node; i++) {
		if ((i + 1) % num_genes_per_node == 0)
			printf(" %d\t", chromosome[i]);

		else
			printf(" %d", chromosome[i]);
	}
	printf("\t\t");
	for (i = 0; i < num_outputs; i++)
		printf(" %d", chromosome[num_nodes * num_genes_per_node + i]);
	printf("\n");
}

/* prints a chromosome to file in raw format,
 so that it can be read by the program
 */
void fprint_a_raw_chromosome(int *chromosome, char name[], int append) {
	int i;
	FILE *fp;
	if (append)
		fp = fopen(name, "a");

	else
		fp = fopen(name, "w");
	for (i = 0; i < num_nodes * num_genes_per_node; i++) {
		if ((i + 1) % num_genes_per_node == 0)
			fprintf(fp, " %d\t", chromosome[i]);

		else
			fprintf(fp, " %d", chromosome[i]);
	}
	fprintf(fp, "\t\t");
	for (i = 0; i < num_outputs; i++)
		fprintf(fp, " %d", chromosome[num_nodes * num_genes_per_node + i]);
	printf("\n");
	fclose(fp);
}

/* prints out array to a file */
void fprint_node_used(int size, int array[MAX_NUM_NODES_PLUS_INPUTS], char name[], int append) {
	int i;
	FILE *fp;
	if (append)
		fp = fopen(name, "a");

	else
		fp = fopen(name, "w");
	fprintf(fp, "\nnode_used is now\n");
	fprintf(fp, "  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20\n");
	for (i = 0; i < size; i++)
		fprintf(fp, "%3d", array[i]);
	fclose(fp);
}

/* determines from the function gene what its arity is.
 It uses the global arrays allowed_functions */
int get_arity(int function_gene) {
	int i, arity = 0;
	for (i = 0; i < num_functions; i++)
		if (function_gene == allowed_functions[i][0])
			arity = allowed_functions[i][1];
	return arity;
}

/* calculates the addresses of nodes used
 and stores them in node_to_process
 returns the number of nodes used
 */
int get_nodes_to_process(int *chromosome, int nodes_to_process[MAX_NUM_NODES]) {
	int i, j, index;
	int num_nodes_to_process;
	int node_genes[MAX_NUM_GENES_PER_NODE];
	int node_used[MAX_NUM_NODES_PLUS_INPUTS];
	int max_size_node_used;
	max_size_node_used = num_nodes + num_inputs;

	/* say all nodes not used */
	for (i = 0; i < max_size_node_used; i++)
		node_used[i] = 0;

	/* all the nodes whose output is given by the output genes are active */
	/* last num_outputs genes are all output genes */
	for (i = num_genes - num_outputs; i < num_genes; i++)
		node_used[chromosome[i]] = 1;
	for (i = max_size_node_used - 1; i >= num_inputs; i--) {
		if (node_used[i]) {

			/* get input addresses and type of this gate */
			index = num_genes_per_node * (i - num_inputs);

			/* write genes for node into array node_genes */
			for (j = 0; j < num_genes_per_node; j++)
				node_genes[j] = chromosome[index + j];

			/* each function has an arity stored in
			 allowed_functions[][1].
			 Find the nodes whose data is used
			 */
			for (j = 0; j < get_arity(node_genes[num_genes_per_node - 1]); j++)
				node_used[node_genes[j]] = 1;
		}
	}

	/* find number of used nodes */
	num_nodes_to_process = 0;
	for (i = num_inputs; i < max_size_node_used; i++) {
		if (node_used[i]) {
			nodes_to_process[num_nodes_to_process] = i;
			num_nodes_to_process++;
		}
	}
	return num_nodes_to_process;
}

/* calculates the addresses of nodes used
 and stores them in node_to_process
 returns the number of nodes used
 */
int get_nodes_to_process(int *chromosome, int nodes_to_process[MAX_NUM_NODES], int node_used[MAX_NUM_NODES_PLUS_INPUTS]) {
	int i, j, index;
	int num_nodes_to_process;
	int node_genes[MAX_NUM_GENES_PER_NODE];
//int node_used[MAX_NUM_NODES_PLUS_INPUTS];
	int max_size_node_used;
	max_size_node_used = num_nodes + num_inputs;

	/* say all nodes not used */
	for (i = 0; i < max_size_node_used; i++)
		node_used[i] = 0;

	/* all the nodes whose output is given by the output genes are active */
	/* last num_outputs genes are all output genes */
	for (i = num_genes - num_outputs; i < num_genes; i++)
		node_used[chromosome[i]] = 1;
	for (i = max_size_node_used - 1; i >= num_inputs; i--) {
		if (node_used[i]) {

			/* get input addresses and type of this gate */
			index = num_genes_per_node * (i - num_inputs);

			/* write genes for node into array node_genes */
			for (j = 0; j < num_genes_per_node; j++)
				node_genes[j] = chromosome[index + j];

			/* each function has an arity stored in
			 allowed_functions[][1].
			 Find the nodes whose data is used
			 */
			for (j = 0; j < get_arity(node_genes[num_genes_per_node - 1]); j++)
				node_used[node_genes[j]] = 1;
		}
	}

	/* find number of used nodes */
	num_nodes_to_process = 0;
	for (i = num_inputs; i < max_size_node_used; i++) {
		if (node_used[i]) {
			nodes_to_process[num_nodes_to_process] = i;
			num_nodes_to_process++;
		}
	}
	return num_nodes_to_process;
}

/**** print chromosome to file and indicate inactive genes with -1 */
void fprint_active_genes(int *chromosome, const char *name) {
	int i, j, index;
	int write_bracket;
	int node_label;
	int num_unused_nodes, num_nodes_active;
	int node_genes[MAX_NUM_GENES_PER_NODE];
	int node_used[MAX_NUM_NODES_PLUS_INPUTS];
	int max_size_node_used;
	int *active_chromosome = NULL;
	FILE *fp;
	active_chromosome = create_1d_space(num_genes);
	max_size_node_used = num_nodes + num_inputs;
	fp = fopen(name, "a");
	for (i = 0; i < num_genes; i++)
		active_chromosome[i] = -1;
	for (i = num_genes - num_outputs; i < num_genes; i++)
		active_chromosome[i] = chromosome[i];

	/* say all nodes not used */
	for (i = 0; i < max_size_node_used; i++)
		node_used[i] = 0;

	/* all the nodes whose output is given by the output genes are active */
	/* last num_outputs genes are all output genes */
	for (i = num_genes - num_outputs; i < num_genes; i++)
		node_used[chromosome[i]] = 1;
	for (i = max_size_node_used - 1; i >= num_inputs; i--) {
		if (node_used[i]) {

			/* get input addresses and type of this gate */
			index = num_genes_per_node * (i - num_inputs);

			/* write genes for node into array node_genes */
			for (j = 0; j < num_genes_per_node; j++) {
				node_genes[j] = chromosome[index + j];
				active_chromosome[index + j] = node_genes[j];
			}

			/* each function has an arity stored in
			 allowed_functions[][1].
			 Find the nodes whose data is used
			 */
			for (j = 0; j < get_arity(node_genes[num_genes_per_node - 1]); j++)
				node_used[node_genes[j]] = 1;
		}
	}
	node_label = num_inputs - 1;
	write_bracket = 1;
	for (i = 0; i < num_nodes * num_genes_per_node; i++) {
		if ((i + 1) % num_genes_per_node == 0) {
			node_label++;
			if (active_chromosome[i] < 0)
				fprintf(fp, "[*]:%d)\t", node_label);

			else
				fprintf(fp, "[%d]:%d)\t", active_chromosome[i], node_label);
			write_bracket = 1;
		} else {
			if (write_bracket == 1)
				fprintf(fp, "(");
			if (active_chromosome[i] < 0)
				fprintf(fp, "*,");

			else
				fprintf(fp, "%d,", active_chromosome[i]);
			write_bracket = 0;
		}
	}
	fprintf(fp, "\t\t");
	for (i = 0; i < num_outputs; i++)
		fprintf(fp, " %d", active_chromosome[num_nodes * num_genes_per_node + i]);
	num_unused_nodes = 0;
	for (i = num_inputs; i < num_inputs + num_nodes; i++)
		if (!node_used[i])
			num_unused_nodes++;
	num_nodes_active = num_nodes - num_unused_nodes;
	fprintf(fp, "\nnumber of active gates is %d\n\n", num_nodes_active);
	fclose(fp);
	free(active_chromosome);
}

/* generate a starting population
 from a chromosome read from a file (cgp.chr)
 */
void read_from_chrom(int **chromosomes) {
	int i, j;
	FILE *fp;
	fp = fopen("cgp.chr", "r");
	if (!fp) {
		puts("Missing file cgp.chr (contains a chromosome)");
		exit(1);
	} else {

		/* make starting population copies of loaded chromosome */
		for (j = 0; j < population_size; j++) {
			if (j == 0) {
				i = 0;

				do {
					if (fscanf(fp, "%d", &chromosomes[j][i])) {
						puts("Can't read the chromosome");
						exit(1);
					}
					i++;
				} while (!feof(fp));
				if (i != num_genes) {
					puts("ERROR. Number of genes in cgp.chr does not match the expected number");
					printf("\nnum_genes required is %d, num_genes read is %d", num_genes, i);
					puts("Check the number of genes in the .par file");
					exit(0);
				}
			} else {
				for (i = 0; i < num_genes; i++)
					chromosomes[j][i] = chromosomes[0][i];
			}
		}
		fclose(fp);
	}
}

/* this decodes the cgp chromosome.
 It is given data_inputs corresponding to a single
 test case and it calculates what the cgp genotype gives
 for the data outputs (cgp_outputs)
 It only processes nodes that are used (nodes_to_process)
 */
void decode_cgp(int *chromosome, data_type data_inputs[MAX_NUM_DATA][MAX_NUM_INPUTS], data_type cgp_outputs[MAX_NUM_OUTPUTS], int num_nodes_to_process, int nodes_to_process[MAX_NUM_NODES],
		int fitness_test) {
	int i, j;
	int first_node_gene, function_type;
	int node_index;
	data_type in[MAX_NUM_GENES_PER_NODE];
	data_type output[MAX_OUTPUT_SIZE];
	data_type constant_inputs[MAX_NUM_CONSTANT_INPUTS] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

	/* load user defined constants into output array */
	for (j = 0; j < num_constant_inputs; j++)
		output[j] = constant_inputs[j];

	/* load test data into output array */
	for (j = num_constant_inputs; j < num_inputs; j++)
		output[j] = data_inputs[fitness_test][j - num_constant_inputs];

	/* only process nodes that are used */
	for (j = 0; j < num_nodes_to_process; j++) {

		/* get address of node */
		node_index = nodes_to_process[j] - num_inputs;

		/* get address of first used gene in node */
		first_node_gene = num_genes_per_node * node_index;
		for (i = 0; i < num_genes_per_node - 1; i++) { /* get input data to node */
			in[i] = output[chromosome[first_node_gene + i]];
		}
		function_type = chromosome[first_node_gene + num_genes_per_node - 1]; /* get node function */
		output[node_index + num_inputs] = node_type(in, function_type); /* compute output of node and store */
	}

	/* process outputs */
	for (j = 0; j < num_outputs; j++)
		cgp_outputs[j] = output[chromosome[num_genes - num_outputs + j]]; /* get output from nodes referenced in output genes */
}

/**
 *  should be used with data from a data_set structure
 *
 * this decodes the cgp chromosome.
 * It is given data_inputs corresponding to a single
 * test case and it calculates what the cgp genotype gives
 * for the data outputs (cgp_outputs)
 * It only processes nodes that are used (nodes_to_process)
 */
void decode_cgp2(int *chromosome, data_set *d, data_type cgp_outputs[MAX_NUM_OUTPUTS], int num_nodes_to_process, int nodes_to_process[MAX_NUM_NODES], int fitness_test) {
	int i, j;
	int first_node_gene, function_type;
	int node_index;
	data_type in[MAX_NUM_GENES_PER_NODE];
	data_type output[MAX_OUTPUT_SIZE];
	data_type constant_inputs[MAX_NUM_CONSTANT_INPUTS] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

	/* load user defined constants into output array */
	for (j = 0; j < num_constant_inputs; j++)
		output[j] = constant_inputs[j];

	/* load test data into output array */
	for (j = num_constant_inputs; j < num_inputs; j++)
		// output[j] = data_inputs[fitness_test][j - num_constant_inputs];
		output[j] = d->din[d->no_inputs * fitness_test + j - num_constant_inputs];

	/* only process nodes that are used */
	for (j = 0; j < num_nodes_to_process; j++) {

		/* get address of node */
		node_index = nodes_to_process[j] - num_inputs;

		/* get address of first used gene in node */
		first_node_gene = num_genes_per_node * node_index;
		for (i = 0; i < num_genes_per_node - 1; i++) { /* get input data to node */
			in[i] = output[chromosome[first_node_gene + i]];
		}
		function_type = chromosome[first_node_gene + num_genes_per_node - 1]; /* get node function */
		output[node_index + num_inputs] = node_type(in, function_type); /* compute output of node and store */
	}

	/* process outputs */
	for (j = 0; j < num_outputs; j++)
		cgp_outputs[j] = output[chromosome[num_genes - num_outputs + j]]; /* get output from nodes referenced in output genes */
}

#ifdef DATA_IS_UNSIGNED_INT
/* counts how many bits the cgp calculated 32-bit output
 has in common with compressed truth table output read from read_data()
 */
double correctness_test(data_type data_output, data_type cgp_output) {
	register int i;
	int result = 0;
	data_type temp;
	temp = ~data_output ^ cgp_output; /* xnor the calculated output with the desired */
	for (i = 0; i < bit_width; i++) /* examine the 32 bits */
		result = result + getbit(temp, i);
	return (double) result;
}
#endif				/*  */

#ifdef DATA_IS_INT

/* checks whether the output integer cgp produces
 is the same as the desired output integer
 Using HITS based fitness is a very crude measurement of correctness!
 */
double correctness_test(data_type data_output, data_type cgp_output)
{
	double result = 0.0;
	double x;
	x = fabs(cgp_output - data_output);

	/* hits based fitness */
#ifdef HITS_BASED_FITNESS
	if (x == 0)
	result = 1;

#else				/*  */
	result = 1.0 / (1.0 + x); /* error based fitness */

#endif				/*  */
	return result;
}

#endif				/*  */

#ifdef DATA_IS_DOUBLE

/* checks how close the output double cgp produces
 is to the desired output integer
 returns a 1 each time the absolute difference between
 the desired output and the cgp output is less than a
 user defined error. This is a hits based fitness measure
 suitable for symbolic regression problems
 */
double correctness_test(data_type data_output, data_type cgp_output) {
	double result = 0.0;
	double x;
	x = fabs(cgp_output - data_output);

	return x;

	/* hits based fitness */
//#ifdef HITS_BASED_FITNESS
//	if (x < ERROR_THRESHOLD)
//	result = 1;
//#else				/*  */
//	result = 1.0 / (1.0 + x); /* error based fitness */
//#endif				/*  */
//	return result;
}

#endif				/*  */

/* evaluate the fitness at a single test point */
double evaluate_cgp_outputs(data_type cgp_outputs[MAX_NUM_OUTPUTS], int test) {
	int i;
	double fit = 0.0, x;
	for (i = 0; i < num_outputs; i++) {
		x = correctness_test(data_outputs[test][i], cgp_outputs[i]);
		fit = fit + x;
	}

	return fit;
}

/* this is the EA fitness function
 */
double fitness(int *chromosome, long *processed_nodes) {
	int fitness_test;
	int num_nodes_to_process;
	int nodes_to_process[MAX_NUM_NODES];
	double fit = 0.0;
	data_type cgp_outputs[MAX_NUM_OUTPUTS];

	/* find out how many nodes there are in the phenotype */
	num_nodes_to_process = get_nodes_to_process(chromosome, nodes_to_process);

	/* apply all fitness tests */
	for (fitness_test = 0; fitness_test < num_tests; fitness_test++) {
		decode_cgp(chromosome, data_inputs, cgp_outputs, num_nodes_to_process, nodes_to_process, fitness_test);
		fit = fit + evaluate_cgp_outputs(cgp_outputs, fitness_test);
	}

// MAE

	fit = fit / (((float) num_tests) * ((float) num_outputs));

	fit = -fit;

	/* if a perfect solution is found and shrink_phenotype is set to 1
	 this adds how many unused nodes there are to the fitness score.
	 this results in phenotypes shrinking
	 */
	if ((fit >= perfect) && (shrink_phenotype))
		fit = fit + num_nodes - num_nodes_to_process;

// increase the number of processed nodes
	*processed_nodes += num_nodes_to_process * num_tests;
	return fit;
}

/* evaluate the fitness at a single test point */
double evaluate_cgp_outputs(data_type cgp_outputs[MAX_NUM_OUTPUTS], data_set *d, int test) {
	int i;
	double fit = 0.0, x;
	for (i = 0; i < d->no_outputs; i++) {
//		x = correctness_test(data_outputs[test][i], cgp_outputs[i]);
		x = fabs(cgp_outputs[i] - d->dout[test * d->no_outputs + i]);
		fit = fit + x;
	}

	return fit;
}

/* this is the EA fitness function
 */
double fitness2(int *chromosome, long *processed_nodes, data_set *d) {
	int fitness_test;
	int num_nodes_to_process;
	int nodes_to_process[MAX_NUM_NODES];
	double fit = 0.0;
	data_type cgp_outputs[MAX_NUM_OUTPUTS];

//	for (int j = 0; j < num_tests; j++)
//		for (int i = 0; i < num_inputs - num_constant_inputs; i++) {
//			printf("test: %d input: %d training: %f test: %f\n", j, i, data_inputs[j][i], d->din[d->no_inputs * j + i]);
//		}

	/* find out how many nodes there are in the phenotype */
	num_nodes_to_process = get_nodes_to_process(chromosome, nodes_to_process);

	/* apply all fitness tests */
	for (fitness_test = 0; fitness_test < d->no_tests; fitness_test++) {
		decode_cgp2(chromosome, d, cgp_outputs, num_nodes_to_process, nodes_to_process, fitness_test);
		fit = fit + evaluate_cgp_outputs(cgp_outputs, d, fitness_test);
	}

	// MAE

	fit = fit / (((float) d->no_tests) * ((float) d->no_outputs));

	fit = -fit;

	/* if a perfect solution is found and shrink_phenotype is set to 1
	 this adds how many unused nodes there are to the fitness score.
	 this results in phenotypes shrinking
	 */
	if ((fit >= perfect) && (shrink_phenotype))
		fit = fit + num_nodes - num_nodes_to_process;

	// increase the number of processed nodes
	*processed_nodes += num_nodes_to_process * d->no_tests;
	return fit;
}

/* This calculates the limits that are used in the calculation
 of allowed gene values (alleles) */
void get_gene_limits(int column, int *limit_min, int *limit) {
	int limit_max;
	limit_max = num_inputs + column * num_rows;
	if (column < levels_back)
		*limit_min = 0;

	else
		*limit_min = num_inputs + (column - levels_back) * num_rows;
	*limit = limit_max - (*limit_min);
}

/* returns a random valid connection gene that
 obeys the constraints imposed by levels_back.
 Also allows program inputs to disobey levels_back */
int get_connection_gene(int limit_min, int limit) {
	int limit_plus, rand_num;
	int gene;
	if (limit_min == 0)
		gene = newrand(limit);

	else { /* allows inputs to disobey levels_back */

		limit_plus = limit + num_inputs;
		rand_num = newrand(limit_plus);
		if (rand_num < limit)
			gene = rand_num + limit_min;

		else
			gene = rand_num - limit;
	}
	return gene;
}

/* returns a random valid function gene */
int get_function_gene(void) {
	return allowed_functions[newrand(num_functions)][0];
}

/* returns a random valid output gene */
int get_output_gene(void) {
	int limit_min, limit;
	int output_gene;
	limit_min = num_inputs + (num_cols - levels_back) * num_rows;
	limit = levels_back * num_rows;
	output_gene = newrand(limit) + limit_min;
	return output_gene;
}

/* checks to see if the gene is not an output gene */
int is_not_output_gene(int gene) {
	return (gene < num_genes_per_node * num_nodes);
}

/* checks to see if the gene is a function gene */
int is_function_gene(int gene, int locus) {
	return (is_not_output_gene(gene) && (locus == (num_genes_per_node - 1)));
}

/* generates a random chromosome. Used by initialise */
void generate_a_random_chromosome(int *chromosome) {
	int i, count = 0;
	int row, col, limit, limit_min;
	for (col = 0; col < num_cols; col++) {
		get_gene_limits(col, &limit_min, &limit);
		for (row = 0; row < num_rows; row++) {

			/* get random connection genes */
			for (i = 0; i < num_genes_per_node - 1; i++)
				chromosome[count + i] = get_connection_gene(limit_min, limit);

			/* get random function gene */
			chromosome[count + num_genes_per_node - 1] = get_function_gene();
			count = count + num_genes_per_node;
		}
	}

	/* get random function genes */
	for (i = 0; i < num_outputs; i++)
		chromosome[count + i] = get_output_gene();
}

/* creates initial population of chromosomes
 either having been generated from a single
 chromosome from a file or by generating
 an entire random population
 */
void initialise(int **chromosomes) {
	int i;
	if (run_from_chrom)
		read_from_chrom(chromosomes);

	else
		/* generate random population */
		for (i = 0; i < population_size; i++)
			generate_a_random_chromosome(chromosomes[i]);
}

/* calculate best population fitness and the best chromosome */
double get_best_chromosome(int **chromosomes, int *best_chromosome, double previous_best_fitness, int gen, int *number_evaluations) {

	int i;
	double fitness_max, fit;
	int best_member;
	long processed_nodes = 0;

	fitness_max = -1.0;
	best_member = 0;
	for (i = 0; i < population_size; i++) {
		if ((i == population_size - 1) && (gen > 1))
			fit = previous_best_fitness;

		else {
			fit = fitness(chromosomes[i], &processed_nodes);
			(*number_evaluations)++;
		}

		if (fit > fitness_max) {
			fitness_max = fit;
			best_member = i;
		}

		/* break out of this as soon as we get a perfect score
		 and shrink_phenotype is not required */
		if ((fit >= perfect) && (shrink_phenotype == 0))
			break;
	}

// if explicit neutrality off, child has to be strictyl better than parent to proceed as parent to next generation
	if ((explicit_neutrality == 0) && (fitness_max == previous_best_fitness))
		return fitness_max;

	/* store the best chromosome */
	for (i = 0; i < num_genes; i++)
		best_chromosome[i] = chromosomes[best_member][i];
	return fitness_max;
}

/* allocate space for 2 dimensional array
 e.g. a population of chromosomes */
int** create_2d_space(int num_horizontals, int num_verticals) {
	int i;
	int **array2d = NULL;

	/* create space for pointers to int pointers */
	array2d = (int**) calloc(num_verticals, sizeof(int*));
	if (array2d == NULL) {
		printf("ERROR.Can not allocate space for %d many int pointers\n", num_verticals);
		exit(0);
	}

	/* create array of pointers to ints  */
	for (i = 0; i < num_verticals; i++) {
		array2d[i] = create_1d_space(num_horizontals);
		if (array2d[i] == NULL) {
			printf("ERROR.Not enough memory for int pointer arrays of length %d\n", num_horizontals);
			exit(0);
		}
	}
	return array2d;
}

/* allocate space for a 1d array of with size items */
int* create_1d_space(int size) {
	int *array1d = NULL;
	array1d = (int*) calloc(size, sizeof(int));
	if (array1d == NULL) {
		printf("ERROR.Not enough memory for a 1d array of length %d\n", size);
		exit(0);
	}
	return array1d;
}

/* release memory */
/* if this is for chromosomes then num_verticals is population_size */
void free_array2d(int num_verticals, int **array2d) {
	int i;

	/* free 1darray of pointers  */
	for (i = 0; i < num_verticals; i++)
		free(array2d[i]);
	free(array2d);
}

void write_generation_to_screen(int gen_index) {
	if (report_interval > 0 && gen_index % report_interval == 0)
		printf("\nGENERATION is %d", gen_index);
}

/* writes generation and fitness of best in population */
void write_progress_info_to_screen(int generation, long nodes_processed, double fit) {
	printf("\nEval. %d Nodes processed %ld Best fitness is now %8.5lf", generation, nodes_processed, fit);

	string histFile;

	ostringstream s1;

	s1 << "histogram_" << problem << "_" << mu << "_" << population_size << ".txt";

	histFile = s1.str();

	ofstream myfile(histFile, std::ios_base::app);

	if (myfile.is_open()) {
		myfile << (fit * -1.0) << "\n";
		myfile.close();
	}

}

/* writes out chromosome to file defined by string prog */
void write_progress_info_to_file(char *prog, int gen, double best_fit, int *best_chromosome) {
	FILE *fp;

// discard for cluster experiments for cec gecco 2017
	return;
	fp = fopen(prog, "a");
	fprintf(fp, "\nGENERATION is %u     Best fitness is now %8.5lf", gen, best_fit);
	fprintf(fp, "\nThe chromosome is\n");
	fclose(fp);
	fprint_a_chromosome(best_chromosome, prog, 1);
	fprint_active_genes(best_chromosome, prog);
}

/* checks to see if the best fitness in the population has improved.
 writes the generation, the new best fitness and the improved chromosome
 to the progress report (if progress report is set to 1)
 */
void check_if_improvement(double best_fit, double *previous_best_fit, int *best_gen, int gen, char prog[MAX_NUM_LETTERS], int *best_chromosome) {
	if (best_fit > *previous_best_fit) { /* we have an improvement */
		if (progress_report) {
			write_progress_info_to_screen(gen, -1, best_fit);
			write_progress_info_to_file(prog, gen, best_fit, best_chromosome);
		}
		*best_gen = gen; /* update the best generation */
		*previous_best_fit = best_fit; /* update previous best fitness */
	}
}

/* report on results of evolutionary run in cgp.txt */
void write_result_of_EA_to_file(int run, int bgen, double best_fit, int *best_chromosome) {
	FILE *fp;
	fp = fopen("cgp.txt", "a");
	fprintf(fp, "Run %d and gen %d achieved fitness %6.2lf\n", run, bgen, best_fit);
	fprintf(fp, "Here is the chromosome\n");
	fclose(fp);
	fprint_a_chromosome(best_chromosome, "cgp.txt", 1);
	fprint_active_genes(best_chromosome, "cgp.txt");
}

float rand01() {
	float op_tmp = RAND_MAX;
	float op = rand();
	return (op / op_tmp);
}

void setup_report_files(int run, char prog[MAX_NUM_LETTERS]) {
	char runstring[MAX_NUM_LETTERS];
//	FILE *fp;
	sprintf(runstring, "%d", run); /* store run as characters */
	if (progress_report > 0) {
		strcpy(prog, "cgp");
		strcat(prog, runstring);
		strcat(prog, ".prg"); /* create .prg file name */

		// discard all debug writes to files, implemented for the cluster-based evaluation for CEC/GECCO 2017
		//fp = fopen(prog, "w"); /* create empty .prg file */
		//fclose(fp);
	}
}

double polyEvaluate(int order, double x) {
	double res = 0.0;

	for (int i = 0; i <= order; i++) {
		res += pow(x, i);
	}

	return res;
}

double pol3(double x) {
	return x * x * x + x * x + x;
}

double koza1(double x) {
	return x * x * x * x + x * x * x + x * x + x;
}

double koza3(double x) {
	return x * x * x * x * x * x - 2 * x * x * x * x + x * x;
}

double koza2(double x) {
	return x * x * x * x * x - 2 * x * x * x + x;
}

double nguyen3(double x) {
	return x * x * x * x * x + x * x * x * x + x * x * x + x * x + x;
}

double nguyen4(double x) {
	return x * x * x * x * x * x + x * x * x * x * x + x * x * x * x + x * x * x + x * x + x;
}

double nguyen5(double x) {
	return sin(x * x) * cos(x) - 1;
}

double nguyen6(double x) {
	return sin(x) * sin(x + x * x);
}

double nguyen7(double x) {
	return log(x + 1) + log(x * x + 1);
}

double nguyen8(double x) {
	return sqrt(x);
}

double nguyen9(double x, double y) {
	return sin(x) + sin(y * y);
}

double nguyen10(double x, double y) {
	return 2.0 * sin(x) * cos(y);
}

double keijzer4(double x) {
	return x * x * x * exp(-1.0 * x) * cos(x) * sin(x) * (sin(x) * sin(x) * cos(x) - 1);
}

double keijzer6(double x) {
	double sum = 0;
	double fx = floor(x);
	int i;
	for (i = 1; i < fx + 1; i++)  // up to and including floor(x)
			{
		sum += (1.0 / i);
	}
	return sum;
}

double korns12(double x, double w) {
	return 2.0 - (2.1 * (cos(9.8 * x) * sin(1.3 * w)));
}

double vladisleva4(double xs[]) {
	double sum = 0;
	int i;
	for (i = 0; i < 5; i++) {
		sum += (xs[i] - 3) * (xs[i] - 3);
	}

	return 10.0 / (5.0 + sum);
}

double vladislavleva6(double x, double y) {
	return 6.0 * sin(x) * cos(y);
}

double pagie1(double x, double y) {
	double res = (1.0 / (1.0 + pow(x, -4.0))) + (1.0 / (1.0 + pow(y, -4.0)));
	return res;
}

double test_function(double x) {
	return x * x + x;
}

int get_num_erc(int problem) {
	if (problem == 11 || problem == 12) {
		return 1;
	} else {
		return 2;
	}
}

void calculate_data2(int problem) {

	int i, j;
	for (i = 0; i < num_tests; i++) {
		for (j = 0; j < num_outputs; j++) {
			if (problem == 0) {
				data_outputs[i][j] = koza1(data_inputs[i][0]);
			} else if (problem == 1) {
				data_outputs[i][j] = koza2(data_inputs[i][0]);
			} else if (problem == 2) {
				data_outputs[i][j] = koza3(data_inputs[i][0]);
			} else if (problem == 3) {
				data_outputs[i][j] = nguyen4(data_inputs[i][0]);
			} else if (problem == 4) {
				data_outputs[i][j] = nguyen5(data_inputs[i][0]);
			} else if (problem == 5) {
				data_outputs[i][j] = nguyen6(data_inputs[i][0]);
			} else if (problem == 6) {
				data_outputs[i][j] = nguyen7(data_inputs[i][0]);
			} else if (problem == 7) {
				data_outputs[i][j] = nguyen8(data_inputs[i][0]);
			} else if (problem == 8) {
				data_outputs[i][j] = nguyen9(data_inputs[i][0], data_inputs[i][1]);
			} else if (problem == 9) {
				data_outputs[i][j] = nguyen10(data_inputs[i][0], data_inputs[i][1]);
			} else if (problem == 10) {
				data_outputs[i][j] = pagie1(data_inputs[i][0], data_inputs[i][1]);
			} else if (problem == 11) {
				data_outputs[i][j] = keijzer6(data_inputs[i][0]);
			} else if (problem == 12) {
				data_outputs[i][j] = keijzer4(data_inputs[i][0]);
			} else if (problem == 13) {
				data_outputs[i][j] = pol3(data_inputs[i][0]);
			} else if (problem == 14) {
				data_outputs[i][j] = polyEvaluate(order, data_inputs[i][0]);
			} else if (problem == 15) {
				data_outputs[i][j] = nguyen3(data_inputs[i][0]);
			} else {
				exit(1);
			}
		}
	}
}

