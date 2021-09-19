/* cgp.h
 Julian F. Miller (c), 2009
 version 1.1 of first public release on 20-July-2009
 Dept. of Electronics, University of York, UK
 */

#ifndef __CGP_H__
#define __CGP_H__

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "cgp.h"
#include "mutation.h"
#include "optimization_algorithms.h"
#include "statistical_methods.h"


/* invoke alternative data types by uncommenting the one you want
 and commenting out the others
 */

#define DATA_IS_UNSIGNED_INT
//#define DATA_IS_INT
// #define DATA_IS_DOUBLE

/* this if statement defines the type of data that
 is processed by cgp genotypes
 if you have a user defined data type
 that is not on this list then you need to add a case
 to the if statement
 Note this will mean that you will have to define a new function set
 this is defined in the c program file node_function.c
 */

/* this defines how the maximum constant inputs are input to cgp */
#define MAX_NUM_CONSTANT_INPUTS 3

#ifdef DATA_IS_UNSIGNED_INT
typedef unsigned long data_type;
#endif
#ifdef DATA_IS_INT
typedef int data_type;
#endif
#ifdef DATA_IS_DOUBLE
typedef double data_type;
#endif

/* this global defines exactly how many constant inputs are used
 in cgp. This varies for the data_type chosen. For unsigned int
 data type num_constant_inputs is set to zero */
extern int num_constant_inputs;

#define PI							3.1415926
#define MAX_NUM_ROWS				200
#define MAX_NUM_COLS				200
#define MAX_NUM_INPUTS				16
#define MAX_NUM_OUTPUTS				256
#define MAX_NUM_DATA				2056
#define MAX_NUM_NODES				MAX_NUM_ROWS*MAX_NUM_COLS
#define MAX_NUM_NODES_PLUS_INPUTS	MAX_NUM_NODES+MAX_NUM_INPUTS+MAX_NUM_CONSTANT_INPUTS
#define MAX_NUM_GENES_PER_NODE		4
#define MAX_NUM_GENES				MAX_NUM_GENES_PER_NODE*MAX_NUM_NODES
#define MAX_OUTPUT_SIZE				MAX_NUM_INPUTS+MAX_NUM_CONSTANT_INPUTS+MAX_NUM_NODES+MAX_NUM_OUTPUTS
#define MAX_NUM_RUNS				500
#define ERROR_THRESHOLD				0.0001

/*
 #define HITS_BASED_FITNESS
 */

#define MAX_NUM_CHROMOSOMES			99999999  /* max population size */
#define MAX_NUM_LETTERS				100
#define MAX_NUM_FUNCTIONS			20
#define MAXNUM						4294967295

/* these are all read from the .par file and never change after that */
extern int population_size;
extern int mu;
extern double per_cent_mutate;
extern int genes_mutate;
extern long long max_number_evaluations;
extern int num_generations;
extern int num_fitness_evaluations;
extern int num_runs_total;
extern int num_rows;
extern int num_cols;
extern int levels_back;
extern int progress_report;
extern int report_interval;
extern unsigned global_seed;
extern int save_best_chrom;
extern int run_from_chrom;
extern int shrink_phenotype;

// global mutation variables
extern int mutate_active_wires;
extern int mutate_inactive_wires;
extern int mutate_active_node_func;
extern int mutate_inactive_node_func;
extern int randomize_inactive_genes; // randomizes all inactive genes after a mutation
extern int log_neighbor_fitness_values;

/*  global constants calculated in get_parameters */
extern int num_functions;
extern int num_genes;
extern int num_nodes, num_genes_per_node;
extern int number[MAX_NUM_FUNCTIONS];
extern int order;
extern char node_types[MAX_NUM_FUNCTIONS][20];

//enum fitness_type_def fitness_type;
//enum algorithm_def algorithm;

/* this stores the node function address in
 allowed_functions[][0] and its arity in
 allowed_functions[][1] */
extern int allowed_functions[MAX_NUM_FUNCTIONS][2];

/* data defining the computational problem read from .dat file */
extern int num_inputs, num_outputs, num_tests, num_erc;
extern data_type data_inputs[MAX_NUM_DATA][MAX_NUM_INPUTS];
extern data_type data_outputs[MAX_NUM_DATA][MAX_NUM_INPUTS];

/* calculated global constants  */
extern double perfect;
extern int bit_width;

/* macros */
#define pow2(x) (1<<x)                                   /* returns 2 to power x (x>=0) */
#define getbit(decimal,nthbit) ((decimal>>nthbit) & 01)  /* gets nth bit */

extern int problem;
extern int eval_method;

// SA parameter struct
struct sa_parameters {
	double tstart;				// start temperature
	double tstop;				// stop temperature
	int cooling_scheme;			// cooling scheme
	int steps_in_temp_level;			// how many children should be evolved and tested before updating the temperature
	double rejection_factor;			// how many iterations in a temperature level have been not successful
};

extern struct sa_parameters sa_parameters_global;

/* function prototypes */

void validate_command_line(int argc, char *argv[], char parfile[], char datafile[]);

void get_parameters(char *parfile, char *datafile);

void write_cgp_info(FILE *fp, char command[], char *parameterfile, char *datafile);

int newrand(int range);

data_type myfscanf(FILE *fp);

void read_data(const char *datafile);

void define_perfect(void);

void fprint_a_chromosome(int *chromosome, const char name[], int append);

void print_a_chromosome(int *chromosome);

void fprint_a_raw_chromosome(int *chromosome, char name[], int append);

void fprint_node_used(int size, int array[MAX_NUM_NODES_PLUS_INPUTS], char name[], int append);

int get_arity(int function_gene);

int get_nodes_to_process(int *chromosome, int nodes_to_process[MAX_NUM_NODES]);

void fprint_active_genes(int *chromosome, const char *name);

void read_from_chrom(int **chromosomes);

data_type node_type(data_type in[MAX_NUM_GENES_PER_NODE], int function_gene);

void decode_cgp(int *chromosome, data_type data_inputs[MAX_NUM_DATA][MAX_NUM_INPUTS], data_type node_outputs[MAX_OUTPUT_SIZE], data_type cgp_outputs[MAX_NUM_OUTPUTS], int num_nodes_to_process,
		int nodes_to_process[MAX_NUM_NODES], int fitness_test);

double correctness_test_boolean(data_type data_output, data_type cgp_output);

double evaluate_cgp_outputs(data_type cgp_outputs[MAX_NUM_OUTPUTS], int test);

double fitness(int *chromosome, long *num_nodes_processed);

double best_inner_node_fitness(int *chromosome, long *processed_nodes);

double fast_best_inner_node_fitness(int *chromosome, long *processed_nodes);

void get_gene_limits(int column, int *limit_min, int *limit);

int get_connection_gene(int limit_min, int limit);

int get_function_gene(void);

int get_output_gene(void);

int is_not_output_gene(int gene);

int is_function_gene(int gene, int locus);

void generate_a_random_chromosome(int *chromosome);

void initialise(int **chromosomes);

double get_best_chromosome(int **chromosomes, int *best_chromosome, double previous_best_fitness, int gen, int *number_evaluations);

int** create_2d_space(int num_horizontals, int num_verticals);

int* create_1d_space(int size);

void free_array2d(int num_verticals, int **array2d);

void write_generation_to_screen(int gen_index);

void write_progress_info_to_screen(int generation, long nodes_processed, double fit);

void write_result_of_EA_to_file(int run, int bgen, double best_fit, int *best_chromosome);

void check_if_improvement(double best_fit, double *previous_best_fit, int *best_gen, int gen, char prog[MAX_NUM_LETTERS], int *best_chromosome);

void setup_report_files(int run, char prog[MAX_NUM_LETTERS]);

double get_standard_deviation(int num_items, double average, double items[MAX_NUM_RUNS]);

int parse_parameters(int argcc, char *argv[], char **envp);

void calculate_data(int problem);

void write_progress_info_to_file(char prog[MAX_NUM_LETTERS], int gen, double best_fit, int *best_chromosome);

void usage(const char *self);

float rand01();

#endif
