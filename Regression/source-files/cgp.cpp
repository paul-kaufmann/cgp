/* cgp.c
 Julian F. Miller (c), 2009
 version 1.1 of first public release on 20-July-2009
 Dept. of Electronics, University of York, UK

 The code was extended by Paul Kaufmann 16 September 2016
 1) almost all parameters have to be specified now at command line and not through the configuration file,
 the configuration file have to contain only the functional table definition
 2) simulated annealing is added
 3) the random seed is always read from /dev/urandom
 4) variable definition in cgp.h moved to cgp-functions.c all other file have to declare global variables as extern
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "cgp.h"
#include "mutation.h"
#include "optimization_algorithms.h"
#include "statistical_methods.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

char *parfile;
char *datafile;

/* these are all read from the .par file and never change after that */
int num_constant_inputs;
int population_size;
int mu;
double per_cent_mutate;
int mutate_genes;
long long max_number_evaluations;
int num_generations;
int num_fitness_evaluations;
int num_runs_total;
int num_rows;
int num_cols;
int levels_back;
int progress_report;
int report_interval;
unsigned global_seed;
int save_best_chrom;
int run_from_chrom;
int shrink_phenotype;
int explicit_neutrality;

int mutate_active_wires;
int mutate_inactive_wires;
int mutate_active_node_func;
int mutate_inactive_node_func;
int randomize_inactive_genes; // randomizes all inactive genes after a mutation

/*  global constants calculated in get_parameters */
int num_functions;
int num_genes;
int num_nodes, num_genes_per_node;
int number[MAX_NUM_FUNCTIONS];
int order;
char node_types[MAX_NUM_FUNCTIONS][20];

/* this stores the node function address in
 allowed_functions[][0] and its arity in
 allowed_functions[][1] */
int allowed_functions[MAX_NUM_FUNCTIONS][2];

/* data defining the computational problem read from .dat file */
int num_inputs, num_outputs, num_tests, num_erc;
data_type data_inputs[MAX_NUM_DATA][MAX_NUM_INPUTS];
data_type data_outputs[MAX_NUM_DATA][MAX_NUM_INPUTS];

/* calculated global constants  */
double perfect;
int bit_width;

int problem;
int eval_method;

struct sa_parameters sa_parameters_global;

data_set *test_data_set;

void usage(const char *self) {
	printf("usage: %s PARfile PLUfile <options>\n", self);
	printf("\n");
	printf("-i <file>           test data for symbolic regression benchmarks\n");
	printf("\n");

	printf("-a <algorithm>      %d = late acceptance HC, %d = simulated annealing\n", LAHC, SIMULATED_ANNEALING);
	printf("                    %d = random walk\n", RANDOM_WALK);
	printf("                    %d = (mu + lambda=paopulation size) ES\n", PLUS_ES);
	printf("                    %d = Hill Climber\n", HC);
	printf("                    %d = Evolutionary Algorithm\n", EA);
	printf("                    %d = Steady State EA\n", SSEA);
	printf("-t <value>          number of runs (repetitions)\n");
	printf("-r <value>          number of CGP rows\n");
	printf("-c <value>          number of CGP columns\n");
	printf("-l <value>          levels back\n");
	printf("-n <value>          stop condition: maximal number of fitness evaluations\n");
	printf("\n");

	printf("Mutation parameters\n");
	printf("-g <0/1>            1=allow / 0=disallow mutating active wires, default \"-g 1\"\n");
	printf("-h <0/1>            1=allow / 0=disallow mutating inactive wires, default \"-h 1\"\n");
	printf("-j <0/1>            1=allow / 0=disallow mutating active node functions, default \"-j 1\"\n");
	printf("-k <0/1>            1=allow / 0=disallow mutating inactive node functions, default \"-k 1\"\n");
	printf("-s <0/1>            randomize inactive genes after a mutation, default \"-s 0\"\n");

	printf("\nHill and Late Acceptance Hill Climber parameters\n");
	printf("                    no parameters\n");
	printf("                    population size is 1\n");
	printf("                    mutation rate is 1 active gene\n");
	printf("                    list length of LAHC is 20\n");

	printf("\nSimulated Annealing parameters\n");
	printf("-1 <value>          start temperature\n");
	printf("-2 <value>          stop temperature\n");
	printf("-3 <value>          cooling scheme\n");
	printf("\n");

	printf("\nmu + lambda(=population_size) ES\n");
	printf("-q <value>          mu - number of parent individual\n");
	printf("-p <value>          lambda=population size\n");
	printf("-m <value>          mutation rate in per cent\n");
	printf("-e <0/1>            interprete mutation rate as number of genes to mutate \"-e 1\", default \"-e 0 \"\n");
	printf("\n");

	printf("Return value: the program returns the costs in the last line (second token)\n");
	printf("If there is some error, MAX_INT is returned as costs, see last line of this text\n");
	printf("-b <value>          %d = return median number of fitness evaluations (default)\n", MEDIAN_EVALUATIONS);
	printf("                    %d = return median fitness\n", MEDIAN_FITNESS);
	printf("                    %d = return processed nodes\n", NODES_PROCESSED);
	printf("                    %d = return Computational Effort\n", CE_FITNESS);
	printf("                    %d = return Computational Effort[processed nodes]\n", CE_NODES_PROCESSED);
	printf("\n");

	printf("-f <value>          specify the index of the symbolic regression benchmark");

	printf("\n\nBest %14.2f\n", (float) INT_MAX);
	exit(1);
}

int parse_parameters(int argcc, char **argvv, char **envp) {

	if (argcc < 3)
		usage(*argvv);

	parfile = argvv[1];
	datafile = argvv[2];

	algorithm_type = LAHC;
	population_size = -1;
	mu = -1;
	per_cent_mutate = 0.0;
	mutate_genes = 0;
	max_number_evaluations = -1;
	num_runs_total = -1;
	num_rows = -1;
	num_cols = -1;
	levels_back = -1;
	progress_report = 1;
	report_interval = 0;
	save_best_chrom = 0;
	run_from_chrom = 0;
	shrink_phenotype = 0;
	fitness_type = MEDIAN_EVALUATIONS;
	explicit_neutrality = 1;

	mutate_active_wires = 1;
	mutate_inactive_wires = 1;
	mutate_active_node_func = 1;
	mutate_inactive_node_func = 1;
	randomize_inactive_genes = 0;

	sa_parameters_global.cooling_scheme = -1;
	sa_parameters_global.rejection_factor = -1;
	sa_parameters_global.steps_in_temp_level = -1;
	sa_parameters_global.tstart = -1;
	sa_parameters_global.tstop = -1;

	// BSD's getopt can't handle mixed parameters/options command
	// so if in BSD then skip two first parameters and go directly to
	// parsing options
#ifdef __APPLE__
	int nr_parameters = 2;
#else
	int nr_parameters = 0;
#endif

	//parse command line parameters
	char opt;
	while ((opt = getopt(argcc - nr_parameters, argvv + nr_parameters, "a:t:r:c:l:n:p:m:1:2:3:4:5:b:q:dg:h:j:k:s:f:e:i:")) != -1) {
		switch (opt) {
		case 'i':
			test_data_set = read_data2(optarg);
			break;
		case 'a':
			algorithm_type = static_cast<algorithm_type_def>(atoi(optarg));
			break;
		case 't':
			num_runs_total = atoi(optarg);
			break;
		case 'r':
			num_rows = atoi(optarg);
			break;
		case 'c':
			num_cols = atoi(optarg);
			break;
		case 'l':
			levels_back = atoi(optarg);
			break;
		case 'n':
			max_number_evaluations = atoi(optarg);
			break;

			// ES parameter
		case 'p':
			population_size = atoi(optarg);
			break;
		case 'q':
			mu = atoi(optarg);
			break;
		case 'm':
			per_cent_mutate = atof(optarg);
			break;
		case 'e':
			mutate_genes = atoi(optarg);
			break;

			// SA parameter, mutation parsed in ES section
		case '1':
			sa_parameters_global.tstart = atof(optarg);
			break;
		case '2':
			sa_parameters_global.tstop = atof(optarg);
			break;
		case '3':
			sa_parameters_global.cooling_scheme = atoi(optarg);
			break;
		case 'b':
			fitness_type = static_cast<fitness_type_def>(atoi(optarg));
			break;
		case 'd':
			explicit_neutrality = 0;
			break;

		case 'g':
			mutate_active_wires = atoi(optarg);
			break;
		case 'h':
			mutate_inactive_wires = atoi(optarg);
			break;
		case 'j':
			mutate_active_node_func = atoi(optarg);
			break;
		case 'k':
			mutate_inactive_node_func = atoi(optarg);
			break;
		case 's':
			randomize_inactive_genes = atoi(optarg);
			break;

			// Symbolic regression function
		case 'f':
			problem = atoi(optarg);
			break;

		default:
			usage(*argvv);
			exit(1);
		}
	}

	if (algorithm_type < 0 || algorithm_type > SSEA) {
		printf("select a valid algorithm, currently selected: %d\n", algorithm_type);
		exit(1);
	}

	if (fitness_type < 0 || fitness_type > CE_NODES_PROCESSED) {
		printf("fitness type selected incorrectly: %d\n", fitness_type);
		exit(1);
	}

	if (num_rows < 1 || num_cols < 1 || levels_back < 1) {
		printf("should be not smaller than zero\n");
		printf("number of rows: %d\n", num_rows);
		printf("number of columns: %d\n", num_cols);
		printf("levels back: %d\n", levels_back);
		exit(1);
	}

	levels_back = MIN(levels_back, num_cols);

	num_nodes = num_rows * num_cols;
	num_functions = 0;
	int max_arity = 0, arity = 0;

	FILE *fp;

	printf("\n********* Reading parameters defined in %s *********\n", parfile);
	fp = fopen(parfile, "r");
	if (!fp) {
		printf("Missing file: %s\n", parfile);
		exit(1);
	}

	int i;
	for (i = 0; i < MAX_NUM_FUNCTIONS; i++) {
		/* number[] holds whether the function is used or not */
		if (fscanf(fp, "%d%d%s", &number[i], &arity, (char*) &node_types[i]) == EOF) {
			printf("Formatting error while reading the parameter file: %s\n", parfile);
			exit(1);
		}
		if (number[i]) {
			allowed_functions[num_functions][0] = i;
			allowed_functions[num_functions][1] = arity;
			num_functions++;

			if (arity > max_arity)
				max_arity = arity;
		}
	}
	fclose(fp);
	printf("********* Reading parameters in %s done *********\n", parfile);

	/* each node is assigned max_arity connection genes */
	num_genes_per_node = max_arity + 1;

	/* get input data */
	read_data(datafile);

	/* if data is double calculate the output
	 * for the data points
	 */

	//calculate_data(problem);
	//#ifdef DATA_IS_DOUBLE
//	calculate_data(problem);
	//#endif
	/* calculate the perfect score and
	 for the boolean case what the bit width is

	 */
	define_perfect();

	if (population_size > MAX_NUM_CHROMOSOMES) {
		printf("Too large a population size (<= %d)\n", MAX_NUM_CHROMOSOMES);
		exit(0);
	}

	if (num_genes > MAX_NUM_GENES) {
		printf("To many genes selected (<= %d)\n", MAX_NUM_GENES);
		exit(0);
	}

	if (num_runs_total < 1) {
		puts("Number of runs of EA must be at least 1");
		exit(0);
	} else if (num_runs_total > MAX_NUM_RUNS) {
		printf("Number of runs of EA must be less than %d\n", MAX_NUM_RUNS);
		exit(0);
	}

	if (num_genes < 10) {
		puts("Number of genes/bits must be at least 10");
		exit(0);
	}

	if ((progress_report < 0) || (progress_report > 1)) {
		puts("Progress report parameter must be 0 or 1");
		exit(0);
	}

	/*if (problem == 0) {
	 puts("Function number must be greater than 0");
	 exit(0);
	 }*/

	/* assigned global constants */
	num_nodes = num_rows * num_cols;

	// set population size to 1 in case of SA
	if (algorithm_type == SIMULATED_ANNEALING || algorithm_type == LAHC || algorithm_type == RANDOM_WALK)
		population_size = 1;

	// check, if mandatory parameters are set
	//if (algorithm_type == PLUS_ES && (per_cent_mutate <= 0.0 || per_cent_mutate > 100.0 || population_size < 1 || population_size < mu)) {
	if (algorithm_type == PLUS_ES && (per_cent_mutate < 0.0 || per_cent_mutate > 100.0 || population_size < 1)) {
		printf("parameterization of evolutionary strategies incorrect\n");
		printf("mutation rate: %lf\nlambda = population size %d\nmu = number of parents %d\n", per_cent_mutate, population_size, mu);
		usage(*argvv);
	}

	if ((algorithm_type == SIMULATED_ANNEALING)
			&& (sa_parameters_global.tstart < 0.0 || sa_parameters_global.tstop < 0.0 || sa_parameters_global.tstart <= sa_parameters_global.tstop || sa_parameters_global.cooling_scheme < 1
					|| sa_parameters_global.cooling_scheme > 8))
		usage(*argvv);

	if (algorithm_type == RANDOM_WALK && (per_cent_mutate < 0.0 || per_cent_mutate > 100.0 || population_size != 1))
		usage(*argvv);

	// initializing the random number generator
	FILE *fi = fopen("/dev/urandom", "r");
	if (fgets((char*) &global_seed, sizeof(global_seed), fi) == NULL) {
		printf("Error while reading from /dev/urandom\n");
		exit(1);
	}
	fclose(fi);

	srand(global_seed);
	printf("global seed: %d\n", global_seed);

	return 1;
}

int main(int argcc, char *argvv[], char **envp) {

	parse_parameters(argcc, argvv, envp);

	write_cgp_info(stdout, argvv[0], parfile, datafile);

	double costs = run_optimizer(num_runs_total);

	printf("\n\nBest %lf\n", costs);

	return costs;
}
