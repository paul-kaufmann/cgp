#include "mutation.h"
#include "cgp.h"
#include "optimizationlog.h"

const char *mutation_type_def2str[] = { "func", "conn", "outp" };

/* calculates how many mutations to do per chromosome */
int get_num_mutant(int num_genes, double per_cent_mutate) {
	return (int) (num_genes * per_cent_mutate / 100.0);
}

/* mutate exactly one gene in a chromosome
 * return the node number that has been mutated (add num_inputs to it)
 *
 * beware:
 * active_nodes[0]       is the first input node/gene, unused here
 * active_nodes[n_i - 1] is the last input node/gene, unused here
 * active_nodes[n_i .. n_i + num_nodes - 1] are the inner nodes, good to mutate
 * active_nodes[n_i + num_nodes .. n_i + num_nodes + n_o -1] are the output nodes, good to mutate
 *
 * but:
 * choromosome[0] is active_nodes[n_i] - the first inner node
 *
 * the indexing of chomosome doesn't index the input nodes
 *
 * return zero if no mutation was done / unsuccessful, repeat mutation
 */
int mutate_a_gene(int *chromosome, int *active_nodes, int *mutation_type) {
	int which_gene, which_locus;
	int limit, limit_min;
	int col;

	while (1) { // try mutating until success
		which_gene = newrand(num_genes);
		which_locus = which_gene % num_genes_per_node;

		if (is_not_output_gene(which_gene)) {
			if (is_function_gene(which_gene, which_locus)) {
				if (num_functions == 1) /* redirect the mutation to a connection */
				{
					// seems to be an incorrect implementation, but we dont have the case of num_functions == 1
					which_locus = newrand(num_genes_per_node - 1);
					which_gene = which_gene - num_genes_per_node - 1 + which_locus;
				}

				// are we mutating active genes and mutation of active genes is disallowed? --> try another gene
				if (active_nodes[num_inputs + (which_gene / num_genes_per_node)] && !mutate_active_node_func)
					continue;

				// are we mutating inactive genes and mutation of inactive genes is disallowed? --> try another gene
				if (!active_nodes[num_inputs + (which_gene / num_genes_per_node)] && !mutate_inactive_node_func)
					continue;

				*mutation_type = FUNCTION_MUTATION;
				chromosome[which_gene] = get_function_gene();
			} else /* it is a connection gene */
			{
				// are we mutating active genes and mutation of active genes is disallowed? --> try another gene
				if (active_nodes[num_inputs + (which_gene / num_genes_per_node)] && !mutate_active_wires)
					continue;

				// are we mutating inactive genes and mutation of inactive genes is disallowed? --> try another gene
				if (!active_nodes[num_inputs + (which_gene / num_genes_per_node)] && !mutate_inactive_wires)
					continue;

				*mutation_type = CONNECTION_MUTATION;
				col = which_gene / (num_genes_per_node * num_rows);

				get_gene_limits(col, &limit_min, &limit);

				chromosome[which_gene] = get_connection_gene(limit_min, limit);
			}

			// compute node number from the function/connection gene number
			return num_inputs + (which_gene / num_genes_per_node);
		} else { /* it is an output gene */

			// are we mutating active genes and mutation of active genes is disallowed? --> try another gene
			if (active_nodes[num_inputs + num_nodes + (which_gene - (num_genes_per_node * num_nodes))] && !mutate_active_wires)
				continue;

			// are we mutating inactive genes and mutation of inactive genes is disallowed? --> try another gene
			if (!active_nodes[num_inputs + num_nodes + (which_gene - (num_genes_per_node * num_nodes))] && !mutate_inactive_wires)
				continue;

			*mutation_type = OUTPUT_MUTATION;
			chromosome[which_gene] = get_output_gene();

			// calculate the node number of an output node from its gene location
			return num_inputs + num_nodes + (which_gene - (num_genes_per_node * num_nodes));
		}
	}
}

/*
 modification of the original get_nodes_to_process: do not compact active node numbers into a
 new array but return the original node_used array, so one can easily check with node index, if
 the node is active

 calculates the addresses of nodes used
 and stores them in node_to_process
 */
void get_used_nodes(int *chromosome, int *nodes_used) {
	int i, j, index;
//	int num_nodes_to_process;
	int node_genes[MAX_NUM_GENES_PER_NODE];
	//int nodes_used[MAX_NUM_NODES_PLUS_INPUTS];
	int max_size_node_used;

	max_size_node_used = num_nodes + num_inputs;

	/* say all nodes not used */
	for (i = 0; i < max_size_node_used; i++)
		nodes_used[i] = 0;

	/* all output nodes are active */
	for (i = max_size_node_used; i < max_size_node_used + num_outputs; i++)
		nodes_used[i] = 1;

	/* all the nodes whose output is given by the output genes are active */
	/* last num_outputs genes are all output genes */
	for (i = num_genes - num_outputs; i < num_genes; i++)
		nodes_used[chromosome[i]] = 1;

	for (i = max_size_node_used - 1; i >= num_inputs; i--) {
		if (nodes_used[i]) {
			/* get input addresses and type of this gate */
			index = num_genes_per_node * (i - num_inputs);

			/* copy genes for node into array node_genes */
			for (j = 0; j < num_genes_per_node; j++)
				node_genes[j] = chromosome[index + j];

			/* each function has an arity stored in
			 allowed_functions[][1].
			 Find the nodes whose data is used
			 */

			for (j = 0; j < get_arity(node_genes[num_genes_per_node - 1]); j++)
				nodes_used[node_genes[j]] = 1;
		}
	}
}

/**
 * randomize genes of all (inner) inactive nodes
 */
void randomize_inactive_nodes(int *chromosome, int *active_nodes) {
	int which_locus;
	int limit, limit_min;
	int col;

	// def: num_genes_per_node = max_arity + 1;
	// def: num_nodes = num_rows * num_cols;
	// def: num_genes = num_genes_per_node * num_nodes + num_outputs;

	for (int which_gene = 0; which_gene < (num_genes_per_node * num_nodes); which_gene++) {

		// are we trying to mutate an active gene? yes --> skip mutating
		if (active_nodes[num_inputs + (which_gene / num_genes_per_node)])
			continue;

		// randomize gene = mutate gene
		which_locus = which_gene % num_genes_per_node;
		if (is_function_gene(which_gene, which_locus))
			chromosome[which_gene] = get_function_gene();
		else { /* it is a connection gene */
			col = which_gene / (num_genes_per_node * num_rows);
			get_gene_limits(col, &limit_min, &limit);
			chromosome[which_gene] = get_connection_gene(limit_min, limit);
		}
	}
}

/**
 *  carry out num_mutations mutations on the chromosome
 */
void mutate_chromosome(int *chromosome, int num_mutations) {
	int i, j;
	int mutation_type;

	//printf("mutating %d genes\n", num_mutations);

	if (num_mutations <= 0) {
		printf("mutation rate to small: %f\n", per_cent_mutate);
		exit(1);
	}

	int active_nodes[MAX_NUM_NODES_PLUS_INPUTS + MAX_NUM_OUTPUTS];
	get_used_nodes(chromosome, active_nodes);

	for (i = 0; i < num_mutations; i++) {
		j = mutate_a_gene(chromosome, active_nodes, &mutation_type);
		if (active_nodes[j] && mutation_type != FUNCTION_MUTATION)
			get_used_nodes(chromosome, active_nodes);
	}

	// randomize genes of all inactive nodes, if requested
	if (randomize_inactive_genes)
		randomize_inactive_nodes(chromosome, active_nodes);
}

/**
 * mutate the chromosome according to the mutation rate
 */
void mutate_chromosome(int *chromosome) {

	int num_mutations = get_num_mutant(num_genes, per_cent_mutate);

	// printf("mutating according to the mutation rate=%f gives for %d mutated genes\n", per_cent_mutate, num_mutations);

	mutate_chromosome(chromosome, num_mutations);
}

/* mutate at least num_mutations active genes in a chromosome */
void mutate_active_chromosome(int *chromosome, int num_mutations, int evaluation_no, optimization_log *log) {

	/* 1. generate a map of active genes of a chromosome
	 * 2. mutate random genes until num_mutations active genes have been hit and mutated
	 * 3. if an active node has been mutated, re-identify active genes
	 */
	int active_nodes[MAX_NUM_NODES_PLUS_INPUTS + MAX_NUM_OUTPUTS];
	int num_active_mutations = 0;
	int i;
	int j;
	int debugs = 0;
	int mutation_type;

	get_used_nodes(chromosome, active_nodes);

	if (debugs) {
		printf("active nodes: ");
		for (i = 0; i < num_inputs + num_nodes + num_outputs; i++)
			if (active_nodes[i])
				printf("%d ", i);
		printf("\nmutations: %d ", num_mutations);
	}

// a test here: mutate until hitting an active node. stop after mutating a single active node.
// investigate, how many active nodes have to be mutated using irace. compare amount of genetic
// material mutated in the original mutation operator and this new one
	if (num_mutations < 1)
		num_mutations = 1;

	while (num_active_mutations < num_mutations) {
		j = mutate_a_gene(chromosome, active_nodes, &mutation_type);
		log->addrecord(evaluation_no, active_nodes[j], mutation_type, -9);

		if (active_nodes[j]) {

			if (debugs)
				printf("%d(a) --> %d\n", j, num_active_mutations);

			num_active_mutations++;
			if (mutation_type != FUNCTION_MUTATION)
				get_used_nodes(chromosome, active_nodes);

			if (debugs) {
				printf("active nodes: ");
				for (i = 0; i < num_inputs + num_nodes + num_outputs; i++)
					if (active_nodes[i])
						printf("%d ", i);
				printf("\n");
			}
		} else {
			if (debugs)
				printf("%d(i) ", j);
		}
	}

// randomize genes of all inactive nodes, if requested
	if (randomize_inactive_genes)
		randomize_inactive_nodes(chromosome, active_nodes);
}

/**
 * implements a mutation operator
 *
 * if per_cent_mutate <= 0.0 then use Goldmann & Punch mutation of one active gene, otherwise use the regular CGP mutation
 *
 */
void mutate(int *chromosome, int evaluation_no, optimization_log *log) {

	if (mutate_genes) // mutate per_cent_mutate genes of a genotype, regardless of active/inactive status
		mutate_chromosome(chromosome, (int) per_cent_mutate);
	else
		mutate_chromosome(chromosome);
}

