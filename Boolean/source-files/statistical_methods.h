/*
 * statistical_methods.h
 *
 *  Created on: Mar 23, 2018
 *      Author: paul
 */

#ifndef STATISTICAL_METHODS_H_
#define STATISTICAL_METHODS_H_

#include "statistical_methods.h"
#include "cgp.h"
#include "mutation.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

enum fitness_type_def {
	MEDIAN_EVALUATIONS, MEDIAN_FITNESS, NODES_PROCESSED, CE_FITNESS, CE_NODES_PROCESSED
};

extern enum fitness_type_def fitness_type;

void print_avg_quartiles_int(FILE *outfile, const char *text, int* elements, int n);
void print_avg_quartiles_double(FILE *outfile, const char *text, double* elements, int n);
double get_standard_deviation(int num_items, double average, double *items);
double get_standard_deviation_int(int num_items, double average, int *items);
//double calc_CE(int evaluations_for_best[], int num_runs_total, int max_number_of_evaluations, int *at_evaluation);
void report_final_results(FILE* outfile, int num_runs_total, double *fitnesses, long long *fitness_evals_to_perfect, int *active_nodes_in_run, long *nodes_processed);

template<typename T> int partition(T a[], int l, int r) {
	T pivot, t;
	int i, j;

	pivot = a[l];
	i = l;
	j = r + 1;

	while (1) {
		do
			++i;
		while (a[i] <= pivot && i <= r);
		do
			--j;
		while (a[j] > pivot);
		if (i >= j)
			break;
		t = a[i];
		a[i] = a[j];
		a[j] = t;
	}
	t = a[l];
	a[l] = a[j];
	a[j] = t;
	return j;
}

template<typename T> void quickSort(T a[], int l, int r) {
	int j;

	if (l < r) {
		// divide and conquer
		j = partition(a, l, r);
		quickSort(a, l, j - 1);
		quickSort(a, j + 1, r);
	}

}

template<class T>
T min(T *elements, int n) {
	T min_element = elements[0];
	for (int i = 1; i < n; i++)
		if (elements[i] < min_element)
			min_element = elements[i];
	return min_element;
}

template<class T>
T max(T *elements, int n) {
	T min_element = elements[0];
	for (int i = 1; i < n; i++)
		if (elements[i] > min_element)
			min_element = elements[i];
	return min_element;
}

template<class T>
double get_standard_deviation(int num_items, double average, T *items) {
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

template<class T>
void quartiles(T* values, int n, T *q1, T *q2, T *q3) {
	T sorted[n];
	memcpy(sorted, values, n * sizeof(T));
	quickSort(sorted, 0, n - 1);
	*q1 = sorted[(n - 1) / 4];
	*q2 = sorted[(n - 1) / 2];
	*q3 = sorted[(3 * (n - 1)) / 4];
}

template<typename T> T average(int num_items, T *items) {
	int i;
	double av = 0.0;

	av = 0.0;
	for (i = 0; i < num_items; i++)
		av = av + items[i];

	av = av / ((double) num_items);

	return av;
}

template<typename T> void print_avg_quartiles(FILE *outfile, const char *text, T* elements, int n) {
	double avg = average(n, elements);
	double st_dev = get_standard_deviation(n, avg, elements);

	T q1, q2, q3; // quartiles
	quartiles(elements, n, &q1, &q2, &q3);
	cout << text << "avg=" << left << setw(16) << setprecision(13) << avg << " std=" << left << setw(16) << setprecision(13) << st_dev << " q1=" << left << setw(16) << setprecision(13) << q1 << " q2=" << left << setw(16) << setprecision(13) << q2
			<< " q3=" << left << setw(16) << setprecision(13) << q3 << endl;
}

template<typename T> double calc_CE(T *evaluations_for_best, int num_runs_total, T max_number_of_evaluations, T *at_evaluation) {

	int succeeded;
	T i;
	int j;
	double P, R;
	double z = 0.99;
	double currentI;
	double CE;

	// reset return values
	CE = FLT_MAX;
	*at_evaluation = -1;

	for (i = 0; i < max_number_of_evaluations; i++) {

		// find number of successful runs for the current generation i
		for (succeeded = j = 0; j < num_runs_total; j++)
			if (evaluations_for_best[j] < i)
				succeeded++;

		// compute CE for the current generation; compute only if enough statistical data has been collected
		// I define "enough" as having at least 10 successful runs ;)
		if (succeeded >= 10) {
			P = ((float) succeeded) / ((float) num_runs_total);
			R = P < 1.0 ? log(1.0 - z) / log(1.0 - P) : 1.0;
			currentI = (i + 1) * R;
			if (currentI < CE) {
				CE = currentI;
				*at_evaluation = i;
			}
		}
	}
	if (CE == FLT_MAX) {
		CE = NAN;
		*at_evaluation = -1;
	}
	return CE;
}

#endif /* STATISTICAL_METHODS_H_ */
