/*
 * optimizationlog.h
 *
 *  Created on: Nov 29, 2018
 *      Author: paul
 */

#ifndef OPTIMIZATIONLOG_H_
#define OPTIMIZATIONLOG_H_

#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

using namespace std;

extern const char *mutation_type_def2str[];

class log_line {
public:
	long evaluation_no;
	int active;
	int mutation_type;
	int impact_on_fitness;
	float fq;

	log_line(long evaliation_no, int active, int mutation_type, int impact_on_fitness) {
		this->evaluation_no = evaliation_no;
		this->active = active;
		this->mutation_type = mutation_type;
		this->impact_on_fitness = impact_on_fitness;
		this->fq = -1;
	}

	void print() {
		cout << evaluation_no << " " << active << " " << mutation_type << " " << impact_on_fitness;
	}
};

class optimization_log {

	vector<log_line*> log;

public:
	optimization_log();

	~optimization_log() {
		for (int i = 0; i < log.size(); ++i)
			delete log[i];
		log.clear();
	}

	void addrecord(long evaluation_no, int active, int mutation_type, int fitness_trend) {
		log.push_back(new log_line(evaluation_no, active, mutation_type, fitness_trend));
	}

	int* getLastImpact() {
		return &log.back()->impact_on_fitness;
	}

	float* getLastFQ() {
		return &log.back()->fq;
	}

	int* getLastMutationType() {
		return &log.back()->mutation_type;
	}
	void print() {
		std::cout << "\n\nlog start [eval_no mutation_no mutated_active/inactive_gene node_type impact_on_fitness]\n";
		for (vector<log_line*>::iterator it = log.begin(); it != log.end(); ++it) {
			if ((*it)->impact_on_fitness == 1)
				cout << setw(10) << (*it)->evaluation_no << " " << setw(12) << (it - log.begin()) << " " << ((*it)->active ? "  active" : "inactive") << " "
						<< mutation_type_def2str[(*it)->mutation_type] << " " << (*it)->fq << " " << (*it)->impact_on_fitness << endl;
		}
		std::cout << "log end\n";
	}
}
;

#endif /* OPTIMIZATIONLOG_H_ */
