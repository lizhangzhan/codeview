#include "sparseLeastSquares.h"

#include <fstream>
#include <sstream>
#include <string>

using namespace std;

SparseLeastSquaresProblem::SparseLeastSquaresProblem(const char* libsvmFilename) {
	ifstream matfile(libsvmFilename);
	if (!matfile.good()) {
		cerr << "error opening matrix file " << libsvmFilename << endl;
		exit(1);
	}

	
	bool intercept = true;

	n=0;
	numFeats=0;
	instance_starts.push_back(0); 

	string s;
	
	while (getline(matfile, s)) {
		n++;

		const char *ch = s.c_str();
		int len = strlen(ch);
		// target value
		float b;
		vector<size_t> is;
		vector<float> vs;

		if (intercept) {
			is.push_back(0);
			vs.push_back(1);
		}

		bool first = true;
		while (*ch!=0) {
			while (isspace(*ch) && *ch!=0) ch++; // skip white space
			string field;
			while (!isspace(*ch) && *ch!=0) field.append(1, *ch++);
			if (*ch == 0) break;
			if (first) {
				b = (float) atof(field.c_str());
				first = false;
			}
			else {
				size_t j=0;
				while (j<field.length() && field[j]!=':') j++;
				field[j] = 0;
				size_t idx = (size_t) atoi(&field[0]);
				float val = (float) atof(&field[j+1]);
				numFeats = max(idx, numFeats);
				if (!intercept) {
					is.push_back(idx-1); // index:value, (input) index starts from 1
				}
				else {
					is.push_back(idx);
				}
				vs.push_back(val);
			}
		}

		for (size_t i=0; i<is.size(); i++) {
			indices.push_back(is[i]);
			values.push_back(vs[i]);
		}
		instance_starts.push_back(indices.size());
		target.push_back(b);
	}

	if (intercept) numFeats++;

	matfile.close();
}


/**
 * 1/2 * sum_i [wx_i - b_i]^2 + 1/2 * lambda * ||w||
 *
 */

double SparseLeastSquaresObjective::Eval(const DblVec& input, DblVec& gradient) {
	static DblVec temp(problem.n);

	if (input.size() != problem.numFeats) {
		cerr << "Error: input is not the correct size." << endl;
		exit(1);
	}

	for (size_t i=0; i<problem.n; i++) {
		temp[i] = -problem.target[i];

		for (size_t j=problem.instance_starts[i]; j < problem.instance_starts[i+1]; j++) {
			size_t idx = problem.indices[j];
			double val = (double) problem.values[j];
			temp[i] += input[idx] * val;
		}
	}
	
	double value = 0.0;
	for (size_t j=0; j<problem.numFeats; j++) {
		value += input[j] * input[j] * l2weight;
		gradient[j] = l2weight * input[j];
	}

	for (size_t i=0; i<problem.n; i++) {
		if (temp[i] == 0) continue;

		value += temp[i] * temp[i];

		for (size_t j=problem.instance_starts[i]; j<problem.instance_starts[i+1]; j++) {
			size_t idx = problem.indices[j];
			double val = (double) problem.values[j];
			gradient[idx] += val * temp[i];
		}

	}

	return 0.5 * value + 1.0;
}
