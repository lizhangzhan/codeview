#pragma once

#include <fstream>
#include <string>

#include "OWLQN.h"

struct SparseLeastSquaresObjective;

class SparseLeastSquaresProblem {
	std::deque<size_t> indices;
	std::deque<float> values;
	std::deque<size_t> instance_starts;

	std::vector<float> target;
	size_t n, numFeats;

	friend struct SparseLeastSquaresObjective;

public:
	//SparseLeastSquaresProblem(size_t m, size_t n) : Amat(m * n), b(m), m(m), n(n) { }

	SparseLeastSquaresProblem(const char* libsvmfile);

	size_t NumFeats() const { return numFeats; }
	size_t NumInstances() const { return n; }
};

struct SparseLeastSquaresObjective : public DifferentiableFunction {
	const SparseLeastSquaresProblem& problem;
	const double l2weight;

	SparseLeastSquaresObjective(const SparseLeastSquaresProblem& p, double l2weight = 0) : problem(p), l2weight(l2weight) { }

	double Eval(const DblVec& input, DblVec& gradient);
};
