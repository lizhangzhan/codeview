// TODO: use pointers in class through the whole program 
// TODO: inspect possible levels, and give data summarize

#define _CRT_SECURE_NO_WARNINGS

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <deque>
#include <algorithm>
#include <math.h>
#include <ctime>
#include <limits>

using namespace std;

void error(string s) {
	cerr << s << endl;
	exit(1);
}


inline bool eq(double a, double b)
{
    return fabs(a - b) < 1e-8;
}

#define STR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str() 

class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str,line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }
    private:
        std::vector<std::string>    m_data;
};

std::istream& operator>>(std::istream& str,CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   


/*
 * DATA STRUCTURES 
 */

enum FeatureType {F_Numeric, F_Category, F_Order};

struct SVec {
	vector<size_t> idx;
	vector<double> vals;
	string name;
	FeatureType type;
	size_t dim;

	double _l2norm;

	SVec() { 
		unsigned long nan[2]={0xffffffff, 0x7fffffff};
		_l2norm = *( double* )nan;
	}

	SVec* copy() {
		SVec *c = new SVec;
		*c = *this;
		return c;
	}

	void calcl2norm() {
		double sum=0;
		for (size_t i=0; i<vals.size(); i++) sum += vals[i]*vals[i];
		_l2norm = sqrt(sum / dim);
	}

	void assign(const vector<double>& vec) {
		dim = vec.size();
		idx.clear();
		vals.clear();
		for (size_t i=0; i<dim; i++) {
			idx.push_back(i);
			vals.push_back(vec[i]);
		}
	}

	double l2norm() const {
		if (_isnan(_l2norm)) error("l2norm nan");
		return _l2norm;
	}

	double sum() const {
		double sum=0;
		for (size_t i=0; i<vals.size(); i++) sum += vals[i];
		return sum;
	}

	double reflectiveCorrelation(const SVec& other) const {
		double xy = mul(other).sum() / dim;
		double xnorm = l2norm();
		double ynorm = other.l2norm();
		if (eq(xnorm, 0) || eq(ynorm, 0)) return 0;
		else return xy / xnorm / ynorm;
	}

	SVec add(const SVec& other, double multiplier=1) const {
		SVec res;
		if (dim != other.dim) error("SVec, add");
		res.name = name + "+" + other.name;
		res.type = type;
		res.dim  = dim;


		size_t i = 0, j = 0;
		while (i < idx.size() && j < other.idx.size()) {
			if (idx[i] == other.idx[j]) {
				res.idx.push_back(idx[i]);
				res.vals.push_back(vals[i] + multiplier * other.vals[j]);
				i++;
				j++;
			} 
			else if (idx[i] < other.idx[j]) {
				res.idx.push_back(idx[i]);
				res.vals.push_back(vals[i]);
				i++;
			}
			else {
				res.idx.push_back(other.idx[j]);
				res.vals.push_back(other.vals[j] * multiplier);
				j++;
			}
		}
		return res;
	}

	SVec mul(const SVec& other) const {
		SVec res;
		if (dim != other.dim) error("SVec, add");
		res.name = name + "*" + other.name;
		res.type = type;
		res.dim  = dim;

		if (isSparse() && other.isSparse()) {
			size_t i = 0, j = 0;
			while (i < idx.size() && j < other.idx.size()) {
				if (idx[i] == other.idx[j]) {
					res.idx.push_back(idx[i]);
					res.vals.push_back(vals[i] * other.vals[j]);
					i++;
					j++;
				}
				else if (idx[i] < other.idx[j]) i++;
				else j++;
			}
		}
		else if (!isSparse() && other.isSparse()) {
			res.idx = other.idx;
			res.vals = other.vals;
			for (size_t i=0; i<res.idx.size(); i++)
				res.vals[i] *= vals[res.idx[i]];
		}
		else if (isSparse() && !other.isSparse()) {
			res.idx = idx;
			res.vals = vals;
			for (size_t i=0; i<res.idx.size(); i++)
				res.vals[i] *= other.vals[res.idx[i]];
		}
		else { // both dense 
			res.idx  = idx;
			res.vals.resize(dim);
			for (size_t i=0; i<vals.size(); i++) 
				res.vals[i] = vals[i] * other.vals[i];
		}
		return res;
	}

	bool isSparse() const {
		return dim>vals.size();
	}
};

struct Dataset {
	// column pointers
	SVec *y;
	vector<SVec *> oriX;
	vector<SVec *> splitX;
	vector<SVec *> crossX;

	// TODO: add row pointer
};

struct Interval {
	// (l,h] or [l,h] intervals
	double low, high;
	string name;
	bool close_on_low;
	Interval(): close_on_low(false) {}
};


struct IntervalTransformer {
	vector<Interval> intervals;
	const SVec *vec;
	/* 
	 * apply the transformer on vec, and add the news cols to cols, and update ds
	 */
	void apply(Dataset &ds, vector<SVec*>& cols) {
		for (size_t i=0; i<intervals.size(); i++) {
			SVec *split = new SVec();
			split->dim = vec->dim;
			split->type = vec->type;
			split->name = vec->name + "_" + intervals[i].name;

			for (size_t j=0; j<vec->vals.size(); j++) {
				if (intervals[i].close_on_low && vec->vals[j] >= intervals[i].low || vec->vals[j] > intervals[i].low) {
					if (intervals[i].high>=vec->vals[j]) {
						split->idx.push_back(vec->idx[j]);
						split->vals.push_back(vec->vals[j] - intervals[i].low + 1); // make the new range to be [1,xx], and leave all others as 0
					}
				}
			}

			cols.push_back(split);
			ds.splitX.push_back(cols[cols.size()-1]);
		}
	}

	
	static IntervalTransformer GetEvenIntervals(const SVec &vec, int parts) {
		vector<double> probs;
		probs.push_back(0);
		for (int i=1; i<=parts; i++) {
			probs.push_back(1.0 / parts * i);
		}
		return GetProbIntervals(vec, probs);
	}

	static IntervalTransformer GetProbIntervals(const SVec &vec, const vector<double> &probs) {
		
		size_t n_vals = vec.vals.size();
		
		vector<int> p;
		//require probs's first to be 0, last to be 1.0
		//if (!eq(probs[0], 0)) p.push_back(0);
		p.push_back(0);
		for (size_t i=1; i<probs.size(); i++) p.push_back((int)(n_vals * probs[i]) - 1);
		//if (!eq(probs[probs.size()-1], 1.0)) p.push_back(n_vals - 1);

		vector<double> vals(vec.vals.begin(), vec.vals.end());
		sort(vals.begin(), vals.end()); // don't count 0 vals in sparse vec

		IntervalTransformer transformer;

		char s[1000];
			
		Interval interval;
		interval.close_on_low = true;
		interval.low = vals[p[0]];
		interval.high = vals[p[1]];
		sprintf(s, "_%.3fto%.3f", probs[0], probs[1]);
		interval.name = string(s);
		
		transformer.intervals.push_back(interval);
		transformer.vec = &vec;
		for (size_t i=1; i<p.size()-1; i++) {
			interval.close_on_low = false;
			interval.low = vals[p[i]];
			interval.high = vals[p[i+1]];
			sprintf(s, "_%.3fto%.3f", probs[i], probs[i+1]);
			interval.name = string(s);

			transformer.intervals.push_back(interval);

		}

		return transformer;
	}

	//TODO: probs pairs
	//static IntervalTransformer GetProbIntervals(const SVec &vec, const vector<pair<double, double> > &probs) {}
};



void addCorssings(Dataset &ds, double topx=0.3) {
	vector<pair<double, pair<SVec*,SVec*> > > corr;
	vector<SVec*> vecs;
	vecs.insert(vecs.end(), ds.oriX.begin(), ds.oriX.end());
	vecs.insert(vecs.end(), ds.splitX.begin(), ds.splitX.end());
	
	ds.y->calcl2norm(); // pre calculate norm of y

	for (size_t i=0; i<vecs.size(); i++) {
		for (size_t j=i; j<vecs.size(); j++) {
			SVec m = vecs[i]->mul(*vecs[j]);
			m.calcl2norm();
			double c = ds.y->reflectiveCorrelation(m);
			corr.push_back(make_pair(fabs(c), make_pair(vecs[i], vecs[j])));
		}
	}

	size_t select = (size_t) (topx * corr.size());
	sort(corr.begin(), corr.end());
	reverse(corr.begin(), corr.end());
	
	for (size_t i=0; i<select; i++) {
		SVec m = corr[i].second.first->mul(*corr[i].second.second);
		ds.crossX.push_back(m.copy());
	}
}


/*
 * Global Data 
 */
vector<SVec*> cols;
Dataset dataset;

void read_csv(ifstream &file, bool header=false) {
	CSVRow row;
	vector<string> colnames;
	bool first = true;
	size_t dim = 0;
	vector<vector<double> > dense_cols;

	if (header) {
		file >> row;
		for (size_t i=0; i<row.size(); i++) colnames.push_back(row[i]);
	}

	while (file >> row) {
		if (first) {
			first = false;
			dim = row.size();
			dense_cols.resize(dim);
		}
		for (int i=0; i<dim; i++) {
			dense_cols[i].push_back(atof(row[i].c_str()));
		}
	}

	// update global variable: cols
	
	for (int i=0; i<dim; i++) {
		cols.push_back(new SVec());
		cols[i]->assign(dense_cols[i]);
		if (header) cols[i]->name = colnames[i];
		else cols[i]->name = "V" + STR(i);
		cols[i]->type = F_Numeric;
	}

	// update global variable: dataset
	dataset.y = cols[0];
	for (int i=1; i<dim; i++) {
		dataset.oriX.push_back(cols[i]);
	}
}

void write_libsvm(const char *name, Dataset& ds, size_t baseIdx=1)
{

	size_t nOri = ds.oriX.size();
	size_t nOriSplit = nOri + ds.splitX.size();
	size_t nOriSplitCross = nOriSplit + ds.crossX.size();

	size_t nrow = ds.y->dim;;
	vector<SVec*> rows(nrow);

	for (size_t i=0; i<nrow; i++) {
		rows[i] = new SVec();
		rows[i]->dim = nOriSplitCross;
		rows[i]->name = "row" + STR(i);
		rows[i]->type = F_Numeric;
	}


	for (size_t i=0; i<ds.oriX.size(); i++) { // each col
		for (size_t j=0; j<ds.oriX[i]->idx.size(); j++) { // each non-zero row
			size_t ii = ds.oriX[i]->idx[j];
			double vv = ds.oriX[i]->vals[j];
			rows[ii]->idx.push_back(i);
			rows[ii]->vals.push_back(vv);
		}
	}

	for (size_t i=0; i<ds.splitX.size(); i++) { // each col
		for (size_t j=0; j<ds.splitX[i]->idx.size(); j++) { // each non-zero row
			size_t ii = ds.splitX[i]->idx[j];
			double vv = ds.splitX[i]->vals[j];
			rows[ii]->idx.push_back(i + nOri);
			rows[ii]->vals.push_back(vv);
		}
	}
	
	for (size_t i=0; i<ds.crossX.size(); i++) { // each col
		for (size_t j=0; j<ds.crossX[i]->idx.size(); j++) { // each non-zero row
			size_t ii = ds.crossX[i]->idx[j];
			double vv = ds.crossX[i]->vals[j];
			rows[ii]->idx.push_back(i + nOriSplit);
			rows[ii]->vals.push_back(vv);
		}
	}


	FILE *fp = fopen(name, "w");
	// output each row
	for (size_t i=0; i<nrow; i++) {
		fprintf(fp, "%f", ds.y->vals[i]);
		for (size_t j=0; j<rows[i]->idx.size(); j++) {
			fprintf(fp, " %d:%f", (baseIdx+rows[i]->idx[j]), rows[i]->vals[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	for (size_t i=0; i<nrow; i++) {
		delete rows[i];
	}
}


int main()
{
	time_t tstart, tend; 
	tstart = time(0);

	//D:\CodeMesh\citadel\data\libsvm\mg_scale

    std::ifstream csvfile("D:\\CodeMesh\\citadel\\data\\libsvm\\mg_scale.csv");

	read_csv(csvfile, false);

	for (size_t i=0; i<dataset.oriX.size(); i++) {
		IntervalTransformer it = IntervalTransformer::GetEvenIntervals(*dataset.oriX[i], 3);
		it.apply(dataset, cols);
	}
	addCorssings(dataset, 0.1);

	write_libsvm("D:\\CodeMesh\\citadel\\data\\libsvm\\mg_scale.split", dataset);


	tend = time(0); 
	cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;
	return 0;
}

