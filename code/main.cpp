#include "Util.h"
#include <iomanip>
#include <boost/bind.hpp>
#include <algorithm>


class mystruct0 {
public :
	mystruct0
		(
		string i,
		std::vector<double> r
		) : _i(i), _r(r) {;}

	string _i;
	std::vector<double> _r;
};

class mystruct {
public:
	mystruct 
		(string mf, double _52wk, double _26wk, double _13wk) : _sym (mf), _52wk (_52wk), _26wk (_26wk), _13wk (_13wk) {}

	string _sym;
	double _52wk;
	double _26wk;
	double _13wk;
};

struct sort_52wk {
	bool operator()(const mystruct &left, const mystruct &right) {
		return left._52wk > right._52wk;
	}
};

struct sort_26wk {
	bool operator()(const mystruct &left, const mystruct &right) {
		return left._26wk > right._26wk;
	}
};

struct sort_13wk {
	bool operator()(const mystruct &left, const mystruct &right) {
		return left._13wk > right._13wk;
	}
};

void
	func0
	(
	void
	)
{
	ifstream ifile("mutual_funds.txt");
	
	string mf;
	unsigned num_mfs = 0;
	std::vector<mystruct> ret_vector;

	while (ifile.good()) {
		std::vector<multivariate_regression_results_t> results;
		getline(ifile, mf);

		string file_str = download_market_data(mf, MARKET_DATA_WEEKLY);

		std::vector<double> mf_prices;
		populate_price_vector(
			file_str,
			NUM_WEEKS_IN_A_YEAR,
			mf_prices);

		if (mf_prices.size() < NUM_WEEKS_IN_A_YEAR)
			continue;

		double _52wk = 100. * (mf_prices[0] - mf_prices[NUM_WEEKS_IN_A_YEAR-1])/mf_prices[NUM_WEEKS_IN_A_YEAR-1];
		double _26wk = 100. * (mf_prices[0] - mf_prices[2*NUM_WEEKS_IN_A_QUARTER-1])/mf_prices[2*NUM_WEEKS_IN_A_QUARTER-1];
		double _13wk = 100. * (mf_prices[0] - mf_prices[NUM_WEEKS_IN_A_QUARTER-1])/mf_prices[NUM_WEEKS_IN_A_QUARTER-1];
		
		ret_vector.push_back(mystruct(mf, _52wk, _26wk, _13wk));
		if (num_mfs++ > 1000)
			break;
	}

	std::sort(
		ret_vector.begin(), 
		ret_vector.end(), 
		sort_52wk());

	ofstream ofile("correlation.txt");

	for (auto iter = ret_vector.begin(); iter != ret_vector.end(); ++iter) {
		cout << iter->_sym << endl;
		ofile << iter->_sym << ":" << iter->_52wk << endl;

		std::vector<multivariate_regression_results_t> results;
		std::vector<string> etfs;
		etfs.push_back("SPY");
		replicate_mutual_fund_with_etfs(
			results,
			iter->_sym,
			etfs,
			1,
			NUM_WEEKS_IN_A_YEAR,
			ofile);
		etfs.push_back("EFA");
		replicate_mutual_fund_with_etfs(
			results,
			iter->_sym,
			etfs,
			1,
			NUM_WEEKS_IN_A_YEAR,
			ofile);
		etfs.push_back("VWO");
		replicate_mutual_fund_with_etfs(
			results,
			iter->_sym,
			etfs,
			1,
			NUM_WEEKS_IN_A_YEAR,
			ofile);
		etfs.push_back("VEA");
		replicate_mutual_fund_with_etfs(
			results,
			iter->_sym,
			etfs,
			1,
			NUM_WEEKS_IN_A_YEAR,
			ofile);
		etfs.push_back("BND");
		replicate_mutual_fund_with_etfs(
			results,
			iter->_sym,
			etfs,
			1,
			NUM_WEEKS_IN_A_YEAR,
			ofile);
		etfs.push_back("VNQ");
		replicate_mutual_fund_with_etfs(
			results,
			iter->_sym,
			etfs,
			1,
			NUM_WEEKS_IN_A_YEAR,
			ofile);

		ofile << "----------------------------------------" << endl;

	}

}

void
	func1
	(
	void
	)
{
	matrix<double> etf_returns(100, NUM_WEEKS_IN_A_YEAR);
	ifstream ifile("etfs.txt");
	unsigned etf_i = 0;

	string line, etf_symbol;
	while (ifile.good()) {
		getline(ifile, line);
		if (line.empty())
			break;

		etf_symbol = get_csv_line_cell_as_string(line, 0);

		string filename = download_market_data(etf_symbol, MARKET_DATA_WEEKLY);
		populate_return_vector(
			filename,
			NUM_WEEKS_IN_A_YEAR,
			etf_i,
			etf_returns);

		++etf_i;
	}

	ifstream mf_file("mutual_funds.txt");
	
	string mf;
	unsigned num_mfs = 0;
	std::vector<mystruct> ret_vector;

	while (mf_file.good()) {
		std::vector<multivariate_regression_results_t> results;
		getline(mf_file, mf);

		string file_str = download_market_data(mf, MARKET_DATA_WEEKLY);

		std::vector<double> mf_prices;
		populate_price_vector(
			file_str,
			NUM_WEEKS_IN_A_YEAR,
			mf_prices);

		if (mf_prices.size() < NUM_WEEKS_IN_A_YEAR)
			continue;

		double _52wk = 100. * (mf_prices[0] - mf_prices[NUM_WEEKS_IN_A_YEAR-1])/mf_prices[NUM_WEEKS_IN_A_YEAR-1];
		double _26wk = 100. * (mf_prices[0] - mf_prices[2*NUM_WEEKS_IN_A_QUARTER-1])/mf_prices[2*NUM_WEEKS_IN_A_QUARTER-1];
		double _13wk = 100. * (mf_prices[0] - mf_prices[NUM_WEEKS_IN_A_QUARTER-1])/mf_prices[NUM_WEEKS_IN_A_QUARTER-1];
		
		ret_vector.push_back(mystruct(mf, _52wk, _26wk, _13wk));
		if (num_mfs++ > 1000)
			break;
	}

	std::sort(
		ret_vector.begin(), 
		ret_vector.end(), 
		sort_52wk());

}

void
	print_rets
	(
	string i,
	std::vector<double> r
	)
{
	for(auto iter = r.begin(); iter != r.end(); ++iter)
		cout << fixed << setprecision(4) << (*iter*100.f) << "%  ";
	cout << 
		i << endl;
}

int
	main
	(
	int argc,
	char ** argv
	)
{
	ifstream mf_file("mutual_funds.txt");
	ofstream results("results.txt");
	string mf;
	std::vector<double> rets(8);
	int num = 30;

	std::vector<mystruct0> rets;

	mf = "SPY";
	for (unsigned q = 0; q < 8; ++q) {
		rets[q] = get_return(
			mf,
			MARKET_DATA_WEEKLY,
			NUM_WEEKS_IN_A_QUARTER,
			q*NUM_WEEKS_IN_A_QUARTER);
	}

	mystruct a(mf, rets);

	rets.push_back(a);
	print_rets(mf, rets);

	try {
		while (mf_file.good()) {
			getline(mf_file, mf);
			if (mf.empty()) break;

			for (unsigned q = 0; q < 8; ++q) {
				rets[q] = get_return(
					mf,
					MARKET_DATA_WEEKLY,
					NUM_WEEKS_IN_A_QUARTER,
					q*NUM_WEEKS_IN_A_QUARTER);
			}

			print_rets(mf, rets);

			if (num-- == 0) break;
		}
	}
	catch (exception e) {

	}

	mf_file.close();

	system("pause");
	return 0;
}