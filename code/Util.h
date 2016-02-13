#ifndef __UTIL_H__
#define __UTIL_H__

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cstdarg>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <limits>
#include <sstream>
#include <ctime>

#include <Windows.h>
#include <tchar.h>
#include <Urlmon.h>
#include <windef.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_statistics_double.h>


#pragma comment(lib, "urlmon.lib")

using namespace boost::numeric::ublas;
using namespace std;

class Date {
public:
	Date 
		(
		int year, 
		int month, 
		int day
		);

	int 
		get_year  
		(
		void
		);

	int 
		get_month 
		(
		void
		);

	int 
		get_day
		(
		void
		);

private:
	int _year;
	int _month;
	int _day;
};

// Custom Instrument Type
class Instrument
{
public:
	Instrument 
		(
		unsigned id,
		string symbol, 
		string name
		);

	Instrument
		(
		void
		);

	Instrument
		(
		const Instrument &rhs
		);

	unsigned 
		get_id 
		(
		void
		);

	string 
		get_symbol 
		(
		void
		);

	string 
		get_name 
		(
		void
		);

	void 
		pretty_print 
		(
		void
		);

	bool
		operator==
		(
		const Instrument& rhs
		);

private:
	unsigned _id;
	string _symbol;
	string _name;
};

// Ordinary Least Squares Regression Metrics
class OrdinaryLeastSquaresRegressionMetrics 
{
public:
	OrdinaryLeastSquaresRegressionMetrics
		(
		Instrument i_instrument,
		double i_std_dev,
		double i_std_err,
		double i_r_square,
		double i_r
		) : 
	_instrument (i_instrument),
		_std_dev    (i_std_dev),
		_std_err    (i_std_err),
		_r_square   (i_r_square),
		_r          (i_r)
	{}

	void pretty_print (void)
	{
		_instrument.pretty_print();
		cout << "OLS Metrics [Std. Dev: " << _std_dev << ", Std. Err: " << _std_err << ", R: " << _r << ", R square: " << _r_square << "]" << endl;
	}

	Instrument _instrument;
	double _std_dev;
	double _std_err;
	double _r_square;
	double _r;
};

// Portfolio Constituent
class PortfolioConstituent
{
public:
	PortfolioConstituent
		(
		Instrument &i_instrument,
		double     i_amount
		);

	PortfolioConstituent
		(
		void
		);

	PortfolioConstituent
		(
		const PortfolioConstituent &rhs
		);

	Instrument
		get_instrument 
		(
		void
		);

	double
		get_amount 
		(
		void
		);

	void
		set_amount 
		(
		const double i_amount
		);

private:
	Instrument _instrument;
	double _amount;
};

class Portfolio
{
public:
	Portfolio 
		(
		void
		);

	void
		insert
		(
		PortfolioConstituent &i_constituent
		);

	double
		get_instrument_weight
		(
		const Instrument &i_instrument
		);

	double
		get_instrument_amount
		(
		const Instrument &i_instrument
		);

	bool
		get_instrument_in_portfolio
		(
		const Instrument &i_instrument
		);

	bool
		iterate
		(
		PortfolioConstituent &o_pc
		);

private:
	Portfolio (const Portfolio&);

	std::vector<PortfolioConstituent> _portfolio;
};

typedef struct {
	double _average;
	double _max;
	double _lp;
	double _cdar95;

	void prettyPrint (void)
	{
		cout << "Draw Down Metrics [Average: "<< _average << ", Max: " << _max << ", LP: " << _lp << ", CDAR95: " << _cdar95 << "]" << endl;
	}
} DrawDownMetrics;

// Static Constants
typedef enum {
	CONTINUOUS = 0, // log
	PERIODIC   = 1, // geometric
	CompoundingMethod_max = 2
} CompoundingMethod;

typedef enum {
	_DATE = 0,
	_OPEN,
	_HIGH,
	_LOW,
	_CLOSE,
	_VOLUME,
	_ADJ_CLOSE
} YahooMarketDataColumn;

const unsigned NUM_TRADING_DAYS       = 250;
const unsigned NUM_WEEKS_IN_A_QUARTER = 13;
const unsigned NUM_WEEKS_IN_A_YEAR    = 52;

const char MARKET_DATA_MONTHLY = 'm';
const char MARKET_DATA_WEEKLY  = 'w';
const char MARKET_DATA_DAILY   = 'd';

static inline
	void
	validate_market_data_frequency
	(
	const char i_market_data_frequency
	)
{
	switch(i_market_data_frequency) {
	case MARKET_DATA_MONTHLY:
	case MARKET_DATA_WEEKLY:
	case MARKET_DATA_DAILY:
		return;

	default:
		break;
	}
	throw exception("Unhandled Market Data Frequency");
}

// Functions
// Inline Functions
inline void _LOG   ( const string &msg ) { cout << msg         << endl; }
inline void _WARN  ( const string &msg ) { cout << "WARNING: " << msg << endl; }
inline void _ABORT ( const string &msg ) { cerr << "ABORT: "   << msg << endl; exit(-1); } 

// CSV Manipulation
inline
	string
	get_csv_line_cell_as_string
	(
	string line,
	unsigned cell
	);

inline
	double
	get_csv_line_cell_as_double
	(
	string line,
	unsigned cell
	);

// String Manipulation
inline
	wstring
	s_to_ws
	(
	const string &s
	);

inline
	string
	construct_market_data_filename
	(
	const string ticker,
	const char frequency
	);

inline
	string
	construct_market_data_url
	(
	const string ticker,
	const char frequency = MARKET_DATA_MONTHLY
	);

string
	download_market_data
	(
	string instrument,
	const char frequency = MARKET_DATA_MONTHLY
	);

void
	populate_return_vector
	(
	string				i_market_data_filename,
	const unsigned		i_num_returns,
	std::vector<double>	&o_returns
	);

void
	populate_price_vector
	(
	string				i_market_data_filename,
	const unsigned		i_num_prices,
	std::vector<double> &o_prices
	);

void
	populate_return_vector
	(
	string				i_market_data_filename,
	const unsigned		i_num_returns,
	const unsigned      i_matrix_row,
	matrix<double>      &i_matrix
	);

inline
	string
	convert_int_to_string
	(
	int i_int
	)
{
	stringstream ss;
	ss << i_int;
	return ss.str();
}

template <class t>
inline t
	list_max
	(
	const std::vector<t> &i_list
	);

template <class t>
inline t
	list_min
	(
	const std::vector<t> &i_list
	);

template <class t>
inline t
	list_mean
	(
	const std::vector<t> &i_list
	);

template <class t>
inline void
	matrix_mean
	(
	const matrix<t> &i_matrix,
	std::vector<t>  &o_means
	);

template <class t>
bool
	invert_matrix
	(
	const matrix<t> &input,
	matrix<t>       &inverse
	);

void
	pretty_print
	(
	gsl_vector *i_vector,
	size_t     i_len
	);

void
	pretty_print
	(
	gsl_matrix *i_matrix,
	size_t      i_rows,
	size_t      i_cols
	);

template <class t>
void
	pretty_print
	(
	const matrix<t> &i_matrix
	);

double
	get_return
	(
	string		&i_instrument,
	const char	i_frequency,
	const int	i_num_periods,
	int         i_offset = 0
	);

void
	ret_to_price 
	(
	std::vector<double> &o_prices,
	std::vector<double>	&i_ret_series,
	const double		i_start_price = 1.f,
	CompoundingMethod	i_method	  = CompoundingMethod::PERIODIC
	);

double
	ret_to_price
	(
	double            i_return,
	double            i_start_price  = 1.f,
	CompoundingMethod i_method       = CompoundingMethod::PERIODIC
	);

void
	get_drawdown
	(
	std::vector<DrawDownMetrics>	&o_data,
	matrix<double>					&i_returns
	);

void
	get_drawdown
	(
	DrawDownMetrics		&o_data,
	std::vector<double>	i_returns
	);

void
	univariate_regression
	(
	std::vector<double>		i_deps,
	std::vector<Instrument>	i_instruments, 
	matrix<double>			i_indeps,
	const unsigned			i_starting_index,
	const unsigned			i_window_size, 
	std::vector<OrdinaryLeastSquaresRegressionMetrics> &o_result
	);

typedef struct {
	double _r_sqr;
	std::vector<double> _coefficients;
} multivariate_regression_results_t;

void
	multivariate_regression
	(
	gsl_vector * i_deps,
	size_t       i_num_data_points,
	gsl_matrix * i_indeps,
	size_t       i_num_factors,
	multivariate_regression_results_t &o_result
	);

void
	replicate_mutual_fund_with_etfs
	(
	std::vector<multivariate_regression_results_t> &o_results,
	string				i_mutual_fund,
	std::vector<string>	i_etfs,
	const unsigned		i_num_regressions = 1,
	const unsigned		i_regression_window = NUM_WEEKS_IN_A_YEAR,
	ostream				&os = std::cout
	);
#endif