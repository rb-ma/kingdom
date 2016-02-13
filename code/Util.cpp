#include "Util.h"

Date::Date
	(
	int year,
	int month,
	int day
	) :
_year  (year),
	_month (month),
	_day   (day)
{
}

int
	Date::get_year
	(
	void
	)
{
	return _year;
}

int
	Date::get_month
	(
	void
	)
{
	return _month;
}

int
	Date::get_day
	(
	void
	)
{
	return _day;
}

Instrument::Instrument
	(
	unsigned id,
	string symbol,
	string name
	) :
_id (id),
	_symbol (symbol),
	_name (name)
{
}

Instrument::Instrument
	(
	void
	)
{}

Instrument::Instrument
	(
	const Instrument &rhs
	)
{
	_id = rhs._id;
	_symbol = rhs._symbol;
	_name = rhs._name;
}

unsigned
	Instrument::get_id
	(
	void
	)
{
	return _id;
}

string
	Instrument::get_symbol
	(
	void
	)
{
	return _symbol;
}

string
	Instrument::get_name
	(
	void
	)
{
	return _name;
}

void
	Instrument::pretty_print
	(
	void
	)
{
	cout << _symbol << " [id=" << _id << ",name=" << _name << "]" << endl;
}

bool
	Instrument::operator==
	(
	const Instrument &rhs
	)
{
	return 
		(_id     == rhs._id)     &&
		(_symbol == rhs._symbol) &&
		(_name   == rhs._name);
}

PortfolioConstituent::PortfolioConstituent
	(
	Instrument &i_instrument,
	double i_amount
	) :
_instrument (i_instrument),
	_amount (i_amount)
{}

PortfolioConstituent::PortfolioConstituent
	()
{}

PortfolioConstituent::PortfolioConstituent
	(
	const PortfolioConstituent &rhs
	)
{
	_instrument = rhs._instrument;
	_amount = rhs._amount;
}

Instrument
	PortfolioConstituent::get_instrument 
	(
	void
	) 
{ 
	return _instrument; 
}

double
	PortfolioConstituent::get_amount 
	(
	void
	) 
{ 
	return _amount; 
}

void
	PortfolioConstituent::set_amount 
	(
	const double i_amount
	) 
{ 
	_amount = i_amount; 
}

Portfolio::Portfolio 
	(
	void
	) {}

void
	Portfolio::insert
	(
	PortfolioConstituent &i_constituent
	)
{
	if (i_constituent.get_amount() <= 0.) return;

	for (auto iter = _portfolio.begin(); iter != _portfolio.end(); ++iter) {
		if (iter->get_instrument() == i_constituent.get_instrument()) { 
			iter->set_amount(i_constituent.get_amount());
			return;
		}
	}

	_portfolio.push_back(i_constituent);
}

double
	Portfolio::get_instrument_amount
	(
	const Instrument &i_instrument
	)
{
	for (auto iter = _portfolio.begin(); iter != _portfolio.end(); ++iter)
		if (iter->get_instrument() == i_instrument) 
			return iter->get_amount();

	return 0.f;
}

double
	Portfolio::get_instrument_weight
	(
	const Instrument &i_instrument
	)
{
	double amount = 0.f;
	double total = 0.f;
	for (auto iter = _portfolio.begin(); iter != _portfolio.end(); ++iter) {
		total += iter->get_amount();
		if (iter->get_instrument() == i_instrument) 
			amount = iter->get_amount();
	}

	if (total == 0.)
		throw exception("No Portfolio Constituents");

	return amount/total;
}

bool
	Portfolio::get_instrument_in_portfolio
	(
	const Instrument &i_instrument
	)
{
	return get_instrument_amount(i_instrument) != 0.f;
}

bool
	Portfolio::iterate
	(
	PortfolioConstituent &o_pc
	)
{
	static auto iter = _portfolio.begin();
	if (iter != _portfolio.end()) {
		iter = _portfolio.begin();
		return false;
	}
	o_pc = *iter;
	return true;
}

string
	get_csv_line_cell_as_string
	(
	string line,
	unsigned cell
	)
{
	stringstream line_stream(line);
	string line_cell;
	unsigned cell_i = 0;
	while (std::getline(line_stream, line_cell, ',')) {
		if (cell_i == cell)
			return line_cell;
		++cell_i;
	}

	throw exception("Unreachable Code");
	return "";
}

double
	get_csv_line_cell_as_double
	(
	string line,
	unsigned cell
	)
{
	return stod(get_csv_line_cell_as_string(line, cell));
}

wstring
	s_to_ws
	(
	const string &s
	)
{
	int slen = static_cast<int>(s.length() + 1);
	int len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slen, 0, 0);
	wchar_t * buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slen, buf, len);
	wstring r(buf);
	delete[] buf;
	return r;
}

template <class t>
t
	list_max
	(
	const std::vector<t> &i_list
	)
{
	if (i_list.empty()) 
		_ABORT("List Is Empty");

	t max = i_list[0];
	for (auto iter = i_list.begin()+1; iter != i_list.end(); ++iter)
		if (*iter > max) max = *iter;

	return max;
}

template <class t>
t
	list_min
	(
	const std::vector<t> &i_list
	)
{
	if (i_list.empty())
		_ABORT("List Is Empty");

	t min = i_list[0];
	for (auto iter = i_list.begin()+1; iter != i_list.end(); ++iter)
		if (*iter < min) min = *iter;

	return min;
}

template <class t>
t
	list_mean
	(
	const std::vector<t> &i_list
	)
{
	if (i_list.empty())
		_ABORT("List Is Empty");

	double i_list_size = static_cast<decltype(i_list_size)>(i_list.size());
	double acc = 0.f;
	for (auto iter = i_list.begin(); iter != i_list.end(); ++iter)
		acc += static_cast<decltype(acc)>(*iter);

	return static_cast<t>(acc/i_list_size);
}

template <class t>
t
	prctile
	(
	std::vector<t>	&i_list,
	double			i_percentile
	)
{
	if (i_list.empty())
		_ABORT("List Is Empty");

	sort(i_list.begin(), i_list.end());
	unsigned index = static_cast<decltype(index)>(ceil(i_percentile * static_cast<double>(i_list.size())));

	if (index >= i_list.size())
		index = static_cast<decltype(index)>(i_list.size() - 1u);

	return i_list[index];
}

string
	construct_market_data_filename
	(
	const string instrument	,
	const char frequency
	)
{
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	int year  = timeinfo->tm_year + 1900;
	int month = timeinfo->tm_mon  + 1;
	int day   = timeinfo->tm_mday;
	
	if (frequency == MARKET_DATA_WEEKLY)
		day -= (timeinfo->tm_wday - 1);
	if (frequency == MARKET_DATA_MONTHLY)
		day = 1;

	string file_prefix = 
		"C:\\Users\\Rohit Banerjee\\Documents\\md\\" + 
		instrument + 
		"_"        + 
		convert_int_to_string(year)  +  
		"_"        + 
		convert_int_to_string(month) + 
		"_";
	
	switch (frequency) {
	case MARKET_DATA_DAILY:
		return file_prefix + convert_int_to_string(day) + "_DAILY.csv"; 
		break;

	case MARKET_DATA_WEEKLY:
		return file_prefix + convert_int_to_string(day) + "_WEEKLY.csv";
		break;

	case MARKET_DATA_MONTHLY:
		return file_prefix + "MONTHLY.csv";
		break;
	}

	throw exception("Unreachable code");
	return "";
}

string
	convert_market_data_url
	(
	const string ticker,
	const char frequency
	)
{
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	int start_year = 2010;
	int year  = timeinfo->tm_year + 1900;
	int month = timeinfo->tm_mon  + 1;
	int day   = timeinfo->tm_mday;
	
	if (frequency == MARKET_DATA_WEEKLY)
		day -= (timeinfo->tm_wday - 1);
	if (frequency == MARKET_DATA_MONTHLY)
		day = 0;

	timeinfo->tm_mday;

	return 
		"http://real-chart.finance.yahoo.com/table.csv?s="
		+ ticker
		+ "&a="
		+ convert_int_to_string(0)
		+ "&b="
		+ convert_int_to_string(1)
		+ "&c="
		+ convert_int_to_string(start_year)
		+ "&d=" 
		+ convert_int_to_string(month - 1)
		+ "&e="
		+ convert_int_to_string(day)
		+ "&f="
		+ convert_int_to_string(year)
		+ "&g="
		+ frequency 
		+ "&ignore=.csv";
}

string
	download_market_data
	(
	string instrument,
	const char frequency
	)
{
	string file_str = construct_market_data_filename(instrument, frequency);
	ifstream ifile(file_str.c_str());
	if (ifile) return file_str;

	wstring url_str  = s_to_ws(convert_market_data_url(instrument, frequency));

	if (S_OK != URLDownloadToFile(0, url_str.c_str(), s_to_ws(file_str).c_str(), 0, NULL)) {
		// throw exception("Unable to download market data file");
		cout << "Could not download market data for " << instrument << " (" << convert_market_data_url(instrument, frequency) << ")" << endl;
		return "";
	}

	return file_str;
}

double
	get_return
	(
	string		&i_instrument,
	const char	i_frequency,
	const int	i_num_periods,
	int			i_offset
	)
{
	download_market_data(i_instrument, i_frequency);

	ifstream data(construct_market_data_filename(i_instrument, i_frequency).c_str());
	if(data.is_open()) {
		string line;
		std::getline(data, line);
		if (data.bad() || data.eof())
			throw new exception("Error reading Market Data file");

		while (data.good() && i_offset--)
			std::getline(data,line);

		double start, end;
		for(int period = 0; period <= i_num_periods; ++period) {
			std::getline(data, line);
			if (data.bad() || data.eof())
				throw new exception("Error reading Market Data file");

			if (period == 0)
				end = get_csv_line_cell_as_double(line, YahooMarketDataColumn::_CLOSE);
			if (period == i_num_periods) 
				start = get_csv_line_cell_as_double(line, YahooMarketDataColumn::_CLOSE);
		}
		return ((end/start) - 1.f);
	}
	throw new exception("Unable to parse Market Data CSV File");
	return -1.f;
}

void
	ret_to_price 
	(
	std::vector<double>	&o_prices,
	std::vector<double>	&i_ret_series,
	const double		i_start_price,
	CompoundingMethod	i_method
	)
{
	o_prices.clear();

	const size_t NUMOBS = i_ret_series.size();

	unsigned obs = 0;
	o_prices.push_back(i_start_price);

	for (obs = 1; obs < NUMOBS; ++obs) {
		double start = o_prices[obs-1];
		double end;

		try {
			end = ret_to_price(
				i_ret_series[obs-1],
				start,
				i_method
				);
		} catch (exception e) {
			_ABORT("Exceptional condition encountered");
		}

		o_prices.push_back(end);
	}
}

double
	ret_to_price
	(
	double            i_return,
	double            i_start_price,
	CompoundingMethod i_method
	)
{
	switch (i_method) {
	case PERIODIC:
		return i_start_price * (1.f + i_return);
		break;

	default:
		throw exception("Unhandled Compounding Method");
	}
}

void
	get_drawdown
	(
	std::vector<DrawDownMetrics>	&o_data,
	matrix<double>					&i_returns
	)
{
	try {
		for (unsigned sim = 0; sim < i_returns.size1(); ++sim) {
			matrix_row<matrix<double> > mr_returns(i_returns, sim);

			DrawDownMetrics instrument_dd_metrics;
			std::vector<double> returns(mr_returns.begin(), mr_returns.end());
			get_drawdown(
				instrument_dd_metrics,
				returns
				);

			o_data.push_back(instrument_dd_metrics);
		}			

	} catch (exception e) {
		_ABORT("Received Exceptional Condition ... Aborting");
	}
}

void
	get_drawdown
	(
	DrawDownMetrics		&o_data,
	std::vector<double>	i_returns	
	)
{
	const size_t NUMOBS    = i_returns.size();
	const size_t NUMPRICES = NUMOBS + 1;

	std::vector<double> asset_prices;
	ret_to_price(
		asset_prices,
		i_returns
		);

	// compute vector of running maximum upto time t
	std::vector<double> running_maximum;
	running_maximum.push_back(asset_prices[0]);
	for (unsigned i = 1; i < asset_prices.size(); ++i)
		running_maximum.push_back(running_maximum[i-1] + asset_prices[i]);

	// compute drawdowns for each time t
	// it is defined as the drop of the asset price from its running maximum
	std::vector<double> draw_down;
	for (unsigned i = 0; i < asset_prices.size(); ++i)
		draw_down.push_back(running_maximum[i] - asset_prices[i]);

	// Capture percentage draw down from running maximum
	std::vector<double> draw_down_pct;
	for (unsigned i = 0; i < asset_prices.size(); ++i)
		draw_down_pct.push_back(min(0.f, -draw_down[i]) / running_maximum[i]);

	// Compute statistics
	o_data._average = list_mean(draw_down_pct);
	o_data._max     = list_max(draw_down_pct);
	o_data._lp      = prctile(draw_down_pct, 5.f);

	std::vector<double> worst_dd;
	for (auto iter = draw_down_pct.begin(); iter != draw_down_pct.end(); ++iter)
		if (*iter < o_data._lp) worst_dd.push_back(*iter);

	o_data._cdar95 = list_mean(worst_dd);
}

template <class t>
bool 
	invert_matrix
	(
	const matrix<t>	&input, 
	matrix<t>			&inverse
	)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<t> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<t> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

template bool invert_matrix(const matrix<double> &input, matrix<double> &inverse);

void
	pretty_print
	(
	gsl_vector * i_vector,
	size_t     i_len
	)
{
	cout << "[";
	for (unsigned i = 0; i < i_len; ++i)
		cout << gsl_vector_get(i_vector, i) << ",";
	cout << "]" << endl;
}

void
	pretty_print
	(
	gsl_matrix * i_matrix,
	size_t       i_rows,
	size_t       i_cols
	)
{
	gsl_vector * row = gsl_vector_alloc(i_cols);
	for (unsigned i = 0; i < i_rows; ++i) {
		gsl_matrix_get_row(row, i_matrix, i);
		pretty_print(row, i_cols);
	}
}

template <class t>
void
	pretty_print
	(
	const matrix<t> &i_matrix
	)
{
	for (unsigned row = 0; row < i_matrix.size1(); ++row) {
		for (unsigned col = 0; col < i_matrix.size2(); ++col) {
			std::cout << i_matrix(row, col) << " ";
		}
		std::cout << std::endl;
	}
}

template void pretty_print(const matrix<double> &i_matrix);
void
	univariate_regression
	(
	std::vector<double>		i_deps,
	std::vector<Instrument>	i_instruments, 
	matrix<double>			i_indeps,
	const unsigned			i_starting_index,
	const unsigned			i_window_size, 
	std::vector<OrdinaryLeastSquaresRegressionMetrics> &o_result
	)
{
	_ASSERT(i_deps.size() >= (i_starting_index + i_window_size));
	_ASSERT(i_indeps.size2() == i_deps.size());
	_ASSERT(i_indeps.size1() == i_instruments.size());

	size_t num_indeps = i_instruments.size();

	// sum of x values
	std::vector<double> sum_x(num_indeps, 0.f);

	// sum of y values
	double sum_y = 0.f;

	// sum of x * y
	std::vector<double> sum_xy(num_indeps, 0.f);

	// sum of x^2
	std::vector<double> sum_xx(num_indeps, 0.f);

	// sum of squared residue
	double sum_res = 0.f;

	// slope of regression line
	std::vector<double> slope(num_indeps, 0.f);

	// y intercept of regression line
	std::vector<double> y_int(num_indeps, 0.f);

	// sum of squared of discrepancies
	std::vector<double> sum_yres(num_indeps, 0.f);

	// mean of y
	double avg_y = 0.f;

	// mean of x
	std::vector<double> avg_x(num_indeps, 0.f);

	// calculate various sums
	for (unsigned y_iter = i_starting_index; y_iter < i_starting_index + i_window_size; ++y_iter) {
		double y = i_deps[y_iter];
		sum_y += y;
		for (unsigned x_iter = 0; x_iter < num_indeps; ++x_iter) {
			double x = i_indeps(x_iter, y_iter);
			sum_x[x_iter]  += x;
			sum_xy[x_iter] += x*y;
			sum_xx[x_iter] += x*x;
		}
	}

	// calculate the means of x and y
	double data_size = static_cast<double>(i_window_size);
	for (unsigned x_iter = 0; x_iter < num_indeps; ++x_iter) {
		avg_x[x_iter] = sum_x[x_iter] / data_size;
	}
	avg_y = sum_y / data_size;

	// slope or a1
	for (unsigned x_iter = 0; x_iter < num_indeps; ++x_iter) {
		slope[x_iter] = 
			( (data_size * sum_xy[x_iter]) - (sum_x[x_iter] * sum_y) ) /
			( (data_size * sum_xx[x_iter]) - (sum_x[x_iter] * sum_x[x_iter]) );
	}

	// y intercept of a0
	for (unsigned x_iter = 0; x_iter < num_indeps; ++x_iter) {
		y_int[x_iter] = avg_y - slope[x_iter] * avg_x[x_iter];
	}

	// calculate squared residues, their sum etc.
	for (unsigned y_iter = i_starting_index; y_iter < (i_starting_index + i_window_size); ++y_iter) {
		double y = i_deps[y_iter];
		sum_res += pow(y - avg_y, 2.f);

		for (unsigned x_iter = 0; x_iter < num_indeps; ++x_iter) {
			double x = i_indeps(x_iter, y_iter);
			sum_yres[x_iter] += pow(y - y_int[x_iter] - (slope[x_iter] * x), 2.f);			
		}
	}

	for (unsigned x_iter = 0; x_iter < num_indeps; ++x_iter) {
		double r_sqr = (sum_res - sum_yres[x_iter]) / sum_res;
		o_result.push_back(OrdinaryLeastSquaresRegressionMetrics(
			i_instruments[x_iter], 
			sqrt(sum_res / (data_size - 1.f)), 
			sqrt(sum_yres[x_iter] / (data_size - 2.f)), 
			r_sqr, 
			sqrt(r_sqr)));
		o_result[x_iter].pretty_print();
	}
}

void
	multivariate_regression
	(
	gsl_vector * i_deps,
	size_t       i_num_data_points,
	gsl_matrix * i_indeps,
	size_t       i_num_factors,
	multivariate_regression_results_t &o_result
	)
{
	double y_sum = 0.;
	for (unsigned i = 0; i < i_num_data_points; ++i)
		y_sum += gsl_vector_get(i_deps, i);

	double y_avg = y_sum / static_cast<double>(i_num_data_points);
	double sum_res = 0.;
	for (unsigned i = 0; i < i_num_data_points; ++i)
		sum_res += pow(gsl_vector_get(i_deps, i) - y_avg, 2.f);

	gsl_vector * c   = gsl_vector_alloc(i_num_factors);
	gsl_matrix * cov = gsl_matrix_alloc(i_num_factors, i_num_factors);

	double chisq;
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(i_num_data_points, i_num_factors);
	gsl_multifit_linear(
		i_indeps,
		i_deps,
		c,
		cov,
		&chisq,
		work);

	double sum_yres = 0;
	for (unsigned i = 0; i < i_num_data_points; ++i) {
		double y_calc = 0.;
		for (unsigned j = 0; j < i_num_factors; ++j) {
			y_calc += gsl_vector_get(c, j) * gsl_matrix_get(i_indeps, i, j);
		}
		sum_yres += pow(gsl_vector_get(i_deps, i) - y_calc, 2.);
	}

	o_result._r_sqr = (sum_res - sum_yres) / sum_res;
	for (unsigned j = 0; j < i_num_factors; ++j)
		o_result._coefficients.push_back(gsl_vector_get(c, j));

	gsl_vector_free(c);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(work);
}

void
	populate_return_vector
	(
	string			i_market_data_filename,
	const unsigned	i_num_returns,
	std::vector<double>	&o_returns
	)
{
	_ASSERT(i_num_returns > 0);
	o_returns.clear();
	ifstream file(i_market_data_filename.c_str());
	if (file.is_open()) {
		string line;
		std::getline(file, line);	
		std::getline(file, line);
		double end = get_csv_line_cell_as_double(line, _ADJ_CLOSE);
		for (unsigned ret = 0; ret < i_num_returns && file.good(); ++ret) {
			std::getline(file, line);
			double start = get_csv_line_cell_as_double(line, _ADJ_CLOSE);
			o_returns.push_back((end/start) - 1.f);
			end = start;
		}
	}
}

void
	populate_price_vector
	(
	string				i_market_data_filename,
	const unsigned		i_num_prices,
	std::vector<double> &o_prices
	)
{
	_ASSERT(i_num_prices > 0);
	o_prices.clear();
	ifstream file(i_market_data_filename.c_str());
	if (file.is_open()) {
		string line;
		std::getline(file, line);	
		for (unsigned i = 0; i < i_num_prices && file.good(); ++i) {
			std::getline(file, line);
			if (!line.empty())
				o_prices.push_back(get_csv_line_cell_as_double(line, _ADJ_CLOSE));
		}
	}
}

void
	populate_return_vector
	(
	string i_market_data_filename,
	const unsigned i_num_returns,
	const unsigned i_matrix_row,
	matrix<double> &i_matrix
	)
{
	_ASSERT(i_num_returns > 0);

	ifstream file(i_market_data_filename.c_str());
	if (file.is_open()) {
		string line;
		// get the first line
		std::getline(file, line);	
		std::getline(file, line);
		double end = get_csv_line_cell_as_double(line, _ADJ_CLOSE);
		for (unsigned ret = 0; ret < i_num_returns && file.good(); ++ret) {
			std::getline(file, line);
			double start = get_csv_line_cell_as_double(line, _ADJ_CLOSE);
			i_matrix(i_matrix_row, ret) = ((end/start) - 1.f);
			end = start;
		}
	}
}

void
	replicate_mutual_fund_with_etfs
	(
	std::vector<multivariate_regression_results_t> &o_results,
	string				i_mutual_fund,
	std::vector<string> i_etfs,
	const unsigned		i_num_regressions,
	const unsigned		i_regression_window,
	ostream				&i_os
	)
{
	_ASSERT(i_regression_window != 0);
	i_os << "Regressing " << i_mutual_fund << " against [ ";
	for (auto iter = i_etfs.begin(); iter != i_etfs.end(); ++iter)
		i_os << (*iter) << " ";
	i_os << "]" << endl;

	// Download Mutual Fund Market Data
	const size_t num_data_points = i_regression_window + (i_num_regressions - 1);
	try {
		string mutual_fund_filename = download_market_data(
			i_mutual_fund, 
			MARKET_DATA_WEEKLY);

		std::vector<double> mutual_fund_returns;
		populate_return_vector(
			mutual_fund_filename, 
			num_data_points, 
			mutual_fund_returns
			);

		// Download ETF Market Data
		std::vector<Instrument> etfs;
		matrix<double> etf_returns(i_etfs.size(), num_data_points);

		for (unsigned i = 0; i < i_etfs.size(); ++i) {
			string etf_symbol = i_etfs[i];
			etfs.push_back(Instrument(0, etf_symbol, ""));

			string etf_md_filename = download_market_data(
				etf_symbol,
				MARKET_DATA_WEEKLY
				);

			populate_return_vector(
				etf_md_filename, 
				num_data_points,
				i,
				etf_returns
				);
		}

		// Rolling OLS
		for (unsigned regression = 0; regression < i_num_regressions; ++regression) {
			gsl_vector * this_reg_mf_returns = gsl_vector_alloc(static_cast<const size_t>(i_regression_window));
			for (unsigned i = regression; i < (regression + i_regression_window); ++i)
				gsl_vector_set(this_reg_mf_returns, i-regression, mutual_fund_returns[i]);

			gsl_matrix * this_reg_etf_returns = gsl_matrix_alloc(i_regression_window, i_etfs.size()+1);
			for (unsigned i = 0; i < i_regression_window; ++i) {
				gsl_matrix_set(this_reg_etf_returns, i, 0, 1.);
				for (unsigned j = 0; j < i_etfs.size(); ++j) {
					gsl_matrix_set(this_reg_etf_returns, i, j+1, etf_returns(j, i+regression));
				}
			}

			multivariate_regression_results_t reg_results;
			multivariate_regression(
				this_reg_mf_returns,
				i_regression_window,
				this_reg_etf_returns,
				i_etfs.size()+1,
				reg_results);

			i_os << "Regression " << regression << endl;
			i_os << "R-squared: " << reg_results._r_sqr << endl;
			i_os << "Constant Term: " << reg_results._coefficients[0] << endl;
			for (unsigned i = 0; i < i_etfs.size(); ++i)
				i_os << etfs[i].get_symbol() << ": " << reg_results._coefficients[i+1] << endl;
			i_os << endl;

			o_results.push_back(reg_results);

			gsl_vector_free(this_reg_mf_returns);
			gsl_matrix_free(this_reg_etf_returns);
		}
	} catch (exception e) {
		cerr << "Caught exception " << endl;
		return;
	}
}
