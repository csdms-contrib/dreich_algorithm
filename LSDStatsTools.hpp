//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDStatsTools
// Land Surface Dynamics StatsTools
//
// A collection of statistical routines for use with the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2013 Simon M. Mudd 2013
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;


#ifndef StatsTools_H
#define StatsTools_H

// computes linear regression
// replaces data in residuals with residuals and returns a 4 element vector, which has slope, intercept, r^2 and
// the Durbin-Watson test statistic which looks for autocorrelation of the residuals
vector<float> simple_linear_regression(vector<float>& x_data, vector<float>& y_data, vector<float>& residuals);
float get_mean(vector<float>& y_data);
float get_mean_ignore_ndv(Array2D<float>& data, float ndv);
float get_SST(vector<float>& y_data, float mean);
float get_variance_ignore_ndv(Array2D<float>& data, float ndv, float mean);
float get_durbin_watson_statistic(vector<float> residuals);
float get_standard_deviation(vector<float>& y_data, float mean);
float get_standard_error(vector<float>& y_data, float standard_deviation);
vector<float> get_common_statistics(vector<float>& y_data);
float get_percentile(vector<float>& data, float percentile);

// these look for linear segments within a data series.
void populate_segment_matrix(int start_node, int end_node, float no_data_value,
								vector<float>& all_x_data, vector<float>& all_y_data, int maximum_segment_length,
								float sigma, Array2D<float>& like_array, Array2D<float>& m_array,
								Array2D<float>& b_array, Array2D<float>& rsquared_array,
								Array2D<float>& DW_array);
void calculate_segment_matrices(vector<float>& all_x_data, vector<float>& all_y_data, int maximum_segment_length,
								float sigma, Array2D<float>& like_array, Array2D<float>& m_array,
								Array2D<float>& b_array, Array2D<float>& rsquared_array,
								Array2D<float>& DW_array);
void find_max_like_of_segments(int minimum_segment_length, Array2D<float>& like_array,
								vector<float>& max_MLE, vector< vector<int> >& segments_for_each_n_segments);
void find_max_AIC_of_segments(int minimum_segment_length, vector<float>& all_x_data, vector<float>& all_y_data,
								Array2D<float>& like_array,
								vector<float>& max_MLE, vector<float>& AIC_of_segments,
								vector<float>& AICc_of_segments, vector< vector<int> >& segments_for_each_n_segments);
void calculate_AIC_of_segments_with_normalized_sigma(float sigma,
								vector<float>& one_sigma_max_MLE, vector<float>& all_x_data,
								vector<float>& AIC_of_segments,vector<float>& AICc_of_segments);

// the below function is the main driver for the segment fitting code.
void best_fit_driver_AIC_for_linear_segments(int minimum_segment_length, float sigma,
											vector<float> all_x_data, vector<float> all_y_data,
											vector<float>& max_MLE);

// this gets the number of segments for several different values of sigma bassed in the vector sigma_values
void get_n_segments_for_various_sigma(vector<float> sigma_values, vector<float> one_sig_max_MLE,
									  vector<float>& all_x_data,
								      vector<int>& best_fit_AIC, vector<int>& best_fit_AICc,
								      vector< vector<float> >& AIC_for_each_n_segments,
								      vector< vector<float> >& AICc_for_each_n_segments);

// this prints full AIC and AICc information to screen
void print_AIC_and_AICc_to_screen(vector<float> sigma_values, vector< vector<int> > segments_for_each_n_segments,
								      vector<int> best_fit_AIC, vector<int> best_fit_AICc,
								      vector< vector<float> > AIC_for_each_n_segments,
								      vector< vector<float> > AICc_for_each_n_segments);

// this prints information about the most likeley segments to screen
void print_to_screen_most_likeley_segment_lengths( vector< vector<int> > segments_for_each_n_segments,
										vector<float> MLE_for_segments);

// this returns the m, b, r2 and DW stats of each segment
void get_properties_of_best_fit_segments(int bestfit_segments_node, vector< vector<int> >& segments_for_each_n_segments,
										 vector<float>& m_values, Array2D<float>& m_array,
										 vector<float>& b_values, Array2D<float>& b_array,
										 vector<float>& r2_values, Array2D<float>& rsquared_array,
										 vector<float>& DW_values, Array2D<float>& DW_array);

// these functions manipulate likelihood matrices and vectors for use with the segment tool
Array2D<float> normalize_like_matrix_to_sigma_one(float sigma, Array2D<float>& like_array);
vector<float> normalize_like_vector_to_sigma_one(float sigma, vector<float> like_vector);
Array2D<float> change_normalized_like_matrix_to_new_sigma(float sigma, Array2D<float>& sig1_like_array);
vector<float> change_normalized_like_vector_to_new_sigma(float sigma, vector<float> sig1_like_vector);

// this uses a moving window to find segments and is incomplete
void find_linear_segments(vector<float>& all_x_data, vector<float>& all_y_data, int segment_length);


// functions for partitioning and permutation (to be used with linear segment finding
int partitions_min( int x, int y);
void partition_print(int t, vector<int>& p);
void partitions_with_minimum_length(int n, int k, int t, int min_length, vector<int>& p);
void partitions_with_minimum_length(int n, int k, int t, int min_length, vector<int>& p,
								vector< vector < vector<int> > >& partitions);
void integer_partition(int n, int k, int t, vector<int>& p);
void partition_driver_to_screen(int n, int minimum_length);
vector< vector < vector<int> > > partition_driver_to_vecvecvec(int k, int minimum_length);
void partition_assign(int t, vector<int>& p, vector< vector < vector<int> > >& partitions);
void partition_vecvecvec_print(vector< vector < vector<int> > >& partitions);
void partition_vecvecvec_print_with_permutation(vector< vector < vector<int> > >& partitions);
void permute_partitioned_integer_vector(vector<int> permute_vector);

// this generates random segments for use in testing the segment finding algorithm
void generate_random_segments(float sigma, int minimum_n_nodes, int mean_segment_length, int segment_range,
							  float dx, float offset_range, float m_range,
							 vector<float>& x_data, vector<float>& y_data,
							 vector<int>& segment_length, vector<float>& slope, vector<float>& intercept);


// maxiumum likihood estimators
float calculate_MLE(vector<float>& measured, vector<float>& modelled, vector<float>& sigma);
float calculate_MLE(vector<float>& measured, vector<float>& modelled, float sigma);
float calculate_MLE_from_residuals(vector<float>& residuals, float sigma);

// a random number generator
float ran3( long *idum );
// Randomly sample from a vector without replacement DTM 21/04/2014
vector<float> sample_without_replacement(vector<float> population_vector, int N);
vector<int> sample_without_replacement(vector<int> population_vector, int N);
// conversion from numbers to strings
string itoa(int num);
string dtoa(float num);


// Log binning module
// two overloaded functions:
//    -> for data stored in a 2D array (e.g. slope-area)
void log_bin_data(Array2D<float>& InputArrayX, Array2D<float>& InputArrayY, float log_bin_width, vector<float>&  MeanX_output, vector<float>& MeanY_output,
                      vector<float>& midpoints_output, vector<float>& StandardDeviationX_output, vector<float>& StandardDeviationY_output,
                      vector<float>& StandardErrorX_output, vector<float>& StandardErrorY_output, vector<int>& num_observations, float NoDataValue);

//    -> for data stored in a 1D vector (e.g. for spectral analysis)
void log_bin_data(vector<float>& InputVectorX, vector<float>& InputVectorY, float log_bin_width,
                  vector<float>&  MeanX_output, vector<float>& MeanY_output,
                      vector<float>& midpoints_output, vector<float>&  StandardDeviationX_output,
                      vector<float>&  StandardDeviationY_output, int NoDataValue);

// Regular binning algoritm for data stored in a 1D vector
void bin_data(vector<float>& InputVectorX, vector<float>& InputVectorY, float bin_width,
                  vector<float>&  MeanX_output, vector<float>& MeanY_output,
                      vector<float>& midpoints_output, vector<float>& MedianY_output,
                      vector<float>&  StandardDeviationX_output, vector<float>&  StandardDeviationY_output,
                      vector<float>& StandardErrorX_output, vector<float>& StandardErrorY_output, 
                      vector<int>& number_observations_output, float& bin_lower_limit, float NoDataValue);

//look for empty bins output from the log binning function and removes them to avoid 
//plotting several empty bins at 0,0 in some cases. SWDG 6/11/13
void RemoveSmallBins(vector<float>&  MeanX_output, vector<float>& MeanY_output,
                      vector<float>& midpoints_output, vector<float>& StandardDeviationX_output, vector<float>& StandardDeviationY_output,
                      vector<float>& StandardErrorX_output, vector<float>& StandardErrorY_output, vector<int>& number_observations, float bin_threshold);

// Load in a vector of data and convert into a histogram with a specified bin width
// that is printed to file containing:
//    Midpoint LowerLim UpperLim Count ProbabilityDensity                      
void print_histogram(vector<float> input_values, float bin_width, string filename);

// This is a much simpler version of the binning software.  It takes two vectors, and
// sorts the values held within the first vector into bins according to their respective
// values in the second vector.  The output is a vector<vector> with the binned dataset.
// and a vector of bin midpoints.  These can then be analysed ahd plotted as desired.
// DTM 14/04/2014
void bin_data(vector<float>& vector1, vector<float>& vector2, float min, float max, float bin_width, vector<float>& mid_points, vector< vector<float> >& binned_data);

// tools for sorting
template<class T> struct index_cmp;
void matlab_float_sort(vector<float>& unsorted, vector<float>& sorted, vector<size_t>& index_map);
void matlab_float_reorder(std::vector<float> & unordered, std::vector<size_t> const & index_map, std::vector<float> & ordered);
void matlab_float_sort_descending(vector<float>& unsorted, vector<float>& sorted, vector<size_t>& index_map);
void matlab_int_sort(vector<int>& unsorted, vector<int>& sorted, vector<size_t>& index_map); // added 27/11/13 SWDG
void matlab_int_reorder(std::vector<int> & unordered, std::vector<size_t> const & index_map, std::vector<int> & ordered);

// Generate spline curves from X and Y vectors of floats
Array2D<float> CalculateCubicSplines(vector<float> X, vector<float> Y);
void PlotCubicSplines(vector<float> X, vector<float> Y, int SplineResolution, vector<float>& Spline_X, vector<float>& Spline_Y);

//Get vector of unique values in an input array of ints
vector<int> Unique(Array2D<int> InputArray, int NoDataValue);

//Get vector of unique values in an input array of floats
vector<float> Unique(Array2D<float> InputArray, int NoDataValue);

// Generate vector of evenly spaced numbers between two points
vector<float> linspace(float min, float max, int n);

// convert degree bearing from north to radians from east
float BearingToRad(float Bearing);

// conversion from degrees to radians
float rad(float degree);

// conversion from radians to degrees
float deg(float radians);

//Method to generate Statistical distribution. - DTM
void get_distribution_stats(vector<float>& y_data, float& mean, float& median, float& UpperQuartile, float& LowerQuartile, float& MaxValue);

// Method to calculate the quadratic mean. - DTM
double get_QuadraticMean(vector<double> input_values, double bin_width);

#endif




