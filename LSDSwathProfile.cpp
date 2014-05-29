//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDSwathProfile.cpp
//------------------------------------------------------------------------------
// This code is used to create swath profiles from either raster or sparse data.
// These two options will be coded up into two different classes, in order to
// deal with different types of data from which to construct the profiles.
// 
// The generalised swath profile framework for constructing transverse profiles
// is derived from the algorithm described by Hergarten et al. 2013: Generalized
// swath pro?les. This will be extended to also produce generalised longitudinal 
// profiles, as well as some other functionality as desired.
//
// The input to the swath profile objects includes the profile itself, which is
// loaded as a PointData structure as described in LSDShapeTools.hpp.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// David T. Milodowski, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Version 0.0.1		17/02/2014
// Prerequisite software packages: TNT, PCL and liblas
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>

// LSDTopotools objects
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDCloudBase.hpp"
#include "LSDShapeTools.hpp"
#include "LSDSwathProfile.hpp"
#include "LSDStatsTools.hpp"

// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>
// liblas
#include <liblas/liblas.hpp>

// TNT
#include "TNT/tnt.h"

using namespace std;
using namespace TNT;

#ifndef LSDSwathProfile_CPP
#define LSDSwathProfile_CPP

//------------------------------------------------------------------------------
// VECTOR GEOMETRY FUNCTIONS
// Functions to solve some geometric problems, such as shortest distance from a
// point to a line etc.

// This function calculates the shortest distance from a point, V_p, to a
// straight line passing through points V_a and V_b, using the vector triple
// product.
float calculate_shortest_distance_to_line(vector<float> V_a, vector<float> V_b, vector<float> V_p)
{
  float d = sqrt(( (V_b[0]-V_a[0])*(V_a[1]-V_p[1]) - (V_b[1]-V_a[1])*(V_a[0]-V_p[0]) )*( (V_b[0]-V_a[0])*(V_a[1]-V_p[1]) - (V_b[1]-V_a[1])*(V_a[0]-V_p[0]) ))
                                                    / sqrt( ((V_b[0]-V_a[0])*(V_b[0]-V_a[0])) + ((V_b[1]-V_a[1])*(V_b[1]-V_a[1])) );
//  cout << "V_a " << V_a[0] << "," << V_a[1] << "\tV_b " << V_b[0] << "," << V_b[1] << "\tV_p " << V_p[0] << "," << V_p[1] << "\td " << d << endl;
  return d;
}
// For a vector line equation V_star = V_a + (V_b-V_a)*t, this function
// calculates t to find the intersection point between the line joining points 
// V_a and V_b and the shortest path from a point V_p and this line.
float calculate_t(vector<float> V_a, vector<float> V_b, vector<float> V_p)
{
  float t = -( (V_b[0]-V_a[0])*(V_a[0]-V_p[0]) + (V_b[1]-V_a[1])*(V_a[1]-V_p[1]) ) / ((V_b[0]-V_a[0])*(V_b[0]-V_a[0]) + (V_b[1]-V_a[1])*(V_b[1]-V_a[1]));
  return t;
}
// Calculate disctance between two points
float calculate_distance_between_two_points(vector<float> V_a, vector<float> V_b)
{
  float distance = sqrt((V_a[0]-V_b[0])*(V_a[0]-V_b[0])+(V_a[1]-V_b[1])*(V_a[1]-V_b[1]));
  return distance;
}

// Use cross product to test whether point lies on left or right hand side of
// a line vector in 2D.  Returns true if the point is on the left.
bool test_point_left(vector<float> V_a, vector<float> V_b, vector<float> V_p)
{
  bool test;
  float cross_product = (V_b[0]-V_a[0])*(V_p[1]-V_a[1]) - (V_b[1]-V_a[1])*(V_p[0]-V_a[0]);
  if(cross_product>0) test = true;
  else test = false; 
}
//------------------------------------------------------------------------------
// CREATE FUNCTIONS
// the create function. This is default and throws an error
void LSDSwath::create()
{
	cout << "LSDSwathProfile line 102 You need to supply a PointData container of profile coordinates, a raster template and the half width of the profile" << endl;
	exit(EXIT_FAILURE);
}
// This function creates the swath profile template for creating a swath profile
// of a raster dataset.
void LSDSwath::create(PointData& ProfilePoints, LSDRaster& RasterTemplate, float HalfWidth)
{
  NPtsInProfile = ProfilePoints.X.size();   // Number of points in profile
  NoDataValue = RasterTemplate.get_NoDataValue();
  NRows = RasterTemplate.get_NRows();
  NCols = RasterTemplate.get_NCols();
  ProfileHalfWidth = HalfWidth;
  float Resolution = RasterTemplate.get_DataResolution();
  // Loop through profile points, calculating cumulative distance along profile
  // with each iteration
  vector<float> DistanceAlongBaseline_temp(NPtsInProfile,NoDataValue);
  float cumulative_distance = 0;
  DistanceAlongBaseline_temp[0]=cumulative_distance;
  for(int i = 1; i<NPtsInProfile; ++i)
  {
    cumulative_distance += sqrt((ProfilePoints.X[i]-ProfilePoints.X[i-1])*(ProfilePoints.X[i]-ProfilePoints.X[i-1])
                                              +(ProfilePoints.Y[i]-ProfilePoints.Y[i-1])*(ProfilePoints.Y[i]-ProfilePoints.Y[i-1]));
    DistanceAlongBaseline_temp[i]=cumulative_distance;
  }
  DistanceAlongBaseline = DistanceAlongBaseline_temp;
  // Read profile data into a LSDCloud object for querying
  vector<float> zeros(NPtsInProfile,0.0);
  LSDCloud ProfileCloud(ProfilePoints, RasterTemplate);
  // For each point in array, find nearest point along the profile and calculate
  // signed distance to profile.  The convention here is that points lying on
  // left hand side of profile as you traverse from start to finish are 
  // considered positive, and those on the right are negative.
  Array2D<float> DistanceToBaseline_temp(NRows,NCols,NoDataValue);
  Array2D<float> ProjectedDistanceAlongBaseline_temp(NRows,NCols,NoDataValue);
  
  float searchPoint_x,searchPoint_y;
  vector<float> SquaredDistanceToBaseline,temp,temp2;
  vector<int> ProfilePointIndex;
  int K=2;

  // Define bounding box of swath profile
  YMin = ProfileCloud.get_YMin();
  YMax = ProfileCloud.get_YMax();
  XMin = ProfileCloud.get_XMin();
  XMax = ProfileCloud.get_XMax();
  int ColStart = int(floor(XMin/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;
  
  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));  
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    cout << flush << "\t\t\t row " << i+1 << "/" << NRows << "\r";
    for(int j = ColStart; j<ColEnd; ++j)
    {
      vector<float> V_p;  // the search point coordinates
      V_p.push_back(float(j));
      V_p.push_back(float(NRows - 1) - float(i));
      // Find nearest two points on profile baseline
      ProfileCloud.NearestNeighbourSearch2D(V_p[0], V_p[1], K, temp, ProfilePointIndex, SquaredDistanceToBaseline);

      if(SquaredDistanceToBaseline.size()<1) cout << "No profile points found - check input files";
      else if(SquaredDistanceToBaseline.size()<2) cout << "Only one point found - check input files";
      else 
      {
        // SCENARIO 1
        //----------------------------------------------------------------------
        // Check that two points are adjacenet on baseline - simplest (and
        // usual) scenario)
        vector<float> V_a, V_b; // first, second nearest points on profile
        int V_a_index, V_b_index;
        if(abs(ProfilePointIndex[0]-ProfilePointIndex[1]) <= 1)
        {
//          cout << "SCENARIO 1" << endl;
          // If yes, then find point along cojoining line that minimises this
          // distance.
          // Note that this point is located at the point given by the vector:
          // (V_a + (V_b - V_a)*t).
          // t is calculated using the equation:
          // t = -(V_b-V_a).(V_a-V_p)/(|V_b - V_a|^2)
          V_a_index = ProfilePointIndex[0];
          V_b_index = ProfilePointIndex[1];
          V_a.push_back(ProfileCloud.get_point_x(V_a_index));
          V_a.push_back(ProfileCloud.get_point_y(V_a_index));
          V_b.push_back(ProfileCloud.get_point_x(V_b_index));
          V_b.push_back(ProfileCloud.get_point_y(V_b_index));       
        }
        
        // SCENARIO 2
        //----------------------------------------------------------------------
        // Two nearest points are not adjacent on the profile.  This is quite
        // common on the inside of a bend on a baseline.  At this point, the
        // strategy used is to find all baseline points located at a radius that
        // is equal to that of the nearest baseline point.  Then, the strategy
        // is to search through all located baseline points, and find the
        // nearest point that neighbours one of these baseline points.  The
        // program then searches through to find the closest point on the 
        // baseline between those points.
        else
        {
          cout << "SCENARIO 2" << endl;
          vector<int> ClosestPointIndex;
          float shortest_distance = NoDataValue;
          float distance_to_point;
          ProfileCloud.RadiusSearch2D(V_p[0], V_p[1], sqrt(SquaredDistanceToBaseline[0]), temp, ClosestPointIndex, temp2);
          for(int i_point_iter = 0; i_point_iter < ClosestPointIndex.size(); ++i_point_iter)
          {
            if((ClosestPointIndex[i_point_iter]!=0) && (ClosestPointIndex[i_point_iter]!=NPtsInProfile-1))
            {
              // Check step backwards
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]-1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]-1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]-1;
              }
              // Check step forwards
              V_b_test.clear();
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]+1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]+1));              
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]+1;
              }
            }
            // conditions for edge of profile
            else if(ClosestPointIndex[i_point_iter]==0)
            {
              // Check step forwards only
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]+1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]+1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]+1;
              }
            }
            else if(ClosestPointIndex[i_point_iter]==NPtsInProfile-1)
            {
              // Check step backwards only
              vector<float> V_a_test, V_b_test;
              V_a_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]));
              V_a_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]));
              V_b_test.push_back(ProfileCloud.get_point_x(ClosestPointIndex[i_point_iter]-1));
              V_b_test.push_back(ProfileCloud.get_point_y(ClosestPointIndex[i_point_iter]-1));
              distance_to_point = calculate_distance_between_two_points(V_p,V_b_test);
              if(shortest_distance == NoDataValue || shortest_distance > distance_to_point)
              {
                shortest_distance = distance_to_point;
                V_a = V_a_test;
                V_b = V_b_test;
                V_a_index = ClosestPointIndex[i_point_iter];
                V_b_index = ClosestPointIndex[i_point_iter]-1;
              }
            }
          }
        }
        // calculate position along baseline vector 
        float t = calculate_t(V_a, V_b, V_p);
        float d;
        if(t>0 && t<=1)
        // find the distance to the nearest point along the straight line
        // segment that links the two points
        {
//          cout << "test1" << endl;
          d = calculate_shortest_distance_to_line(V_a,V_b,V_p);
          if(d<ProfileHalfWidth)
          {
            ProjectedDistanceAlongBaseline_temp[i][j] = DistanceAlongBaseline[V_a_index]+(DistanceAlongBaseline[V_b_index]-DistanceAlongBaseline[V_a_index])*t;
            // determine which side of the profile the point is on using the cross
            // product
            if(V_a_index < V_b_index)
            {
              if(test_point_left(V_b,V_a,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d; 
            }
            else
            {
              if(test_point_left(V_a,V_b,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d; 
            }
          }
          else
          {
            DistanceToBaseline_temp[i][j] = NoDataValue;
            ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
          }
        }
        
        else if(V_a_index==0 || V_a_index==NPtsInProfile-1)
        // Avoid end points
        {
//          cout << "test2" << endl;
          DistanceToBaseline_temp[i][j] = NoDataValue;
          ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
        }
        else
        // Select nearest point (e.g. on outer side of bend)
        {
//          cout << "test3" << endl;
          d = calculate_distance_between_two_points(V_a,V_p);
          if(d<ProfileHalfWidth)
          {
            ProjectedDistanceAlongBaseline_temp[i][j] = DistanceAlongBaseline[V_a_index];
            // determine which side of the profile the point is on using the cross
            // product.  In this case, take the baseline vector as being formed
            // by straight line that joins the points on either side.
            vector<float> V_c(2,0.0);
            if(V_a_index < V_b_index)
            {
              V_c[0]=ProfileCloud.get_point_x(V_a_index-1);
              V_c[1]=ProfileCloud.get_point_y(V_a_index-1);
              if(test_point_left(V_b,V_c,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;         
            }
            else
            {
              V_c[0]=ProfileCloud.get_point_x(V_a_index+1);
              V_c[1]=ProfileCloud.get_point_y(V_a_index+1);
              if(test_point_left(V_c,V_b,V_p)==false) DistanceToBaseline_temp[i][j] = -1*d;
              else DistanceToBaseline_temp[i][j] = d;            
            }
          }
          else
          {
            DistanceToBaseline_temp[i][j] = NoDataValue;
            ProjectedDistanceAlongBaseline_temp[i][j] = NoDataValue;
          }
        }
      }       
    }
  }
  DistanceToBaselineArray = DistanceToBaseline_temp.copy();
  DistanceAlongBaselineArray = ProjectedDistanceAlongBaseline_temp.copy();
}

//------------------------------------------------------------------------------
// SWATH PROFILE GENERATION
// These routines take a swath profile template, comprising the LSDSwath object,
// and then uses this to construct either transverse (normal to profile) or
// longitudinal (parallel to profile) profiles.

// GET_TRANSVERSE_SWATH_PROFILES
// Function takes a raster and calculates transverse swath profiles, based on
// the swath template in the LSDSwath object.  Note that the input raster at
// present must have the same extent as the original template raster used to
// create the LSDSwath object.
// The function returns a vector container of the desired profiles.  The first
// and second profiles in the container are ALWAYS the mean and standard 
// deviation respectively.  The following profiles contain the desired 
// percentile profiles indicated in the input vector "desired_percentiles".
void LSDSwath::get_transverse_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles)
{  
  vector<float> TransverseDistance, RasterValues;
  float Resolution = Raster.get_DataResolution();
  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;
  
  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));  
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    for(int j = ColStart; j<ColEnd; ++j)
    {
      if((DistanceToBaselineArray[i][j]!=NoDataValue) && (Raster.get_data_element(i,j)!=NoDataValue))
      {
        TransverseDistance.push_back(DistanceToBaselineArray[i][j]);
        RasterValues.push_back(Raster.get_data_element(i,j));
      }
    }
  }
  // Sort values if you need percentiles  
  int NumberOfPercentileProfiles = desired_percentiles.size();
  if (NumberOfPercentileProfiles > 0)
  {
    vector<long unsigned int> index_map;
    matlab_float_sort(RasterValues,RasterValues,index_map);
    matlab_float_reorder(TransverseDistance, index_map, TransverseDistance);
  }
  
  // Bin data
  cout << "BIN DATA" << endl;
  vector< vector<float> > binned_RasterValues;
  bin_data(RasterValues, TransverseDistance, -ProfileHalfWidth, ProfileHalfWidth, BinWidth, mid_points, binned_RasterValues);
  cout << "done" << endl;
  // Produce desired profiles from binned_RasterValues.
  int NBins = mid_points.size();
  vector<float> mean_profile_temp(NBins, NoDataValue);
  vector<float> sd_profile_temp(NBins, NoDataValue);
  vector< vector<float> > output_percentile_profiles_temp;
  cout << "CALCULATING MEAN AND SD" << endl;
  for(int i = 0; i < NBins; ++i)
  {
    mean_profile_temp[i]=get_mean(binned_RasterValues[i]);
    sd_profile_temp[i]=get_standard_deviation(binned_RasterValues[i], mean_profile_temp[i]);
  }
  cout << "done" << endl;
  cout << "CALCULATING PERCENTILES" << endl;
  for(int j = 0; j < NumberOfPercentileProfiles; ++j)
  {
    vector<float> percentile(NBins,NoDataValue);
    for(int i = 0; i < NBins; ++i)
    {
      percentile[i] = get_percentile(binned_RasterValues[i], desired_percentiles[j]);
    }
    output_percentile_profiles_temp.push_back(percentile);
  }
  cout << "done" << endl;
  // Load profiles into export vectors
  mean_profile = mean_profile_temp;
  sd_profile = sd_profile_temp;
  output_percentile_profiles = output_percentile_profiles_temp; 
}

// GET_LONGITUDINAL_SWATH_PROFILES
// Function takes a raster and calculates longitudinal swath profiles, based on
// the swath template in the LSDSwath object.  Note that the input raster at
// present must have the same extent as the original template raster used to
// create the LSDSwath object.
// The function returns a vector container of the desired profiles.  The first
// and second profiles in the container are ALWAYS the mean and standard 
// deviation respectively.  The following profiles contain the desired 
// percentile profiles indicated in the input vector "desired_percentiles".
void LSDSwath::get_longitudinal_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles)
{  
  vector<float> LongitudinalDistance, RasterValues;
  float Resolution = Raster.get_DataResolution();
  // Define bounding box of swath profile
  int ColStart = int(floor((XMin)/Resolution));
  int ColEnd = ColStart + int(ceil((XMax-XMin)/Resolution));
  ColStart = ColStart - int(ceil(ProfileHalfWidth/Resolution));
  ColEnd = ColEnd + int(ceil(ProfileHalfWidth/Resolution));
  if (ColStart < 0) ColStart = 0;
  if (ColEnd > NCols) ColEnd = NCols;
  
  int RowEnd = NRows - 1 - int(floor(YMin/Resolution));
  int RowStart = RowEnd - int(ceil((YMax-YMin)/Resolution));
  RowStart = RowStart - int(ceil(ProfileHalfWidth/Resolution));
  RowEnd = RowEnd + int(ceil(ProfileHalfWidth/Resolution));  
  if (RowEnd > NRows) RowEnd = NRows;
  if (RowStart < 0) RowStart = 0;

  for(int i = RowStart; i<RowEnd; ++i)
  {
    for(int j = ColStart; j<ColEnd; ++j)
    {
      if((DistanceToBaselineArray[i][j]!=NoDataValue) && (Raster.get_data_element(i,j)!=NoDataValue))
      {
        LongitudinalDistance.push_back(DistanceToBaselineArray[i][j]);
        RasterValues.push_back(Raster.get_data_element(i,j));
      }
    }
  }
  // Sort values if you need percentiles  
  int NumberOfPercentileProfiles = desired_percentiles.size();
  if (NumberOfPercentileProfiles > 0)
  {
    vector<long unsigned int> index_map;
    matlab_float_sort(RasterValues,RasterValues,index_map);
    matlab_float_reorder(LongitudinalDistance, index_map, LongitudinalDistance);
  }
  
  // Bin data
  vector< vector<float> > binned_RasterValues;
  float start_point = DistanceAlongBaseline[0];
  float end_point = DistanceAlongBaseline[NPtsInProfile-1];
  bin_data(RasterValues, LongitudinalDistance, start_point, end_point, BinWidth, mid_points, binned_RasterValues);
  // Produce desired profiles from binned_RasterValues.
  int NBins = mid_points.size();
  vector<float> mean_profile_temp(NBins, NoDataValue);
  vector<float> sd_profile_temp(NBins, NoDataValue);
  vector< vector<float> > output_percentile_profiles_temp;
  for(int i = 0; i < NBins; ++i)
  {
    mean_profile_temp[i]=get_mean(binned_RasterValues[i]);
    sd_profile_temp[i]=get_standard_deviation(binned_RasterValues[i], mean_profile[i]);
  }
  for(int j = 0; j < NumberOfPercentileProfiles; ++j)
  {
    vector<float> percentile(NBins,NoDataValue);
    for(int i = 0; i < NBins; ++i)
    {
      percentile[i] = get_percentile(binned_RasterValues[i], desired_percentiles[j]);
    }
    output_percentile_profiles_temp.push_back(percentile);
  }
  // Load profiles into export vectors
  mean_profile = mean_profile_temp;
  sd_profile = sd_profile_temp;
  output_percentile_profiles = output_percentile_profiles_temp;
}

//------------------------------------------------------------------------------
// WRITE PROFILES TO FILE
// These routines take a swath profile template, comprising the LSDSwath object,
// and then uses this to construct either transverse (normal to profile) or
// longitudinal (parallel to profile) profiles.
void LSDSwath::write_transverse_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth, string prefix)
{
  string profile_extension = "_trans_profile.txt";
  vector<float> mid_points;
  vector<float> mean_profile;
  vector<float> sd_profile;
  vector< vector<float> > output_percentile_profiles;
  get_transverse_swath_profile(Raster, desired_percentiles, BinWidth, mid_points, mean_profile, sd_profile, output_percentile_profiles);
  // Print profiles to file
  string filename = prefix + profile_extension;
  cout << "\t printing profiles to " << filename << endl;
  ofstream ofs;
  ofs.open(filename.c_str());
  
  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
		exit(EXIT_FAILURE);
  }
  
  ofs << "Midpoint Mean SD ";
  for(int i_perc = 0; i_perc<desired_percentiles.size(); ++i_perc) ofs << desired_percentiles[i_perc] << " ";
  ofs << "\n";
  for(int i = 0; i < mid_points.size(); ++i)
  {
    ofs << mid_points[i] << " " << mean_profile[i] << " " << sd_profile[i] << " " ;
    for(int i_perc = 0; i_perc<desired_percentiles.size(); ++i_perc) ofs << output_percentile_profiles[i_perc][i] << " ";
    ofs << "\n";
  }
  
  ofs.close();
}
void LSDSwath::write_longitudinal_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth, string prefix)
{
  string profile_extension = "_long_profile.txt";
  vector<float> mid_points;
  vector<float> mean_profile;
  vector<float> sd_profile;
  vector< vector<float> > output_percentile_profiles;
  get_longitudinal_swath_profile(Raster, desired_percentiles, BinWidth, mid_points, mean_profile, sd_profile, output_percentile_profiles);
  // Print profiles to file
  string filename = prefix + profile_extension;
  cout << "\t printing profiles to " << filename << endl;
  ofstream ofs;
  ofs.open(filename.c_str());
  
  if(ofs.fail())
  {
    cout << "\nFATAL ERROR: unable to write output_file" << endl;
		exit(EXIT_FAILURE);
  }
  
  ofs << "Midpoint Mean SD ";
  for(int i_perc = 0; i_perc<desired_percentiles.size(); ++i_perc) ofs << desired_percentiles[i_perc] << " ";
  ofs << "\n";
  for(int i = 0; i < mid_points.size(); ++i)
  {
    ofs << mid_points[i] << " " << mean_profile[i] << " " << sd_profile[i] << " " ;
    for(int i_perc = 0; i_perc<desired_percentiles.size(); ++i_perc) ofs << output_percentile_profiles[i_perc][i] << " ";
    ofs << "\n";
  }
  
  ofs.close();
}
#endif