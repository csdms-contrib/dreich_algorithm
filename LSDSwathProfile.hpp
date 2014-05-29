//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDSwathProfile.hpp
//------------------------------------------------------------------------------
// This code houses the LSDCloud object, and associated functions, designed to
// analyse 3D pointcloud data, such as airborne LiDAR, and interface with the
// raster based LSDTopotools.  Currently reads .las files, but this could be
// expanded in the future to include other input types.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// David T. Milodowski, University of Edinburgh
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Version 0.0.1		17/02/2014
// Prerequisite software packages: TNT, PCL and liblas
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDCloudBase.hpp"
#include "LSDShapeTools.hpp"
// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>
// liblas
#include <liblas/liblas.hpp>
using namespace std;
using namespace TNT;

#ifndef LSDSwathProfile_H
#define LSDSwathProfile_H

class LSDSwath
{
  public:   
  LSDSwath()	{ create(); }
  ///@brief create an LSDSwath using a raster as a template.
  ///
  ///@param PointData ProfilePoints -> coordinates of points forming the
  /// baseline of the profile.
  ///@param LSDRaster RasterTemplate -> a raster dataset that is used as a
  /// template for the swath profile i.e. any rasters that you wish to generate
  /// the profile for should have the same characteristics/spatial extent as the
  /// original template.
  ///@param float ProfileHalfWidth
  ///@author DTM
  ///@date 11/04/2014
  ///
  LSDSwath(PointData& ProfilePoints, LSDRaster& RasterTemplate, float HalfWidth) { create(ProfilePoints, RasterTemplate, HalfWidth); }
	
	void get_transverse_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles);
	void get_longitudinal_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles);
  // write profiles to file
  void write_transverse_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth, string prefix);
  void write_longitudinal_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth, string prefix);
  // get functions
  // these get data elements
  int get_NPtsInProfile() const {return NPtsInProfile;}
  Array2D<float> get_DistanceToBaselineArray() const { return DistanceToBaselineArray; }
  Array2D<float> get_DistanceAlongBaselineArray() const { return DistanceAlongBaselineArray; }
  
	protected:
  
  // Swath template
  vector<float> DistanceAlongBaseline;
  Array2D<float> DistanceToBaselineArray;
  Array2D<float> DistanceAlongBaselineArray;

	// metadata
  int NPtsInProfile;
  float ProfileHalfWidth;
  float NoDataValue;
  int NRows;
  int NCols;
  
  // Bounding Box of profile
  float XMax;
  float XMin;
  float YMax;
  float YMin;

	private:
  void create();
  void create(PointData& ProfilePoints, LSDRaster& RasterTemplate, float ProfileHalfWidth);
            
};

#endif
