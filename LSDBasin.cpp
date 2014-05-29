#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDIndexChannel.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
#include "LSDBasin.hpp"

using namespace std;
using namespace TNT;

#ifndef LSDBasin_CPP
#define LSDBasin_CPP


void LSDBasin::create(int JunctionNumber, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet){

  //NO BOUNDS CHECKING ON JunctionNumber

  //setting all of the instance variables for the given junction

  NRows = ChanNet.get_NRows();
	NCols = ChanNet.get_NCols();
	XMinimum = ChanNet.get_XMinimum();
	YMinimum = ChanNet.get_YMinimum();
	DataResolution = ChanNet.get_DataResolution();
	NoDataValue = ChanNet.get_NoDataValue();

  Junction = JunctionNumber;
  
  vector <int> JunctionVector = ChanNet.get_JunctionVector();
  vector <int> ReceiverVector = ChanNet.get_ReceiverVector();
  
  LSDIndexChannel StreamLinkVector = LSDIndexChannel(Junction, JunctionVector[Junction],
                                                     ReceiverVector[Junction], JunctionVector[ReceiverVector[Junction]], FlowInfo);

  int n_nodes_in_channel = StreamLinkVector.get_n_nodes_in_channel();
  int basin_outlet = StreamLinkVector.get_node_in_channel(n_nodes_in_channel-2);
  BasinNodes = FlowInfo.get_upslope_nodes(basin_outlet);
                                                                                     
  NumberOfCells = int(BasinNodes.size());
  Area = NumberOfCells * (DataResolution*DataResolution);
  
  Beheaded = ChanNet.node_tester(FlowInfo, Junction);

  FlowInfo.retrieve_current_row_and_col(ChanNet.get_Node_of_Junction(Junction), Outlet_i, Outlet_j);
    
  vector<int> StreamOrderVector = ChanNet.get_StreamOrderVector();
  
  BasinOrder = StreamOrderVector[Junction];


  int i_max = 0;
  int i_min = 9999999; //a very large number
  int j_max = 0;
  int j_min = 9999999; //a very large number
  
  int i = 0;
  int j = 0;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    
    if (i > i_max){i_max = i;}
    else if (i < i_min){i_min = i;}
    if (j > j_max){j_max = j;}
    else if (j < j_min){j_min = j;}
    
  }
  
  Centroid_i = i_min + ((i_max - i_min)/2);
  Centroid_j = j_min + ((j_max - j_min)/2);   //how do these handle 0.5s ??


  //finished setting all the instance variables
  
  
  // now we set up empty variables to store properties of the basin
  // these are populated as they are required using the set methods in LSDBasin
    
  SlopeMean = NoDataValue;
  ElevationMean = NoDataValue;
  AspectMean = NoDataValue;
  ReliefMean = NoDataValue;
  PlanCurvMean = NoDataValue;
  ProfileCurvMean = NoDataValue;
  TotalCurvMean = NoDataValue;
  PlanCurvMax = NoDataValue;
  ProfileCurvMax = NoDataValue;
  TotalCurvMax = NoDataValue;
  HillslopeLength_HFR = NoDataValue;
  HillslopeLength_Binned = NoDataValue;
  HillslopeLength_Spline = NoDataValue;
  HillslopeLength_Density = NoDataValue;
  FlowLength = NoDataValue;
  DrainageDensity = NoDataValue;  
  Perimeter_i = vector<int>(1,NoDataValue);
  Perimeter_j =  vector<int>(1,NoDataValue);
  CosmoErosionRate = NoDataValue;
  OtherErosionRate = NoDataValue;
  CHTMean = NoDataValue;
  EStar = NoDataValue;
  RStar = NoDataValue;
   
  //finished creating empty variables 

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate mean basin value.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMean(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  float TotalData = 0;
  int CountNDV = 0;
  float BasinAverage;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    
    //exclude NDV from average
    if (Data.get_data_element(i,j) != NoDataValue){
      TotalData += Data.get_data_element(i,j);
    }
    else {
      ++CountNDV;
    }
  }

  BasinAverage = TotalData/(NumberOfCells-CountNDV);

  return BasinAverage;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate max basin value.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMax(LSDFlowInfo& FlowInfo, LSDRaster Data){

  //could use max_element here? how would that cope with NDVs??

  int i;
  int j;
  float MaxData = -10000000;   //a very small number
  float CurrentData;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    CurrentData = Data.get_data_element(i,j);
    
    //exclude NDV
    if (CurrentData != NoDataValue && CurrentData > MaxData){
      MaxData = CurrentData;     
    }
  }
    
  return MaxData;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate min basin value.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMin(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  float MinData = 100000000; // a large number
  float CurrentData;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    CurrentData = Data.get_data_element(i,j);
    
    //exclude NDV
    if (CurrentData != NoDataValue && CurrentData < MinData){
      MinData = CurrentData;     
    }
  }
    
  return MinData;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate median basin value.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinMedian(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  vector<float> UnsortedData;
  vector<float> SortedData;
  vector<size_t> index_map;
  float Median;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    //exclude NDV
    if (Data.get_data_element(i,j) != NoDataValue){
      UnsortedData.push_back(Data.get_data_element(i,j));     
    }
  }
  
  //get size of dataset
  size_t n = UnsortedData.size() / 2;
  
  //sort all non NDV values
  matlab_float_sort(UnsortedData, SortedData, index_map);
  
  if (n % 2 != 0){
    Median = (SortedData[n] + SortedData[n+1])/2;
  }
  else{
    Median = SortedData[n];
  }
    
  return Median;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate Standard devaition of the basin values.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinStdDev(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  vector<float> DataValues;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    //exclude NDV
    if (Data.get_data_element(i,j) != NoDataValue){
      DataValues.push_back(Data.get_data_element(i,j));     
    }
  }
  
  float mean = get_mean(DataValues);
  float StdDev = get_standard_deviation(DataValues, mean);
    
  return StdDev;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate Standard Error of the basin values.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinStdError(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  vector<float> DataValues;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    //exclude NDV
    if (Data.get_data_element(i,j) != NoDataValue){
      DataValues.push_back(Data.get_data_element(i,j));     
    }
  }
  
  float mean = get_mean(DataValues);
  float StdDev = get_standard_deviation(DataValues, mean);
  float StdError = get_standard_error(DataValues, StdDev);
    
  return StdError;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate basin range.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
float LSDBasin::CalculateBasinRange(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  float MinData = 100000000; // a large number
  float MaxData = -100000000; // a small number
  float CurrentData;

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    CurrentData = Data.get_data_element(i,j);
    
    //exclude NDV
    if (CurrentData != NoDataValue && CurrentData < MinData){
      MinData = CurrentData;     
    }
    if (CurrentData != NoDataValue && CurrentData > MaxData){
      MaxData = CurrentData;     
    }
  }
  float Range = MaxData - MinData;     
  return Range;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate basin range.
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDBasin::CalculateNumDataPoints(LSDFlowInfo& FlowInfo, LSDRaster Data){

  int i;
  int j;
  int count = 0;
  
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
        
    //exclude NDV
    if (Data.get_data_element(i,j) != NoDataValue){
      ++count;     
    }
  }
      
  return count;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate E* and R* values for the basin, using hilltop flow routed hillslope 
// lengths. 
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_EStar_RStar(float CriticalSlope){

    EStar = (2 * (abs(CHTMean)) * HillslopeLength_HFR) / CriticalSlope;
    RStar = ReliefMean / (HillslopeLength_HFR * CriticalSlope);
    
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate flow length for the basin using the D8 flow directions. 
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_FlowLength(LSDIndexRaster& StreamNetwork, LSDFlowInfo& FlowInfo){

  int j;
  int i;
  float LengthSum = 0;
  float two_times_root2 = 2.828427;
  Array2D<int> FlowDir = FlowInfo.get_FlowDirection();


  //Loop over every pixel and record it's stream length and basin ID in two vectors  
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
           
     if (StreamNetwork.get_data_element(i,j) != NoDataValue){
       if ((FlowDir[i][j] % 2) != 0 && (FlowDir[i][j] != -1 )){ //is odd but not -1
         LengthSum += (DataResolution * two_times_root2); //diagonal
       }
       else if (FlowDir[i][j] % 2 == 0){  //is even
         LengthSum +=  DataResolution; //cardinal                   
       }
     }
  }

  FlowLength = LengthSum;
 
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Calculate hillslope lengths from boomerang plots. 
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_HillslopeLengths_Boomerang(LSDRaster& Slope, LSDRaster& DinfArea, LSDFlowInfo& FlowInfo, float log_bin_width, int SplineResolution, float bin_threshold){
  
  int j;
  int i;
  Array2D<float> slope(NRows, NCols, NoDataValue);
  Array2D<float> area(NRows, NCols, NoDataValue);
  
  //create subset arrays for just the basin data - this should be rolled into its own method.
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    
      slope[i][j] = Slope.get_data_element(i,j);
      area[i][j] = DinfArea.get_data_element(i,j);
    
  }

  //do some log binning
  vector<float> Mean_x_out;
  vector<float> Mean_y_out;
  vector<float> Midpoints_out;
  vector<float> STDDev_x_out;
  vector<float> STDDev_y_out;
  vector<float> STDErr_x_out;
  vector<float> STDErr_y_out;
  vector<int> number_observations;
  
  log_bin_data(area, slope, log_bin_width, Mean_x_out, Mean_y_out, Midpoints_out, STDDev_x_out, STDDev_y_out, STDErr_x_out, STDErr_y_out, number_observations, NoDataValue);  
  
  //remove empty bins 
  RemoveSmallBins(Mean_x_out, Mean_y_out, Midpoints_out, STDDev_x_out, STDDev_y_out, STDErr_x_out, STDErr_y_out, number_observations, bin_threshold);
  
  //index value of max slope
  int slope_max_index = distance(Mean_y_out.begin(), max_element(Mean_y_out.begin(), Mean_y_out.end()));

  //hillslope length from the maximum binned values
  HillslopeLength_Binned = Mean_x_out[slope_max_index]/DataResolution;
      
  // Fit splines through the binned data to get the LH
  vector<float> Spline_X;
  vector<float> Spline_Y;
  PlotCubicSplines(Mean_x_out, Mean_y_out, SplineResolution, Spline_X, Spline_Y);

  //index value of max spline slope
  int slope_max_index_spline = distance(Spline_Y.begin(), max_element(Spline_Y.begin(), Spline_Y.end()));

  //hillslope length from spline curve
  HillslopeLength_Spline = Spline_X[slope_max_index_spline]/DataResolution;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Generate data to create boomerang plots. 
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::Plot_Boomerang(LSDRaster& Slope, LSDRaster& DinfArea, LSDFlowInfo& FlowInfo, float log_bin_width, int SplineResolution, float bin_threshold, string Path){
  
  int j;
  int i;
  Array2D<float> slope(NRows, NCols, NoDataValue);
  Array2D<float> area(NRows, NCols, NoDataValue);
  
  //create subset arrays for just the basin data - this should be rolled into its own method.
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    
      slope[i][j] = Slope.get_data_element(i,j);
      area[i][j] = DinfArea.get_data_element(i,j);
    
  }

  //do some log binning
  vector<float> Mean_x_out;
  vector<float> Mean_y_out;
  vector<float> Midpoints_out;
  vector<float> STDDev_x_out;
  vector<float> STDDev_y_out;
  vector<float> STDErr_x_out;
  vector<float> STDErr_y_out;
  vector<int> number_observations;
  
  log_bin_data(area, slope, log_bin_width, Mean_x_out, Mean_y_out, Midpoints_out, STDDev_x_out, STDDev_y_out, STDErr_x_out, STDErr_y_out, number_observations, NoDataValue);  
  
  //remove empty bins 
  RemoveSmallBins(Mean_x_out, Mean_y_out, Midpoints_out, STDDev_x_out, STDDev_y_out, STDErr_x_out, STDErr_y_out, number_observations, bin_threshold);
        
  // Fit splines through the binned data to get the LH
  vector<float> Spline_X;
  vector<float> Spline_Y;
  PlotCubicSplines(Mean_x_out, Mean_y_out, SplineResolution, Spline_X, Spline_Y);

  //set up a filestream object to write the binned data
  ofstream file;

  stringstream ss_bin;
  ss_bin << Path << Junction << "_boom_binned.txt";
  file.open(ss_bin.str().c_str());   //needs a null terminated character array, not a string. See pg 181 of accelerated c++
    
  for(int q = 0; q < int(Mean_x_out.size()); q++){
    file << Mean_x_out[q] << " " << Mean_y_out[q] << " " << STDDev_x_out[q] << " " << STDDev_y_out[q] << " " << STDErr_x_out[q] << " " << STDErr_y_out[q] << endl;
  }
  file.close();
      
  //set up a filestream object to write the spline data
  ofstream SplineFile;

  stringstream ss_spline;
  ss_spline << Path << Junction << "_boom_spline.txt";
  SplineFile.open(ss_spline.str().c_str());   //needs a null terminated character array, not a string. See pg 181 of accelerated c++
    
  for(int q = 0; q < int(Spline_X.size()); q++){ //fixed bug here where I looped over the wrong vector - SWDG 7/11/13
    SplineFile << Spline_X[q] << " " << Spline_Y[q] << endl;

  }
  SplineFile.close();
  
  //set up a filestream object to write the data cloud
  ofstream cloud;

  stringstream ss_cloud;
  ss_cloud << Path << Junction << "_boom_cloud.txt";
  cloud.open(ss_cloud.str().c_str());     //needs a null terminated character array, not a string. See pg 181 of accelerated c++

  for (int i = 1; i < NRows-1; ++i){
    for (int j = 1; j < NCols-1; ++j){
      if(area[i][j] != NoDataValue && slope[i][j] != NoDataValue){
        cloud << area[i][j] << " " << slope[i][j] << endl;
      }
    }
  }
  cloud.close();

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set the mean basin aspect. Does not use the normal basin mean method as angles
// need to be handled differently. 
// Bug fixed in the average calculation when values wrapped around 0
// SWDG 17/2/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_AspectMean(LSDFlowInfo& FlowInfo, LSDRaster Aspect){

  int i;
  int j;
  float avg_r;
  float angle_r;
  float x_component = 0.0;
  float y_component = 0.0;
  int ndv_cell_count = 0;  

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    
    if (Aspect.get_data_element(i,j) != NoDataValue){
    
      angle_r = rad(Aspect.get_data_element(i,j));
      x_component += cos(angle_r);
      y_component += sin(angle_r);
  
    }
    else{
      ++ndv_cell_count;
    }
  
  }
    
  avg_r = atan2(y_component, x_component);
  AspectMean = deg(avg_r);
  
  if (AspectMean < 0){
    AspectMean = 360 + AspectMean;
  }
   
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set the perimeter pixels using a simple edge detection algorithm. This is quite 
// messy and will be improved soon.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_Perimeter(LSDFlowInfo& FlowInfo){

  int i;
  int j;
  vector<int> I;
  vector<int> J;
  int NDVCount = 0;
  Array2D<float> BasinData(NRows, NCols, NoDataValue);

  //create subset arrays for just the basin data - this should be rolled into its own method.
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
      BasinData[i][j] = BasinNodes[q];
    
  }

  for (int q = 0; q < int(BasinNodes.size()); ++q){
    
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    
      NDVCount = 0;
     
        //count border cells that are NDV
        if (BasinData[i-1][j-1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i][j-1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i+1][j-1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i-1][j] == NoDataValue){ ++NDVCount; }
        if (BasinData[i+1][j] == NoDataValue){ ++NDVCount; }
        if (BasinData[i-1][j+1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i][j+1] == NoDataValue){ ++NDVCount; }
        if (BasinData[i+1][j+1] == NoDataValue){ ++NDVCount; }
        
        if (NDVCount >= 4 && NDVCount < 8){  //increase the first value to get a simpler polygon
          //edge pixel
          I.push_back(i);
          J.push_back(j);
        }
    
  }

  //now have 2 vectors of i and j indexes of every point
  Perimeter_i = I;
  Perimeter_j = J;


}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set the four different hillslope length measurements for the basin. 
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_all_HillslopeLengths(LSDFlowInfo& FlowInfo, LSDRaster& HillslopeLengths, LSDRaster& Slope, LSDRaster& DinfArea, float log_bin_width, int SplineResolution, float bin_threshold){

  set_HillslopeLength_HFR(FlowInfo, HillslopeLengths);
  set_HillslopeLengths_Boomerang(Slope, DinfArea, FlowInfo, log_bin_width, SplineResolution, bin_threshold);

  if (DrainageDensity != NoDataValue){ 
    set_HillslopeLength_Density();
  }
  else{
    cout << "\nDrainage Density has not been set, so the hillslope length cannot be set." << endl;
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Set all of the basin parameters with one call.
//
// Runs polyfit to get the elevation derivatives, so can be quite memory intensive. Method
// calls all the setters one by one, to populate all the basin parameters. So a
// basin can be created and all it's properties set with 2 calls. The erosion rates have default 
// parameters of -9999 as these are rarely used variables.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDBasin::set_All_Parameters(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDRaster& CHT, LSDIndexRaster& StreamNetwork,
                                  LSDRaster& HillslopeLengths, LSDRaster& Relief, float window_radius, float log_bin_width,
                                  int SplineResolution, float bin_threshold, float CriticalSlope, float CosmoErosionRate, 
                                  float OtherErosionRate){

  // coefficent matrices for polyfit routine
  Array2D<float> a;
  Array2D<float> b;
  Array2D<float> c;
  Array2D<float> d;
  Array2D<float> e;
  Array2D<float> f;

  Elevation.calculate_polyfit_coefficient_matrices(window_radius, a, b, c, d, e, f);
  LSDRaster TotalCurv = Elevation.calculate_polyfit_curvature (a,b);
  LSDRaster ProfileCurv = Elevation.calculate_polyfit_profile_curvature (a,b,c,d,e);
  LSDRaster PlanCurv = Elevation.calculate_polyfit_planform_curvature (a,b,c,d,e);
  LSDRaster Aspect = Elevation.calculate_polyfit_aspect(d,e);  
  LSDRaster Slope = Elevation.calculate_polyfit_slope(d,e);
  LSDRaster DinfArea = Elevation.D_inf_units(); 
  
  set_SlopeMean(FlowInfo, Slope);
  set_ElevationMean(FlowInfo, Elevation);
  set_ReliefMean(FlowInfo, Relief);
  set_PlanCurvMean(FlowInfo, PlanCurv);
  set_ProfileCurvMean(FlowInfo, ProfileCurv);
  set_TotalCurvMean(FlowInfo, TotalCurv);
  set_PlanCurvMax(FlowInfo, PlanCurv);
  set_ProfileCurvMax(FlowInfo, ProfileCurv);
  set_TotalCurvMax(FlowInfo, TotalCurv);
  set_CHTMean(FlowInfo, CHT);
  set_AspectMean(FlowInfo, Aspect);
  set_FlowLength(StreamNetwork, FlowInfo);
  set_DrainageDensity();
  set_all_HillslopeLengths(FlowInfo, HillslopeLengths, Slope, DinfArea, log_bin_width, SplineResolution, bin_threshold);
  set_Perimeter(FlowInfo);
  set_EStar_RStar(CriticalSlope);
  set_CosmoErosionRate(CosmoErosionRate);
  set_OtherErosionRate(OtherErosionRate);

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Write integer basin parameters into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDBasin::write_integer_data_to_LSDIndexRaster(int Param, LSDFlowInfo FlowInfo){
  
  int i;
  int j; 
  Array2D<int> Output(NRows, NCols, NoDataValue);
  
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    Output[i][j] = Param;
  }

  LSDIndexRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output);

  return OutputRaster;

}                                       

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Write real basin parameters into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDBasin::write_real_data_to_LSDRaster(float Param, LSDFlowInfo FlowInfo){
  
  int i;
  int j; 
  Array2D<float> Output(NRows, NCols, NoDataValue);
  
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    Output[i][j] = Param;
  }

  LSDRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Cookie cut data from an LSDRaster into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDRaster LSDBasin::write_raster_data_to_LSDRaster(LSDRaster Data, LSDFlowInfo FlowInfo){
  
  int i;
  int j; 
  Array2D<float> Output(NRows, NCols, NoDataValue);
  
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    Output[i][j] = Data.get_data_element(i,j);
  }

  LSDRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Cookie cut data from an LSDIndexRaster into the shape of the basin.
// SWDG 12/12/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDBasin::write_raster_data_to_LSDIndexRaster(LSDIndexRaster Data, LSDFlowInfo FlowInfo){
  
  int i;
  int j; 
  Array2D<int> Output(NRows, NCols, NoDataValue);
  
  for (int q = 0; q < int(BasinNodes.size()); ++q){
    FlowInfo.retrieve_current_row_and_col(BasinNodes[q], i, j);
    Output[i][j] = Data.get_data_element(i,j);
  }

  LSDIndexRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution,
                               NoDataValue, Output);

  return OutputRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Function to check whether a node is in the basin
// FJC 21/02/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int LSDBasin::is_node_in_basin(int test_node)
{
  int node_checker = 0;
  for (int i =0; i < int(BasinNodes.size()); i++)
  {
    if (test_node == BasinNodes[i])
    {
      node_checker = 1;
    }
  }  
  return node_checker;
}

#endif