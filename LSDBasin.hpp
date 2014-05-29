//beginning of the LSDBasin object

#include <vector>
#include <string>
#include <algorithm>
#include "TNT/tnt.h"
#include "LSDFlowInfo.hpp"
#include "LSDRaster.hpp"
#include "LSDIndexRaster.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDStatsTools.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDBasin_H
#define LSDBasin_H

///@brief Object to store information about drainage basins and generate basin average metrics..
class LSDBasin
{

  public:
  
  /// @brief Default constructor method used to create a basin object.
  /// @param Junction outlet junction of the basin to be constructed.
  /// @param FlowInfo LSDFlowInfo object.
  /// @param ChanNet Channel network object.
  /// @author SWDG
  /// @date 11/12/12
  LSDBasin(int Junction, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet)
											{ create(Junction, FlowInfo, ChanNet); }


  /// @return Number of rows as an integer.
	int get_NRows() const				{ return NRows; }
	/// @return Number of columns as an integer.
  int get_NCols() const				{ return NCols; }
  /// @return Minimum X coordinate as an integer.
	float get_XMinimum() const			{ return XMinimum; }
	/// @return Minimum Y coordinate as an integer.
	float get_YMinimum() const			{ return YMinimum; }
	/// @return Data resolution as an integer.
	float get_DataResolution() const	{ return DataResolution; }
	/// @return No Data Value as an integer.
	int get_NoDataValue() const			{ return NoDataValue; }
	
	/// @return Nodes of the basin.
	vector<int> get_BasinNodes() const { return BasinNodes; }
	
	/// @return Number of cells in the basin.
	int get_NumberOfCells() const { return NumberOfCells; }
	
  /// @return Area of basin in spataial units.
  float get_Area() const { return Area; }
  
  /// @return Junction Number of the basin.
  int get_Junction() const { return Junction; }
  
  /// @return Stream order of basin.
  int get_BasinOrder() const { return BasinOrder; }
  
  /// @return Boolean value of whether a basin is beheaded or not.										
  bool get_Beheaded() const { return Beheaded; }
  
  /// @return i index of outlet pixel.
  int get_Outlet_i() const { return Outlet_i; }

  /// @return j index of outlet pixel.
  int get_Outlet_j() const { return Outlet_j; }
  
  ///@return i index of Centroid pixel.
  int get_Centroid_i() const { return Centroid_i; }

  ///@return j index of Centroid pixel.
  int get_Centroid_j() const { return Centroid_j; }
  
  //Getters of basin parameters
  
  /// @return Mean slope.
  float get_SlopeMean() const { return SlopeMean; }
  /// @return Mean elevation.
  float get_ElevationMean() const { return ElevationMean; }
  /// @return Mean aspect.
  float get_AspectMean() const { return AspectMean; }
  /// @return Mean relief.
  float get_ReliefMean() const { return ReliefMean; }
  /// @return Mean plan curvature.
  float get_PlanCurvMean() const { return PlanCurvMean; }
  /// @return Mean profile curvature.
  float get_ProfileCurvMean() const { return ProfileCurvMean; }
  /// @return Mean total curvature.
  float get_TotalCurvMean() const { return TotalCurvMean; }
  /// @return Max plan curvature.
  float get_PlanCurvMax() const { return PlanCurvMax; }
  /// @return Max profile curvature.
  float get_ProfileCurvMax() const { return ProfileCurvMax; }
  /// @return Max total curvature.
  float get_TotalCurvMax() const { return TotalCurvMax; }
  /// @return Hillslope length from hilltop flow routing.
  float get_HillslopeLength_HFR() const { return HillslopeLength_HFR; }
  /// @return Hillslope length from boomerang bins.
  float get_HillslopeLength_Binned() const { return HillslopeLength_Binned; }
  /// @return Hillslope length from boomerang splines.
  float get_HillslopeLength_Spline() const { return HillslopeLength_Spline; }
  /// @return Hillslope length from drainage density.
  float get_HillslopeLength_Density() const { return HillslopeLength_Density; }
  /// @return Flow length.
  float get_FlowLength() const { return FlowLength; }
  /// @return Drainage Density.
  float get_DrainageDensity() const { return DrainageDensity; }
  /// @return Basin Perimeter's i index.
  vector<int> get_Perimeter_i() const { return Perimeter_i; }
  /// @return Basin Perimeter's j index.
  vector<int> get_Perimeter_j() const { return Perimeter_j; }
  /// @return Cosmo erosion rate.
  float get_CosmoErosionRate() const { return CosmoErosionRate; }
  /// @return Other eroision rate.
  float get_OtherErosionRate() const { return OtherErosionRate; }
  /// @return Mean hilltop curvature.
  float get_CHTMean() const { return CHTMean; }
  /// @return E* value.
  float get_EStar() const { return EStar; }
  /// @return R* value.
  float get_RStar() const { return RStar; }  
  
  /// @brief Calculate the mean value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the mean of.
  /// @return Mean value.
  /// @author SWDG
  /// @date 11/12/13
  float CalculateBasinMean(LSDFlowInfo& FlowInfo, LSDRaster Data);
  
  /// @brief Calculate the max value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the max of.
  /// @return Max value.
  /// @author SWDG
  /// @date 11/12/13
  float CalculateBasinMax(LSDFlowInfo& FlowInfo, LSDRaster Data);
  
  /// @brief Calculate the min value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the minimum of.
  /// @return Min value.
  /// @author SWDG
  /// @date 17/2/14
  float CalculateBasinMin(LSDFlowInfo& FlowInfo, LSDRaster Data);  
  
  /// @brief Calculate the median value of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the median of.
  /// @return Median value.
  /// @author SWDG
  /// @date 17/2/14
  float CalculateBasinMedian(LSDFlowInfo& FlowInfo, LSDRaster Data);

  /// @brief Calculate the Standard Deviation of values of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the standard deviation of.
  /// @return Standard deviation value.
  /// @author SWDG
  /// @date 17/2/14  
  float CalculateBasinStdDev(LSDFlowInfo& FlowInfo, LSDRaster Data);
  
  /// @brief Calculate the Standard error of values of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the standard error of.
  /// @return Standard error value.
  /// @author SWDG
  /// @date 17/2/14  
  float CalculateBasinStdError(LSDFlowInfo& FlowInfo, LSDRaster Data);
  
  /// @brief Calculate the range of values of an LSDRaster which falls inside a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to find the range of.
  /// @return Range value.
  /// @author SWDG
  /// @date 17/2/14
  float CalculateBasinRange(LSDFlowInfo& FlowInfo, LSDRaster Data);
  
  /// @brief Calculate the number of data points of an LSDRaster which fall inside a basin.
  ///
  /// @details Useful for checking that an average value is not just taken from a single data point.
  /// @param FlowInfo Flowinfo object.
  /// @param Data Values to count.
  /// @return Number of data points.
  /// @author SWDG
  /// @date 17/2/14
  int CalculateNumDataPoints(LSDFlowInfo& FlowInfo, LSDRaster Data);
  
  /// @brief Set the mean slope of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Slope Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_SlopeMean(LSDFlowInfo& FlowInfo, LSDRaster Slope){ SlopeMean = CalculateBasinMean(FlowInfo, Slope); }

  /// @brief Set the mean Elevation of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Elevation Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_ElevationMean(LSDFlowInfo& FlowInfo, LSDRaster Elevation) { ElevationMean = CalculateBasinMean(FlowInfo, Elevation); }

  /// @brief Set the mean Relief of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Relief Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_ReliefMean(LSDFlowInfo& FlowInfo, LSDRaster Relief) { ReliefMean = CalculateBasinMean(FlowInfo, Relief); }

  /// @brief Set the mean PlanCurve of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param PlanCurv Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_PlanCurvMean(LSDFlowInfo& FlowInfo, LSDRaster PlanCurv) { PlanCurvMean = CalculateBasinMean(FlowInfo, PlanCurv); }

  /// @brief Set the mean ProfCurv of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param ProfileCurv Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_ProfileCurvMean(LSDFlowInfo& FlowInfo, LSDRaster ProfileCurv) { ProfileCurvMean = CalculateBasinMean(FlowInfo, ProfileCurv); }

  /// @brief Set the mean TotalCurv of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param TotalCurv Values to find the mean of.
  /// @author SWDG
  /// @date 11/12/13
  void set_TotalCurvMean(LSDFlowInfo& FlowInfo, LSDRaster TotalCurv) { TotalCurvMean = CalculateBasinMean(FlowInfo, TotalCurv); }

  /// @brief Set the max PlanCurve of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param PlanCurv Values to find the max of.
  /// @author SWDG
  /// @date 11/12/13
  void set_PlanCurvMax(LSDFlowInfo& FlowInfo, LSDRaster PlanCurv) { PlanCurvMax = CalculateBasinMax(FlowInfo, PlanCurv); }

  /// @brief Set the max ProfCurv of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param ProfileCurv Values to find the max of.
  /// @author SWDG
  /// @date 11/12/13
  void set_ProfileCurvMax(LSDFlowInfo& FlowInfo, LSDRaster ProfileCurv) { ProfileCurvMax = CalculateBasinMax(FlowInfo, ProfileCurv); }

  /// @brief Set the max TotalCurv of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param TotalCurv Values to find the max of.
  /// @author SWDG
  /// @date 11/12/13
  void set_TotalCurvMax(LSDFlowInfo& FlowInfo, LSDRaster TotalCurv) { TotalCurvMax = CalculateBasinMax(FlowInfo, TotalCurv); }
 
  /// @brief Set the mean hilltop curvature of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param CHT Values to find the mean of.
  /// @author SWDG
  /// @date 12/12/13
  void set_CHTMean(LSDFlowInfo& FlowInfo, LSDRaster CHT) { CHTMean = CalculateBasinMean(FlowInfo, CHT); }

  /// @brief Set the Cosmogenic erosion rate.
  /// @param ErosionRate Erosion rate - No sanity check on this value.
  /// @author SWDG
  /// @date 11/12/13
  void set_CosmoErosionRate(float ErosionRate) { CosmoErosionRate = ErosionRate; }
  
  /// @brief Set the Other erosion rate.
  /// @param ErosionRate Erosion rate - No sanity check on this value.
  /// @author SWDG
  /// @date 11/12/13
  void set_OtherErosionRate(float ErosionRate) { OtherErosionRate = ErosionRate; }
 
  /// @brief Calculate E* and R* values for the basin, using hilltop flow routed hillslope lengths. 
  /// @param CriticalSlope slope threshold value, typically 0.4.
  /// @author SWDG
  /// @date 12/12/13
  void set_EStar_RStar(float CriticalSlope);
 
  /// @brief Calculate flow length for the basin using the D8 flow directions. 
  /// @param StreamNetwork the channel network.
  /// @param FlowInfo Flowinfo object.
  /// @author SWDG
  /// @date 12/12/13
  void set_FlowLength(LSDIndexRaster& StreamNetwork, LSDFlowInfo& FlowInfo);

  /// @brief Set basin drainage density.
  /// @author SWDG
  /// @date 12/12/13
  void set_DrainageDensity() { DrainageDensity = FlowLength/Area ; }
  
  /// @brief Set the mean HillslopeLength from hilltop flow routing.
  /// @param FlowInfo Flowinfo object.
  /// @param HillslopeLengths Values to find the mean of.
  /// @author SWDG
  /// @date 12/12/13
  void set_HillslopeLength_HFR(LSDFlowInfo& FlowInfo, LSDRaster HillslopeLengths) { HillslopeLength_HFR = CalculateBasinMean(FlowInfo, HillslopeLengths); }
 
  /// @brief Set mean HillslopeLengths from boomerang plots from both splines and binned data.
  /// @param Slope LSDRaster of slope.
  /// @param DinfArea D-infinity Flowarea LSDRaster.
  /// @param FlowInfo Flowinfo object.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf.
  /// @param SplineResolution Number of values between each point for the spline curve.
  /// @param bin_threshold Threshold fraction of values needed to keep a bin. 
  /// @author SWDG
  /// @date 12/12/13
  void set_HillslopeLengths_Boomerang(LSDRaster& Slope, LSDRaster& DinfArea, LSDFlowInfo& FlowInfo, float log_bin_width, int SplineResolution, float bin_threshold);
 
  /// @brief Set the mean HillslopeLength from drainage density.
  /// @author SWDG
  /// @date 12/12/13
  void set_HillslopeLength_Density() { HillslopeLength_Density = (1 / (2 * DrainageDensity)); }
  
  /// @brief Generate text files containing data to plot boomerangs.
  ///
  /// @details Writes 3 files to the output path, coded with the basin's unique
  /// junction number which can the be read with python and plotted.
  /// @param Slope LSDRaster of slope.
  /// @param DinfArea D-infinity Flowarea LSDRaster.
  /// @param FlowInfo Flowinfo object.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf.
  /// @param SplineResolution Number of values between each point for the spline curve.
  /// @param bin_threshold Threshold fraction of values needed to keep a bin.
  /// @param Path The output path where the data files will be written to, including the final slash. 
  /// @author SWDG
  /// @date 12/12/13  
  void Plot_Boomerang(LSDRaster& Slope, LSDRaster& DinfArea, LSDFlowInfo& FlowInfo, float log_bin_width, int SplineResolution, float bin_threshold, string Path);
  
  /// @brief Set the mean aspect of a basin.
  /// @param FlowInfo Flowinfo object.
  /// @param Aspect Values to find the mean of.
  /// @author SWDG
  /// @date 17/2/14
  void set_AspectMean(LSDFlowInfo& FlowInfo, LSDRaster Aspect);
  
  /// @brief Set the perimeter pixels using a simple edge detection algorithm. 
  ///
  /// @details This is quite messy and will be improved soon.
  /// @param FlowInfo Flowinfo object.
  /// @author SWDG
  /// @date 12/12/13
  void set_Perimeter(LSDFlowInfo& FlowInfo);
  
  /// @brief Set the four different hillslope length measurements for the basin.
  /// @param FlowInfo Flowinfo object.
  /// @param HillslopeLengths LSDRaster of hillslope lengths from the hilltop flow routing method.
  /// @param Slope LSDRaster of slope.
  /// @param DinfArea D-infinity Flowarea LSDRaster.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf.
  /// @param SplineResolution Number of values between each point for the spline curve.
  /// @param bin_threshold Threshold fraction of values needed to keep a bin.
  /// @author SWDG
  /// @date 12/12/13
  void set_all_HillslopeLengths(LSDFlowInfo& FlowInfo, LSDRaster& HillslopeLengths, LSDRaster& Slope, 
                                LSDRaster& DinfArea, float log_bin_width, int SplineResolution, float bin_threshold);  

  /// @brief Set all of the basin parameters with one call.
  ///
  /// @details Runs polyfit to get the elevation derivatives, so can be quite memory intensive. Method
  /// calls all the setters one by one, to populate all the basin parameters. So a
  /// basin can be created and all it's properties set with 2 calls. The erosion rates have default 
  /// parameters of -9999 as these are rarely used variables.
  /// @param Elevation LSDRaster of filled elevation values.
  /// @param FlowInfo Flowinfo object.
  /// @param CHT LSDRaster of hilltop curvatures.
  /// @param StreamNetwork LSDIndexRaster of the stream network.
  /// @param HillslopeLengths LSDRaster of hillslope lengths from the hilltop flow routing method.
  /// @param Relief LSDRaster of the hilltop relief.
  /// @param window_radius Radius in spatial units for the polyft routine.
  /// @param log_bin_width Width (in log space) of the bins, with respect to D_inf. Default value is 0.1.
  /// @param SplineResolution Number of values between each point for the spline curve. Default value is 10000.
  /// @param bin_threshold Threshold fraction of values needed to keep a bin. Default value is 0.
  /// @param CriticalSlope Slope threshold used for E* R* calculations. Default value is 0.4.
  /// @param CosmoErosionRate Erosion rate from cosmo.
  /// @param OtherErosionRate Erosion rate from another source.
  /// @author SWDG
  /// @date 12/12/13
  void set_All_Parameters(LSDRaster& Elevation, LSDFlowInfo& FlowInfo, LSDRaster& CHT, LSDIndexRaster& StreamNetwork,
                          LSDRaster& HillslopeLengths, LSDRaster& Relief, float window_radius, float log_bin_width, 
                          int SplineResolution, float bin_threshold, float CriticalSlope, 
                          float CosmoErosionRate = -9999, float OtherErosionRate = -9999);

  /// @brief Cookie cut data from an LSDIndexRaster into the shape of the basin.
  /// @param Data LSDIndexRaster data to be written.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of the data in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDIndexRaster write_raster_data_to_LSDIndexRaster(LSDIndexRaster Data, LSDFlowInfo FlowInfo);
  
  /// @brief Cookie cut data from an LSDRaster into the shape of the basin.
  /// @param Data LSDRaster data to be written.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of the data in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_raster_data_to_LSDRaster(LSDRaster Data, LSDFlowInfo FlowInfo);
  
  /// @brief check whether a test node is in the basin or not
  /// @param test_node node to test
  /// @return integer which is 1 if node is in the basin or 0 if it is not
  /// @author FJC
  /// @date 21/02/14
  int is_node_in_basin(int test_node);
  
  /// @brief Write a real value to an LSDRaster in the shape of the basin.
  /// @param Param real value to be written
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of the data in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_real_data_to_LSDRaster(float Param, LSDFlowInfo FlowInfo);
 
  /// @brief Write an integer value to an LSDIndexRaster in the shape of the basin.
  /// @param Param integer value to be written
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of the data in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDIndexRaster write_integer_data_to_LSDIndexRaster(int Param, LSDFlowInfo FlowInfo);
  
  /// @brief Write Junction values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of Junction values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDIndexRaster write_Junction(LSDFlowInfo FlowInfo) { return write_integer_data_to_LSDIndexRaster(Junction, FlowInfo); }
  
  /// @brief Write NumberOfCells values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of NumberOfCells values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDIndexRaster write_NumberOfCells(LSDFlowInfo FlowInfo) { return write_integer_data_to_LSDIndexRaster(NumberOfCells, FlowInfo); }
  
  /// @brief Write BasinOrder values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDIndexRaster of BasinOrder values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDIndexRaster write_BasinOrder(LSDFlowInfo FlowInfo) { return write_integer_data_to_LSDIndexRaster(BasinOrder, FlowInfo); }
  
  /// @brief Write Area values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of Area values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_Area(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(Area, FlowInfo); }
  
  /// @brief Write SlopeMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of SlopeMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_SlopeMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(SlopeMean, FlowInfo); }
  
  /// @brief Write ElevationMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ElevationMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_ElevationMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ElevationMean, FlowInfo); }
  
  /// @brief Write AspectMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of AspectMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_AspectMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(AspectMean, FlowInfo); }
  
  /// @brief Write ReliefMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ReliefMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_ReliefMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ReliefMean, FlowInfo); }
  
  /// @brief Write PlanCurvMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of PlanCurvMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_PlanCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(PlanCurvMean, FlowInfo); }
  
  /// @brief Write ProfileCurvMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ProfileCurvMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_ProfileCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ProfileCurvMean, FlowInfo); }
  
  /// @brief Write TotalCurvMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of TotalCurvMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_TotalCurvMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(TotalCurvMean, FlowInfo); }
  
  /// @brief Write PlanCurvMax values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of PlanCurvMax values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_PlanCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(PlanCurvMax, FlowInfo); }
  
  /// @brief Write ProfileCurvMax values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of ProfileCurvMax values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_ProfileCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(ProfileCurvMax, FlowInfo); }
  
  /// @brief Write TotalCurvMax values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of TotalCurvMax values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_TotalCurvMax(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(TotalCurvMax, FlowInfo); }
  
  /// @brief Write HillslopeLength_HFR values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of HillslopeLength_HFR values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_HillslopeLength_HFR(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(HillslopeLength_HFR, FlowInfo); }
  
  /// @brief Write HillslopeLength_Binned values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of HillslopeLength_Binned values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_HillslopeLength_Binned(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(HillslopeLength_Binned, FlowInfo); }
  
  /// @brief Write HillslopeLength_Spline values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of HillslopeLength_Spline values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_HillslopeLength_Spline(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(HillslopeLength_Spline, FlowInfo); }
  
  /// @brief Write HillslopeLength_Density values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of HillslopeLength_Density values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_HillslopeLength_Density(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(HillslopeLength_Density, FlowInfo); }
  
  /// @brief Write FlowLength values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of FlowLength values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_FlowLength(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(FlowLength, FlowInfo); }
  
  /// @brief Write DrainageDensity values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of DrainageDensity values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_DrainageDensity(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(DrainageDensity, FlowInfo); }
  
  /// @brief Write CosmoErosionRate values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of CosmoErosionRate values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_CosmoErosionRate(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(CosmoErosionRate, FlowInfo); }
  
  /// @brief Write OtherErosionRate values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of OtherErosionRate values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_OtherErosionRate(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(OtherErosionRate, FlowInfo); }
  
  /// @brief Write CHTMean values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of CHTMean values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_CHTMean(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(CHTMean, FlowInfo); }
  
  /// @brief Write EStar values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of EStar values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_EStar(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(EStar, FlowInfo); }
  
  /// @brief Write RStar values into the shape of the basin.
  /// @param FlowInfo Flowinfo object.
  /// @return LSDRaster of RStar values in the shape of the basin.
  /// @author SWDG
  /// @date 12/12/13 
  LSDRaster write_RStar(LSDFlowInfo FlowInfo) { return write_real_data_to_LSDRaster(RStar, FlowInfo); }
   
  protected:
  
  //These instance variables are set at initialisation
  
  ///Number of rows.
  int NRows;
  ///Number of columns.
	int NCols;
	///Minimum X coordinate.
  float XMinimum;
	///Minimum Y coordinate.
	float YMinimum;
	///Data resolution.
	float DataResolution;
	///No data value.
	int NoDataValue;  
  /// Junction Number of the basin, serves as a unique ID of a basin.
  int Junction;
  ///Vector of all nodes in basin.
  vector<int> BasinNodes;
  /// Number of DEM cells. 
  int NumberOfCells;
  /// Area in spatial units of the basin.
  float Area;
  /// Boolean to show if the basin is beheaded or not.
  bool Beheaded;
  /// i index of the outlet pixel
  int Outlet_i;
  /// j index of the outlet pixel
  int Outlet_j;
  /// The Strahler order of the basin
  int BasinOrder;
  ///The i index of the centroid of the basin
  int Centroid_i;
  ///The j index of the centroid of the basin
  int Centroid_j;
  
  
  //These instance variables are used to store calculated basin parameters 
  
  /// Mean basin slope.
  float SlopeMean;
  /// Mean basin elevation.
  float ElevationMean;
  /// Mean basin aspect.
  float AspectMean;
  /// Mean basin relief.
  float ReliefMean;
  /// Mean basin planform curvature.
  float PlanCurvMean;
  /// Mean basin profile curvature.
  float ProfileCurvMean;
  /// Mean basin total curvature.
  float TotalCurvMean;
  /// Max basin planform curvature.
  float PlanCurvMax;
  /// Max basin profile curvature.
  float ProfileCurvMax;
  /// Max basin total curvature.
  float TotalCurvMax;
  /// Mean hillslope length from hilltop flow routing.
  float HillslopeLength_HFR;
  /// Mean hillslope length from binned boomerang plots.
  float HillslopeLength_Binned;
  /// Mean hillslope length from spline curves on boomerang plots.
  float HillslopeLength_Spline;
  /// Mean hillslope length from drainage density.
  float HillslopeLength_Density;
  /// Basin flowlength.
  float FlowLength;
  /// Basin drainage density.
  float DrainageDensity;
  /// Basin Perimeter's j index.
  vector<int> Perimeter_i;
  /// Basin Perimeter's j index.
  vector<int> Perimeter_j;
  /// Cosmo erosion rate.
  float CosmoErosionRate;
  /// Other erosion rate.
  float OtherErosionRate;
  /// Mean basin hilltop curvature.
  float CHTMean;
  /// Basin E* value.
  float EStar;
  /// Basin R* value.
  float RStar;
 											
  private:
	void create(int Junction, LSDFlowInfo& FlowInfo, LSDJunctionNetwork& ChanNet);

};

#endif
