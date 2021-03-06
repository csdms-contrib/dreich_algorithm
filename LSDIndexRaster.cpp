//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDIndexRaster
// Land Surface Dynamics IndexRaster
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for manipulating
//  and analysing raster data, with a particular focus on topography
//
// The IndexRaster object stores only integer values and is used mostly
//  for storing indices into raster data.
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

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDIndexRaster.cpp
// cpp file for the LSDIndexRaster object
// LSD stands for Land Surface Dynamics
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This object is written by
// Simon M. Mudd, University of Edinburgh
// David Milodowski, University of Edinburgh
// Martin D. Hurst, British Geological Survey
// <your name here>
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Version 0.0.1		20/08/2012
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "TNT/tnt.h"
#include "LSDIndexRaster.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDIndexRaster_CPP
#define LSDIndexRaster_CPP

// operators
// SMM 2012
LSDIndexRaster& LSDIndexRaster::operator=(const LSDIndexRaster& rhs)
 {
  if (&rhs != this)
   {
    create(rhs.get_NRows(),rhs.get_NCols(),rhs.get_XMinimum(),rhs.get_YMinimum(),
           rhs.get_DataResolution(), rhs.get_NoDataValue(), rhs.get_RasterData() );
   }
  return *this;
 }

// the create function. This is default and throws an error
// SMM 2012
void LSDIndexRaster::create()
{
	cout << "You need to initialize with a filename!" << endl;
	exit(EXIT_FAILURE);
}

// this creates a raster using an infile
// SMM 2012
void LSDIndexRaster::create(string filename, string extension)
{
	read_raster(filename,extension);
}

// this creates a raster filled with no data values
// SMM 2012
void LSDIndexRaster::create(int nrows, int ncols, float xmin, float ymin,
            float cellsize, int ndv, Array2D<int> data)
{
	NRows = nrows;
	NCols = ncols;
	XMinimum = xmin;
	YMinimum = ymin;
	DataResolution = cellsize;
	NoDataValue = ndv;

	RasterData = data.copy();

	if (RasterData.dim1() != NRows)
	{
		cout << "dimesntion of data is not the same as stated in NRows!" << endl;
		exit(EXIT_FAILURE);
	}
	if (RasterData.dim2() != NCols)
	{
		cout << "dimesntion of data is not the same as stated in NRows!" << endl;
		exit(EXIT_FAILURE);
	}

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function reads a DEM
// One has to provide both the filename and the extension
// the '.' between the filename and extension is not included
// for example, if the full filename is test.asc
// then
// filename = "test"
// and
// ext = "asc"
// The full filename coult also be "test.01.asc"
// so filename would be "test.01"
// and ext would again be "asc"
// SMM 2012
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexRaster::read_raster(string filename, string extension)
{
	string string_filename;
	string dot = ".";
	string_filename = filename+dot+extension;
	cout << "The filename is " << string_filename << endl;


	if (extension == "asc")
	{
		// open the data file
		ifstream data_in(string_filename.c_str());

		//Read in raster data
		string str;			// a temporary string for discarding text

		// read the georeferencing data and metadata
		data_in >> str >> NCols >> str >> NRows
			    >> str >> XMinimum >> str >> YMinimum
		   		>> str >> DataResolution
			    >> str >> NoDataValue;

		cout << "Loading asc file; NCols: " << NCols << " NRows: " << NRows << endl
		     << "X minimum: " << XMinimum << " YMinimum: " << YMinimum << endl
		     << "Data Resolution: " << DataResolution << " and No Data Value: "
		     << NoDataValue << endl;

		// this is the array into which data is fed
		Array2D<int> data(NRows,NCols,NoDataValue);

		// read the data
		for (int i=0; i<NRows; ++i)
		{
			for (int j=0; j<NCols; ++j)
			{
				data_in >> data[i][j];
			}
		}
		data_in.close();

		// now update the objects raster data
		RasterData = data.copy();
	}
	else if (extension == "flt")
	{
		// float data (a binary format created by ArcMap) has a header file
		// this file must be opened first
		string header_filename;
		string header_extension = "hdr";
		header_filename = filename+dot+header_extension;

		ifstream ifs(header_filename.c_str());
		if( ifs.fail() )
		{
			cout << "\nFATAL ERROR: the header file \"" << header_filename
				 << "\" doesn't exist" << std::endl;
			exit(EXIT_FAILURE);
		}
		else
		{
			string str;
			ifs >> str >> NCols >> str >> NRows
				>> str >> XMinimum >> str >> YMinimum
				>> str >> DataResolution
				>> str >> NoDataValue;
		}
		ifs.close();

		cout << "Loading asc file; NCols: " << NCols << " NRows: " << NRows << endl
			 << "X minimum: " << XMinimum << " YMinimum: " << YMinimum << endl
		     << "Data Resolution: " << DataResolution << " and No Data Value: "
		     << NoDataValue << endl;

		// this is the array into which data is fed
		Array2D<int> data(NRows,NCols,NoDataValue);

		// now read the DEM, using the binary stream option
		ifstream ifs_data(string_filename.c_str(), ios::in | ios::binary);
		if( ifs_data.fail() )
		{
			cout << "\nFATAL ERROR: the data file \"" << string_filename
			     << "\" doesn't exist" << endl;
			exit(EXIT_FAILURE);
		}
		else
		{
			float temp;
			for (int i=0; i<NRows; ++i)
			{
				for (int j=0; j<NCols; ++j)
				{
					ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
					data[i][j] = int(temp);
				}
			}
		}
		ifs_data.close();

		// now update the objects raster data
		RasterData = data.copy();
	}
	else
	{
		cout << "You did not enter and approprate extension!" << endl
				  << "You entered: " << extension << " options are .flt and .asc" << endl;
		exit(EXIT_FAILURE);
	}


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// write_raster
// this function writes a raster. One has to give the filename and extension
// currently the options are for .asc and .flt files
// SMM 2012
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void LSDIndexRaster::write_raster(string filename, string extension)
{
	string string_filename;
	string dot = ".";
	string_filename = filename+dot+extension;
	cout << "The filename is " << string_filename << endl;

	// this first bit of logic is for the asc file.
	if (extension == "asc")
	{
		// open the data file
		ofstream data_out(string_filename.c_str());

		if( data_out.fail() )
		{
			cout << "\nFATAL ERROR: unable to write to " << string_filename << endl;
			exit(EXIT_FAILURE);
		}

		data_out <<  "ncols         " << NCols
				<< "\nnrows         " << NRows
				<< "\nxllcorner     " << setprecision(14) << XMinimum
				<< "\nyllcorner     " << setprecision(14) << YMinimum
				<< "\ncellsize      " << DataResolution
				<< "\nNODATA_value  " << NoDataValue << endl;

		for (int i=0; i<NRows; ++i)
		{
			for (int j=0; j<NCols; ++j)
			{
				data_out << setprecision(6) << RasterData[i][j] << " ";
			}
			if (i != NRows-1) data_out << endl;
		}
		data_out.close();

	}
	else if (extension == "flt")
	{
		// float data (a binary format created by ArcMap) has a header file
		// this file must be opened first
		string header_filename;
		string header_extension = "hdr";
		header_filename = filename+dot+header_extension;

		ofstream header_ofs(header_filename.c_str());
		string str;
		header_ofs <<  "ncols         " << NCols
			<< "\nnrows         " << NRows
			<< "\nxllcorner     " << setprecision(14) << XMinimum
			<< "\nyllcorner     " << setprecision(14) << YMinimum
			<< "\ncellsize      " << DataResolution
			<< "\nNODATA_value  " << NoDataValue
			<< "\nbyteorder     LSBFIRST" << endl;
		header_ofs.close();

		// now do the main data
		ofstream data_ofs(string_filename.c_str(), ios::out | ios::binary);
		float temp;
		for (int i=0; i<NRows; ++i)
		{
			for (int j=0; j<NCols; ++j)
			{
				temp = float(RasterData[i][j]);
				data_ofs.write(reinterpret_cast<char *>(&temp),sizeof(temp));
			}
		}
		data_ofs.close();
	}
	else
	{
		cout << "You did not enter and approprate extension!" << endl
				  << "You entered: " << extension << " options are .flt and .asc" << endl;
		exit(EXIT_FAILURE);
	}

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Calculate the minimum bounding rectangle for an LSDIndexRaster Object and crop out
// all the surrounding NoDataValues to reduce the size and load times of output
// rasters.
//
// Ideal for use with chi analysis tools which output basin and chi m value rasters
// which can be predominantly no data. As an example, a 253 Mb file can be reduced to
// ~5 Mb with no loss or resampling of data.
//
// Returns A trimmed LSDIndexRaster object.
//
// SWDG 22/08/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDIndexRaster::RasterTrimmer(){

  //minimum index value in a column
  int a = 0;
  int min_col = 100000; //a big number

  for (int row = 0; row < NRows; ++row){
    a = 0;
    while (RasterData[row][a] == NoDataValue && a < NCols-1){
      ++a;
    }
    if (min_col > a){
      min_col = a;
    }
  }

  //maximum index value in a column
  a = NCols - 1;
  int max_col = 0; //a small number

  for (int row = 0; row < NRows; ++row){
    a = NCols - 1;
    while (RasterData[row][a] == NoDataValue && a > 0){
      --a;
    }
    if (max_col < a){
      max_col = a;
    }
  }

  //minimum index value in a row
  a = 0;
  int min_row = 100000; //a big number

  for (int col = 0; col < NCols; ++col){
    a = 0;
    while (RasterData[a][col] == NoDataValue && a < NRows - 1){
      ++a;
    }
    if (min_row > a){
      min_row = a;
    }
  }

  //maximum index value in a row
  a = NRows - 1;
  int max_row = 0; //a small number

  for (int col = 0; col < NCols; ++col){
    a = NRows - 1;
    while (RasterData[a][col] == NoDataValue && a > 0){
      --a;
    }
    if (max_row < a){
      max_row = a;
    }
  }

  // create new row and col sizes taking account of zero indexing
  int new_row_dimension = (max_row-min_row) + 1;
  int new_col_dimension = (max_col-min_col) + 1;

  Array2D<int>TrimmedData(new_row_dimension, new_col_dimension, NoDataValue);

  //loop over min bounding rectangle and store it in new array of shape new_row_dimension x new_col_dimension
  int TrimmedRow = 0;
  int TrimmedCol = 0;
  for (int row = min_row - 1; row < max_row; ++row){
    for(int col = min_col - 1; col < max_col; ++col){
      TrimmedData[TrimmedRow][TrimmedCol] = RasterData[row][col];
      ++TrimmedCol;
    }
    ++TrimmedRow;
    TrimmedCol = 0;
  }

  //calculate lower left corner coordinates of new array
  float new_XLL = ((min_col - 1) * DataResolution) + XMinimum;
  float new_YLL = YMinimum + ((NRows - (max_row + 0)) * DataResolution);

  LSDIndexRaster TrimmedRaster(new_row_dimension, new_col_dimension, new_XLL,
                          new_YLL, DataResolution, NoDataValue, TrimmedData);

  return TrimmedRaster;

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Make LSDIndexRaster object using a 'template' raster and an Array2D of data.
// SWDG 2/9/13
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LSDIndexRaster LSDIndexRaster::LSDRasterTemplate(Array2D<int> InputData){

  //do a dimensions check and exit on failure
  if (InputData.dim1() == NRows && InputData.dim2() == NCols){
    LSDIndexRaster OutputRaster(NRows, NCols, XMinimum, YMinimum, DataResolution, NoDataValue, InputData);
    return OutputRaster;
  }
  else{
   	cout << "Array dimensions do not match template LSDIndexRaster object" << endl;
		exit(EXIT_FAILURE);
  }

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function implements the thinning algorithm described in Rosenfeld and
// Kak (1982).  It takes a binary map and turns it into a continuous single 
// thread skeleton.  At present, pixels at the limits of the raster are
// automatically set to 0.  If it is necessary to extend the skeleton to the 
// edge, this should be a straightforward operation afterwards.
//
// Added by DTM 28/10/2013
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LSDIndexRaster LSDIndexRaster::thin_to_single_thread_network()
{
  Array2D<int> Skeleton(NRows,NCols,NoDataValue);
  Array2D<int> SkeletonUpdate(NRows,NCols,NoDataValue);
  // Step 1:- outlne borders
  // Corners of array
  // SE
  if(RasterData[0][0]!=NoDataValue)
  { 
    if(RasterData[0][0]==0) Skeleton[0][0]=0;
    else if((RasterData[1][0]==0) && (RasterData[0][1]==0)) Skeleton[0][0]=1;
    else if((RasterData[1][0]==0) || (RasterData[0][1]==0)) Skeleton[0][0]=2;
    else Skeleton[0][0]=3;
  }
  // SW
  if(RasterData[0][NCols-1]!=NoDataValue)
  { 
    if(RasterData[0][NCols-1]==0) Skeleton[0][NCols-1]=0;
    else if((RasterData[1][NCols-1]==0) && (RasterData[0][NCols-2]==0)) Skeleton[0][NCols-1]=1;
    else if((RasterData[1][NCols-1]==0) || (RasterData[0][NCols-2]==0)) Skeleton[0][NCols-1]=2;
    else Skeleton[0][0]=3;
  }
  // NE
  if(RasterData[NRows-1][0]!=NoDataValue)
  { 
    if(RasterData[NRows-1][0]==0) Skeleton[NRows-1][0]=0;
    else if((RasterData[NRows-2][0]==0) && (RasterData[NRows-1][1]==0)) Skeleton[NRows-1][0]=1;
    else if((RasterData[NRows-2][0]==0) || (RasterData[NRows-1][1]==0)) Skeleton[NRows-1][0]=2;
    else Skeleton[NRows-1][0]=3;
  }
  // SE
  if(RasterData[NRows-1][NCols-1]!=NoDataValue)
  { 
    if(RasterData[NRows-1][NCols-1]==0) Skeleton[NRows-1][NCols-1]=0;
    else if((RasterData[NRows-2][NCols-1]==0) && (RasterData[NRows-1][NCols-2]==0)) Skeleton[NRows-1][NCols-1]=1;
    else if((RasterData[NRows-2][NCols-1]==0) || (RasterData[NRows-1][NCols-2]==0)) Skeleton[NRows-1][NCols-1]=2;
    else Skeleton[NRows-1][NCols-1]=3;
  }
  // Edges of Array
  for (int i = 1; i < NRows-1; ++i)
  {
    //North
    if(RasterData[i][0]!=NoDataValue)
    {
      if(RasterData[i][0]==0) Skeleton[i][0]=0;
      else if((RasterData[i+1][0]==0) && (RasterData[i-1][0]==0)) Skeleton[i][0]=1;
      else if((RasterData[i+1][0]==0) || (RasterData[i-1][0]==0) || (RasterData[i][1]==0)) Skeleton[i][0]=2;
      else Skeleton[i][0]=3;
    }
    // South
    if(RasterData[i][NCols-1]!=NoDataValue)
    { 
      if(RasterData[i][NCols-1]==0) Skeleton[i][NCols-1]=0;
      else if((RasterData[i+1][NCols-1]==0) && (RasterData[i-1][NCols-1]==0)) Skeleton[i][NCols-1]=1;
      else if((RasterData[i+1][NCols-1]==0) || (RasterData[i-1][NCols-1]==0) || (RasterData[i][NCols-2]==0)) Skeleton[i][NCols-1]=2;
      else Skeleton[i][NCols-1]=3;
    }
  }
  for (int j = 1; j < NRows-1; ++j)
  {  
    // East
    if(RasterData[0][j]!=NoDataValue)
    {
      if(RasterData[0][j]==0) Skeleton[0][j]=0;
      else if((RasterData[0][j+1]==0) && (RasterData[0][j-1]==0)) Skeleton[0][j]=1;
      else if((RasterData[0][j+1]==0) || (RasterData[0][j-1]==0) || (RasterData[1][j]==0)) Skeleton[0][j]=2;
      else Skeleton[0][j]=3;
    }  
    // West
    if(RasterData[NRows-1][j]!=NoDataValue)
    {
      if(RasterData[NRows-1][j]==0) Skeleton[NRows-1][j]=0;
      else if((RasterData[NRows-1][j+1]==0) && (RasterData[NRows-1][j-1]==0)) Skeleton[NRows-1][j]=1;
      else if((RasterData[NRows-1][j+1]==0) || (RasterData[NRows-1][j-1]==0) || (RasterData[NRows-2][j]==0)) Skeleton[NRows-1][j]=2;
      else Skeleton[NRows-1][j]=3;
    }
  }
  // for rest of array
  for (int i = 1; i < NRows-1; ++i)
  {
    for (int j = 1; j< NCols-1; ++j)
    {
      if(RasterData[i][j] == 1)
      {
        // Skeleton
        if ( ((RasterData[i+1][j] == 0) && (RasterData[i-1][j] == 0)) || (RasterData[i][j+1] == 0 && (RasterData[i][j-1] == 0)) ) Skeleton[i][j] = 1;
        // Non skeleton boundaries
        else if(RasterData[i-1][j] == 0 || RasterData[i+1][j] == 0 || RasterData[i][j-1] == 0 || RasterData[i][j-1] == 0) Skeleton = 2;
        // Non boundary
        else Skeleton = 3;
      }
      else Skeleton = 0;
    }
  }
  //----------------------------------------------------------------------------
  // Step 2:- Now loop through the array again, progressively searching for 
  // N,E,S,W boundaries
  bool IsSkeleton = false;
  SkeletonUpdate = Skeleton.copy();
  while (IsSkeleton == false)
  {
    IsSkeleton = true;
    // NORTH BOUNDARIES
    // For the edges - E 
    for (int i = 1; i < NRows-1; ++i) // Northern most row can't have a Northern boundary and it makes little sense to do the Southern most row as there is no row beneath!
    { 
      int j=0;
      if((Skeleton[i][j] == 2 && Skeleton[i-1][j] == 0))
      {
        IsSkeleton = false;
//         // Logic for case where you are at an end, but the feature is two 
//         // pixels wide - do not want to truncate skeleton here!
//         if((Skeleton[i][j+1]==0) && (Skeleton[i+1][j+1]==0) && (Skeleton[i+1][j]==2)) SkeletonUpdate[i][j] = 2;
//         else
//         {
          SkeletonUpdate[i][j] = 0;
          // Check to see whether pixel to South is a border - if yes, make it
          // a skeleton pixel
          if (Skeleton[i+1][j] == 2) SkeletonUpdate[i+1][j] = 1;
          // Otherwise it will be the boundary pixel for the next round
          else SkeletonUpdate[i+1][j] = 2;
//         }
      }
      // W
      j=NCols-1;
      if((Skeleton[i][j] == 2 && Skeleton[i-1][j] == 0))
      {
        IsSkeleton = false;
//         if((Skeleton[i][j-1]==0) && (Skeleton[i+1][j-1]==0) && (Skeleton[i+1][j]==2)) SkeletonUpdate[i][j] = 2; 
//         else
//         {
          SkeletonUpdate[i][j] = 0;
          if (Skeleton[i+1][j] == 2) SkeletonUpdate[i+1][j] = 1;
          else SkeletonUpdate[i+1][j] = 2;
//         }
      }      
    }

    // for the rest of the array    
    for (int i = 1; i < NRows-1; ++i)
    {
      for (int j = 1; j< NCols-1; ++j)
      {            
        if((Skeleton[i][j] == 2) && (Skeleton[i-1][j] == 0))
        {
          IsSkeleton = false;
          // Logic for case where you are at an end, but the feature is two 
          // pixels wide - do not want to truncate skeleton here!
          if( ((Skeleton[i][j+1]==0) && (Skeleton[i+1][j+1]==0) && (Skeleton[i+1][j]==2))
           || ((Skeleton[i][j-1]==0) && (Skeleton[i+1][j-1]==0) && (Skeleton[i+1][j]==2)) )
          {
            SkeletonUpdate[i][j] = 2;
          }  
          else
          {
            SkeletonUpdate[i][j] = 0;
            // Check to see whether pixel to South is a border - if yes, make it
            // a skeleton pixel
            if (Skeleton[i+1][j] == 2) SkeletonUpdate[i+1][j] = 1;
            // Otherwise it will be the boundary pixel for the next round
            else SkeletonUpdate[i+1][j] = 2;
          }
        }
      }
    }
    Skeleton = SkeletonUpdate.copy();
    
    // SOUTH BOUNDARIES
    // For the edges - E and W
    for (int i = 1; i < NRows - 1; ++i)
    { 
      int j=0;
      if((Skeleton[i][j] == 2) && (Skeleton[i+1][j] == 0))
      {
        IsSkeleton = false;
        SkeletonUpdate[i][j] = 0;
        // Check to see whether pixel to South is a border - if yes, make it
        // a skeleton pixel
        if (Skeleton[i-1][j] == 2) SkeletonUpdate[i-1][j] = 1;
        // Otherwise it will be the boundary pixel for the next round
        else SkeletonUpdate[i-1][j] = 2;
      }
      j=NCols-1;
      if((Skeleton[i][j] == 2 && Skeleton[i+1][j] == 0))
      {
        IsSkeleton = false;
        SkeletonUpdate[i][j] = 0;
        if (Skeleton[i-1][j] == 2) SkeletonUpdate[i-1][j] = 1;
        else SkeletonUpdate[i-1][j] = 2;
      }      
    }

    // for the rest of the array    
    for (int i = 1; i < NRows-1; ++i)
    {
      for (int j = 1; j< NCols-1; ++j)
      {            
        if((Skeleton[i][j] == 2 && Skeleton[i+1][j] == 0))
        {
          IsSkeleton = false;
          // Logic for case where you are at an end, but the feature is two 
          // pixels wide - do not want to truncate skeleton here!
          if( ((Skeleton[i][j+1]==0) && (Skeleton[i-1][j+1]==0) && (Skeleton[i-1][j]==2))
           || ((Skeleton[i][j-1]==0) && (Skeleton[i-1][j-1]==0) && (Skeleton[i-1][j]==2)) )
          {
            SkeletonUpdate[i][j] = 2;
          }  
          else
          {
            SkeletonUpdate[i][j] = 0;
            // Check to see whether pixel to North is a border - if yes, make it
            // a skeleton pixel
            if (Skeleton[i-1][j] == 2) SkeletonUpdate[i-1][j] = 1;
            // Otherwise it will be the boundary pixel for the next round
            else SkeletonUpdate[i-1][j] = 2;
          }
        }
      }
    }
    Skeleton = SkeletonUpdate.copy();
    // EAST-FACING BOUNDARIES
    // For the edges - N and S
    for (int j = 1; j < NCols; ++j)
    { 
      // South edge
      int i=0;
      if((Skeleton[i][j] == 2 && Skeleton[i][j-1] == 0))
      {
        IsSkeleton = false;
        // Logic for case where you are at an end, but the feature is two 
        // pixels wide - do not want to truncate skeleton here!
        SkeletonUpdate[i][j] = 0;
        // Check to see whether pixel to West is a border - if yes, make it
        // a skeleton pixel
        if (Skeleton[i][j+1] == 2) SkeletonUpdate[i][j+1] = 1;
        // Otherwise it will be the boundary pixel for the next round
        else SkeletonUpdate[i][j+1] = 2;
      }
      // North edge
      i=NRows-1;
      if((Skeleton[i][j] == 2 && Skeleton[i][j-1] == 0))
      {
        IsSkeleton = false;
        SkeletonUpdate[i][j] = 0;
        if (Skeleton[i][j+1] == 2) SkeletonUpdate[i][j+1] = 1;
        else SkeletonUpdate[i][j+1] = 2;
      }      
    }
    
    // for the rest of the array    
    for (int i = 1; i < NRows-1; ++i)
    {
      for (int j = 1; j< NCols-1; ++j)
      {
        // For rest of the array              
        if((Skeleton[i][j] == 2 && Skeleton[i][j-1] == 0))
        {
          IsSkeleton = false;
          // Logic for case where you are at an end, but the feature is two 
          // pixels wide - do not want to truncate skeleton here!
          if( ((Skeleton[i+1][j]==0) && (Skeleton[i+1][j+1]==0) && (Skeleton[i][j+1]==2))
           || ((Skeleton[i-1][j]==0) && (Skeleton[i-1][j+1]==0) && (Skeleton[i][j+1]==2)) )
          {
            SkeletonUpdate[i][j] = 2;
          }  
          else
          {
            SkeletonUpdate[i][j] = 0;
            // Check to see whether pixel to West is a border - if yes, make it
            // a skeleton pixel
            if (Skeleton[i][j+1] == 2) SkeletonUpdate[i][j+1] = 1;
            // Otherwise it will be the boundary pixel for the next round
            else SkeletonUpdate[i][j+1] = 2;
          }
        }
      }
    }
    Skeleton = SkeletonUpdate.copy();
    // WEST-FACING BOUNDARIES
    // For the edges - N and S
    for (int j = 1; j < NCols; ++j)
    { 
      // South edge
      int i=0;
      if((Skeleton[i][j] == 2 && Skeleton[i][j+1] == 0))
      {
        IsSkeleton = false;
        // Logic for case where you are at an end, but the feature is two 
        // pixels wide - do not want to truncate skeleton here!
        SkeletonUpdate[i][j] = 0;
        // Check to see whether pixel to West is a border - if yes, make it
        // a skeleton pixel
        if (Skeleton[i][j-1] == 2) SkeletonUpdate[i][j-1] = 1;
        // Otherwise it will be the boundary pixel for the next round
        else SkeletonUpdate[i][j-1] = 2;
      }
      // North edge
      i=NRows-1;
      if((Skeleton[i][j] == 2 && Skeleton[i][j+1] == 0))
      {
        IsSkeleton = false;
        SkeletonUpdate[i][j] = 0;
        if (Skeleton[i][j-1] == 2) SkeletonUpdate[i][j-1] = 1;
        else SkeletonUpdate[i][j-1] = 2;
      }      
    }
    // For rest of the array 
    for (int i=1; i<NRows-1; ++i)
    {
      for (int j=1; j<NCols-1; ++j)
      {                     
        if((Skeleton[i][j] == 2 && Skeleton[i][j+1] == 0))
        {
          IsSkeleton = false;
          // Logic for case where you are at an end, but the feature is two 
          // pixels wide - do not want to truncate skeleton here!
          if( ((Skeleton[i+1][j]==0) && (Skeleton[i+1][j-1]==0) && (Skeleton[i][j-1]==2))
           || ((Skeleton[i-1][j]==0) && (Skeleton[i-1][j-1]==0) && (Skeleton[i][j-1]==2)) )
          {
            SkeletonUpdate[i][j] = 2;
          }  
          else
          {
            SkeletonUpdate[i][j] = 0;
            // Check to see whether pixel to West is a border - if yes, make it
            // a skeleton pixel
            if (Skeleton[i][j-1]==2) SkeletonUpdate[i][j-1] = 1;
            // Otherwise it will be the boundary pixel for the next round
            else SkeletonUpdate[i][j-1] = 2;
          }
        }
      }
    }
    Skeleton = SkeletonUpdate.copy();    
  }
  
  // finally, loop through, and remove any remaining 3-pixels, which should only
  // be skeleton pixels that are in this arrangement:
  //                            0 1 0
  //                            1 3 1
  //                            0 1 0
  for (int i = 0; i < NRows; ++i)
  {
    for (int j = 0; j< NCols; ++j)
    {
      if (Skeleton[i][j] == 3) Skeleton[i][j] = 1;
    }
  }
  
  LSDIndexRaster skeleton_raster(NRows,NCols,XMinimum,YMinimum,DataResolution,NoDataValue,Skeleton);
	return skeleton_raster;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Method to resample an LSDIndexRaster to a lower resolution. 
// OutputResolution is the resolution in spatial units to be resampled to.
// Returns an LSDRaster resampled to the OutputResolution.
// SWDG 17/3/14
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=   
LSDIndexRaster LSDIndexRaster::Resample(float OutputResolution){

  if (OutputResolution < DataResolution){
    cout << "Your resample resolution of " << OutputResolution << " is lower that the current data resolution " << DataResolution << endl;
	  exit(EXIT_FAILURE);
  } 

  int NewNRows = (NRows*DataResolution/OutputResolution); 
  int NewNCols = (NCols*DataResolution/OutputResolution);

  Array2D<int> Resampled(NewNRows, NewNCols, NoDataValue);
  
  int centre_i;
  int centre_j;   
  
  float ResolutionRatio = OutputResolution/DataResolution;
  
  for (int i = 0; i < NewNRows; ++i){
    for (int j = 0; j < NewNCols; ++j){
    
      //find the centre of the new grid in the old grid units
      centre_i = (i*ResolutionRatio) + (ResolutionRatio/2);
      centre_j = (j*ResolutionRatio) + (ResolutionRatio/2);
     
      Resampled[i][j] = RasterData[centre_i][centre_j];           
      
    }
  }                              

  LSDIndexRaster OutputRaster(NewNRows,NewNCols,XMinimum,YMinimum,OutputResolution,NoDataValue,Resampled);
  return OutputRaster;

}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#endif
