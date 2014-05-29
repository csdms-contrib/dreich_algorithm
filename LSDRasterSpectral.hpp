//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDRasterSpectral
// Land Surface Dynamics StatsTools
//
// An object for manipulating rasters developed for the University of Edinburgh
//  Land Surface Dynamics group topographic toolbox. This is a derivative class
// from LSDRaster, for use specifically with spectral analysis.
//
// These tools have been seperated from the LSDRaster class mainly because
//  they require the FFTW library and are therefore less portable than
//  the standard LSDRaster object.
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

/** @file LSDRasterSpectral.hpp
@author Simon M. Mudd, University of Edinburgh
@author David Milodowski, University of Edinburgh
@author Martin D. Hurst, British Geological Survey
@author Stuart W. D. Grieve, University of Edinburgh
@author Fiona Clubb, University of Edinburgh

@version Version 0.0.1
@brief This object performs spectral analysis.
@details It is seperate from LSDRaster simply because it requires the FFTW package so this can be removed from compilation to retain portability.

@date 02/04/2013
*/

//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
using namespace std;
using namespace TNT;

#ifndef LSDRasterSpectral_H
#define LSDRasterSpectral_H

/// @brief This object performs spectral analysis.
class LSDRasterSpectral: public LSDRaster
{
  public:
    //LSDRasterSpectral()		    	{ create(); }

    /// @brief Create an LSDRasterSpectral from a file.
    /// Uses a filename and file extension
    /// @return LSDRasterSpectral
    /// @param filename A String, the file to be loaded.
    /// @param extension A String, the file extension to be loaded.
    /// @author SMM
    /// @date 18/12/2012
    LSDRasterSpectral(string filename, string extension)	{ create(filename, extension); }
    
    /// @brief Create an LSDRasterSpectral from memory.
    /// @return LSDRasterSpectral
    /// @param nrows An integer of the number of rows.
    /// @param ncols An integer of the number of columns.
    /// @param xmin A float of the minimum X coordinate.
    /// @param ymin A float of the minimum Y coordinate.
    /// @param cellsize A float of the cellsize.
    /// @param ndv An integer of the no data value.
    /// @param data An Array2D of floats in the shape nrows*ncols,
    /// containing the data to be written.
    /// @author SMM
    /// @date 18/12/2012
    LSDRasterSpectral(int nrows, int ncols, float xmin, float ymin,
	          float cellsize, float ndv, Array2D<float> data)
			{ create(nrows, ncols, xmin, ymin, cellsize, ndv, data); }
	
    /// @brief Create an LSDRasterSpectral from an LSDRaster object.
    /// @param An_LSDRaster LSDRaster object.
    /// @return LSDRasterSpectral.
    /// @author SMM
    /// @date 18/12/2012
    LSDRasterSpectral(LSDRaster& An_LSDRaster)		{ create(An_LSDRaster); }

    /// @brief Create an LSDRasterSpectral object that has dimensions 2^raster_order.
    /// @param raster order, that is the order of the raster dimension where 
    /// the dimension is 2^raster_order.
    /// @param The size of the cells
    /// @param no data value
    /// @author SMM
    /// @date 18/02/2014
    LSDRasterSpectral(int raster_order, float cellsize, float ndv)
                     { create(raster_order, cellsize, ndv); }

    /// Assignment operator.
    LSDRasterSpectral& operator=(const LSDRasterSpectral& LSDR);

    // Fourier helper functions
    // these functions are used to manipulate fourier transformed data
    
    /// @brief This returns the frequency values of an UNSHIFTED DFT along the rows.
    /// @return A float vector containing the frequency 
    /// @author SMM
    /// @date 19/02/2014
    vector<float> get_row_direction_frequencies_unshifted(); 

    /// @brief This returns the frequency values of an UNSHIFTED DFT along the columns.
    /// @return A float vector containing the frequency 
    /// @author SMM
    /// @date 19/02/2014
    vector<float> get_col_direction_frequencies_unshifted(); 

    /// @brief This calucaltes a scaling array for scaling an unshifted DFT
    /// by the factor 1/f^beta.
    /// @param beta the fractal exponent
    /// @return a float array with the scaling factor 1/f^beta
    /// @author SMM
    /// @date 19/02/2014
    Array2D<float> get_frequency_scaling_array(float beta);

    /// @brief This creates a fractal surface using the spectral method.
    /// @details The method works as follows:\n
    ///  1) Generate a random surface.\n
    ///  2) Perform DFT on this random surface.\n
    ///  3) Scale the tranform (both real and imaginary parts) by 1/f^beta.\n
    ///  4) Perform the inverse DFT.\n 
    ///
    ///  This results in a pseudo fractal surface that can be used in comarison 
    ///  with real topography.
    /// @param beta value which is the scaling exponent
    /// @author SMM
    /// @date 20/02/2014
    void generate_fractal_surface_spectral_method(float beta);

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // FAST FOURIER TRANSFORM MODULE
    //------------------------------------------------------------------------------

    /// @brief Computes the forward fast fourier transform of a 2D discrete dataset.
    /// @param InputArray = zeta_padded (padded DEM).
    /// @param transform_direction = -1.
    /// @param OutputArrayReal = Real 2D spectrum.
    /// @param OutputArrayImaginary = Imaginary 2D spectrum.
    /// @author David Milodowski
    /// @date 18/12/2012
    void dfftw2D_fwd(Array2D<float>& InputArray, Array2D<float>& OutputArrayReal, Array2D<float>& OutputArrayImaginary,
	                 int transform_direction);

    /// @brief Computes the inverse fast fourier transform of a 2D discrete dataset.
    /// @param InputArrayReal = Real component of 2D spectrum.
    /// @param InputArrayImaginary = Imaginary component of 2D spectrum.
    /// @param OutputArray = reconstructed DEM.
    /// @param transform_direction = 1.
    /// @author David Milodowski
    /// @date 18/12/2012
    void dfftw2D_inv(Array2D<float>& InputArrayReal, Array2D<float>& InputArrayImaginary,
  	                 Array2D<float>& OutputArray, int transform_direction);

    /// @brief Detrend Data.
    ///
    /// @details Fit plane by least squares regression and use coefficients to determine local slope ax + by + c = z.
    /// @param zeta Input elevation data.
    /// @param zeta_detrend Output detrended elevation data.
    /// @param trend_plane Output array of trend plane.
    /// @author David Milodowski
    /// @date 18/12/2012
    void detrend2D(Array2D<float>& zeta, Array2D<float>& zeta_detrend, Array2D<float>& trend_plane);

    /// @brief Hann Window Module.
    ///
    /// @details Use 2D elliptical Hann (raised cosine) window on data matrix, to reduce spectral leakage and retain good frequency resolution.
    /// @param zeta_detrend Detrended elevation data
    /// @param zeta_Hann2D Output windowed data.
    /// @param Hann2D Output Hann window.
    /// @author David Milodowski
    /// @date 18/12/2012
    void window_data_Hann2D(Array2D<float>& zeta_detrend, Array2D<float>& zeta_Hann2D, Array2D<float>& Hann2D);

    /// @brief SHIFT ORIGIN OF SPECTRUM IN FOURIER DOMAIN.
    ///
    /// @details The output of the DFT algorithm must be rearranged to place the zero wavenumber element near the center of the array.
    /// @param spectrum_real
    /// @param spectrum_imaginary
    /// @param spectrum_real_shift
    /// @param spectrum_imaginary_shift
    /// @author David Milodowski
    /// @date 18/12/2012
    void shift_spectrum(Array2D<float>& spectrum_real,  Array2D<float>& spectrum_imaginary,
  	                    Array2D<float>& spectrum_real_shift, Array2D<float>& spectrum_imaginary_shift);

    /// @brief DE-SHIFT ORIGIN OF SPECTRUM.
    ///
    /// @details Inverse process of shift_spectrum() to return filtered spectrum to 
    /// original format required for the inverse fourier transform algorithm.
    /// @param FilteredSpectrumReal.
    /// @param FilteredSpectrumImaginary.
    /// @param FilteredSpectrumReal_deshift.
    /// @param FilteredSpectrumImaginary_deshift
    /// @author David Milodowski
    /// @date 18/12/2012
    void shift_spectrum_inv(Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
  	                        Array2D<float>& FilteredSpectrumReal_deshift, Array2D<float>& FilteredSpectrumImaginary_deshift);

    /// @brief CALCULATE THE DFT PERIODOGRAM.
    ///
    /// @details Multiply fourier analysis output by complex conjugate and normalises.
    /// @param spectrum_real_shift
    /// @param spectrum_imaginary_shift
    /// @author David Milodowski
    /// @date 18/12/2012
    void calculate_2D_PSD(Array2D<float>& spectrum_real_shift, Array2D<float>& spectrum_imaginary_shift);

    /// @brief GET RADIAL POWER SPECTRUM.
    ///
    /// @details Collapse 2D PSD into a radial PSD.
    /// @author David Milodowski
    /// @date 18/12/2012
    void calculate_radial_PSD();

    /// @brief COMPUTE DISCRETE FAST FOURIER TRANSFORM OF A REAL, 2-DIMENSIONAL DATASET.
    ///
    /// @details Computes the 2D and radial power spectra of a 2D array.
    /// @param file_id File identifier to prefix output files
    /// @param LogBinWidth Width of the logarithmically spaced bins. For topography, suggest this is 0.1 to start.
    /// @author David Milodowski
    /// @date 18/12/2012
    void fftw2D_spectral_analysis(char* file_id, float LogBinWidth);

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // FUNCTIONS TO ADD WEIGHTS TO FOURIER SPECTRA (FOR USE IN SPECTRA FILTERS)
    //------------------------------------------------------------------------------

    /// @brief BANDPASS FILTER.
    ///
    /// @details Filter array to band between frequency bands f1 and f2.  The bandpass filter
    /// is a gaussian filter centred at (f1+f2)/2 and with a SD of |f2-f1|/6.
    /// @param RawSpectrumReal
    /// @param RawSpectrumImaginary
    /// @param FilteredSpectrumReal
    /// @param f1
    /// @param f2
    /// @author David Milodowski
    /// @date 18/12/2012
    void bandpass_filter(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
  	                     Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
  	                     float f1, float f2);

    /// @brief LOWPASS FILTER.
    ///
    /// @details Filter array to retain frequencies below f1.  The filter edge is a radial gaussian function with a SD of |f2-f1|/3.
    /// @param RawSpectrumReal
    /// @param RawSpectrumImaginary
    /// @param FilteredSpectrumReal
    /// @param FilteredSpectrumImaginary
    /// @param f1
    /// @param f2
    /// @author David Milodowski
    /// @date 18/12/2012
    void lowpass_filter(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
  	                    Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
  	                    float f1, float f2);

    /// @brief HIGHPASS FILTER.
    ///
    /// @details Filter array to retain frequencies above f1.  The filter edge is a radial gaussian function with a SD of |f2-f1|/3.
    /// @param RawSpectrumReal
    /// @param RawSpectrumImaginary
    /// @param FilteredSpectrumReal
    /// @param FilteredSpectrumImaginary
    /// @param f1
    /// @param f2
    /// @author David Milodowski
    /// @date 18/12/2012
    void highpass_filter(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
  	                     Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary,
  	                     float f1, float f2);

    /// @brief WIENER FILTER.
    ///
    /// @details The Wiener filter is a spectral filter that removes noise from an image or DEM.
    ///  Weights in filter given by amplitudes of noise and signal: \n\n
    ///        phi(f) = |S(f)|^2/(|S(f)|^2 + |N(f)|^2)
    /// @param RawSpectrumReal
    /// @param RawSpectrumImaginary
    /// @param FilteredSpectrumReal
    /// @param FilteredSpectrumImaginary
    /// @param WSS Summed square of the weighting coefficients.
    /// @author David Milodowski
    /// @date 18/12/2012
    void wiener_filter(Array2D<float>& RawSpectrumReal, Array2D<float>& RawSpectrumImaginary,
  	                   Array2D<float>& FilteredSpectrumReal, Array2D<float>& FilteredSpectrumImaginary);

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // MAIN FUNCTIONS USING SPECTRAL FILTERS
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    /// @brief FAST FOURIER TRANSFORM FILTER FOR A REAL, 2-DIMENSIONAL DATASET.
    ///
    /// Note that FLow <= FHigh
    ///
    /// There are three types of filters depending on the intentions of the user
    ///
    /// BANDPASS FILTER (FilterType = 1) \n
    /// Filter array to band between frequency bands f1 and f2.  The bandpass filter
    /// is a gaussian filter centred at (f1+f2)/2 and with a SD of |f2-f1|/6.
    /// \n\n
    /// LOWPASS FILTER (FilterType = 2)\n
    /// Filter array to retain frequencies below f1.  The filter edge is a radial
    /// gaussian function with a SD of |f2-f1|/3.  f1 is the frequency below which
    /// the filter starts to taper; f2 is the frequency at which the filter tapers to
    /// zero. If f1 = f2, the edge is effectively a step function.
    /// \n\n
    /// HIGHPASS FILTER (FilterType = 3) \n
    /// Filter array to retain frequencies above f2.  The filter edge is a radial
    /// gaussian function with a SD of |f2-f1|/3.  f2 is the frequency below which
    /// the filter starts to taper; f1 is the frequency at which the filter tapers to
    /// zero. If f1 = f2, the edge is effectively a step function.
    /// \n\n
    /// A second type of bandpass filter is possible by combining the highpass and
    /// lowpass filters.
    ///
    /// @param FilterType
    /// @param FLow
    /// @param FHigh
    /// @author David Milodowski
    /// @date 18/12/2012
    LSDRaster fftw2D_filter(int FilterType, float FLow, float FHigh);

    /// @brief WIENER FILTER FOR A REAL, 2-DIMENSIONAL DATASET.
    ///
    /// The Wiener filter is a spectral filter that removes noise from an image or
    /// DEM.  Essentially, it works on the principle that the observed spectrum
    /// contains the superposition of the real signal and an additional noise signal,
    /// which we want to remove.  If we know, or can make a reasonable guess at the
    /// noise, N(f), and signal, S(f), parts of the spectrum then we can remove the
    /// noise using the filter:
    /// \n\n
    ///        phi(f) = |S(f)|^2/(|S(f)|^2 + |N(f)|^2)
    /// \n\n
    /// For topography; at long wavelengths the topographic signal obeys an
    /// approximate power law relationship between amplitude and frequency,
    /// decreasing as the frequency increases (and wavelength decreases).  Noise
    /// typically dominates the high frequency part of the spectrum.  Thus at high
    /// frequencies the spectrum is dominated by noise, and the filter weight goes to
    /// zero.  In contrast, at low frequencies, the signal dominates and the filter
    /// weight goes to 1.
    /// \n\n
    /// The optimal wiener filter is described in more detail in Numerical Recipes,
    /// 13.3, p149.
    /// \n\n
    /// The exact structure of the noise is worth thinking about.  White noise, which
    /// is random, has equal power across all wavelengths.  In the instance of
    /// topography, noise can be created by a whole range of sources, from rock
    /// exposure, to pit and mound topography, to unfiltered vegetation etc.  It is
    /// likely that these sources will not produce purely white noise, but rather
    /// will show an element of structure.  This program makes two assumptions about
    /// the noise: i) it dominates the signal at high frequencies (close to the
    /// Nquist frequency) and ii) we can reasonably model this using a linear fit in
    /// log-log space - i.e. it obeys some form of power law function between
    /// frequency and amplitude.  Note that if the noise in the signal is really
    /// white noise, then the power law function for the noise would simply have an
    /// exponent of zero.  I prefer this formulation because it permits the
    /// characterisation of the noise model without assuming that the noise has a
    /// particular structure (white noise, pink noise etc.)
    /// @author David Milodowski
    /// @date 18/12/2012
    LSDRaster fftw2D_wiener();

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // FUNCTIONS TO PRINT RADIAL SPECTRA
    //------------------------------------------------------------------------------
    void print_radial_spectrum(float bin_width, string file_id);

    
  protected:
    int Lx;
    int Ly;
    float WSS;
    Array2D<float> P_DFT;
    vector<float> RadialFrequency;
    vector<float> RadiallyAveragedPSD;
	
  private:
    void create();
    void create(string filename, string extension);
    void create(int ncols, int nrows, float xmin, float ymin,
                float cellsize, float ndv, Array2D<float> data);
    void create(int raster_order, float cellsize, float ndv);
    void create(LSDRaster& An_LSDRaster);

};

#endif
