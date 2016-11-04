/**
 * \file FDimAna.h
 *
 * \ingroup ProjectF
 *
 * \brief Class def header for a class FDimAna
 *
 * @author jcrespo
 */

/** \addtogroup ProjectF

    @{*/

#ifndef LARLITE_FDIMANA_H
#define LARLITE_FDIMANA_H

#include "Analysis/ana_base.h"

#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>
#include <TGraph.h>

namespace larlite {
  /**
     \class FDimAna
     User custom analysis class made by SHELL_USER_NAME
   */
  class FDimAna : public ana_base{

  public:

    /// Default constructor
    FDimAna(){ _name="FDimAna"; _fout=0;}

    /// Default destructor
    virtual ~FDimAna(){}

    /** IMPLEMENT in FDimAna.cc!
	Initialization method to be called before the analysis event loop.
    */
    virtual bool initialize();

    /** IMPLEMENT in FDimAna.cc!
	Analyze a data event-by-event
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FDimAna.cc!
	Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
    ;

  private:

    /// Bin widths in m
    std::vector<double> binWidths;

    /// Distance widths in m
    std::vector<double> distWidths;

    /// Vector of histograms at different scales
    std::vector<TH2D*> hitHistos;

    /// Vector of counting boxes at different scales
    std::vector<double> boxNumbers;

    /// Vector with correlation F dimension from scaling fit and its error, scale transition and fit chi square
    std::vector<double> CDim;

    /// Vector with exponents for the generalized F dimension
    std::vector<double> QExp;

    /// Matrix with generalized F dimension from scaling fit and its error, scale cutoff and fit chi square
    std::vector< std::vector<double> > QDim;

    /// Vector with correlation F dimensions from local slopes
    std::vector<double> CDimSlope;

    /// Matrix with generalized F dimensions from local slopes
    std::vector< std::vector<double> > QDimSlope;

    /// MC truth shower energy
    double showerTruthE;

    /// MC truth track energy
    double trackTruthE;

    /// MC deposited shower energy
    double showerDepE;

    /// MC deposited track energy
    double trackDepE;

    /// Sum of hit integrals
    double hitSumInt;

    /// Vector with list of creation processes of MC tracks or showers different than primary
    std::vector<std::string> processTruth;

    /// Hit minimum charge
    double _chargeMin;

    /// Profile of generalized F dimension
    TProfile* QDimProf;

    /// Profile of correlation F dimension
    TProfile* CDimProf;

    /// Reset output variables to be filled in the tree
    void resetOutTreeVars();

    /// Return sum of bin probabilities raised to the power of q
    double probabilitySumQ(TH2* h2, double q);

    /// Linear function to fit scaling
    static double scalingFitLin(double* x, double* par);

    /// Logistic function to fit scaling
    static double scalingFitLogi(double* x, double* par);

    /// Compute mean F dimension and a cutoff from local slopes and store the results in a vector passed by reference
    void meanFDim(const std::vector<double>& slopes, const std::vector<double>& widths, std::vector<double>& FDim );

    /// Compute F dimension and a cutoff from linear scaling and store the results in a vector passed by reference
    void linFDim(const double& saturation, TGraph*, const std::vector<double>& widths, std::vector<double>& FDim );

    /// Compute F dimension and a scale limit and store the results in a vector passed by reference
    void logiFDim(TGraph*, const std::vector<double>& widths, std::vector<double>& FDim );

    /// Tree to store output variables related to FDimAna
    TTree* outTree;
  };
}
#endif

//**************************************************************************
//
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group
