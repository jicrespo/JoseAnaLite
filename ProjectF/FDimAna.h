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

    /// Vector of histograms at different scales
    std::vector<TH2D*> hitHistos;

    /// Vector of counting boxes at different scales
    std::vector<double> boxNumbers;

    ///// F dimension from scaling fit
    //    double FDim;
    /// Vector with box counting F dimension from scaling fit, scale cutoff and fit chi square computed using several methods (see enum in FDimAna.cxx)
    std::vector<double> BDim;

    /// Vector with correlation F dimension from scaling fit, scale cutoff and fit chi square computed using several methods (see enum in FDimAna.cxx)
    std::vector<double> CDim;

    /// Vector with information F dimension from scaling fit, scale cutoff and fit chi square computed using several methods (see enum in FDimAna.cxx)
    std::vector<double> IDim;

    /// Vector with box counting F dimensions from local slopes
    std::vector<double> BDimSlope;

    /// Vector with correlation F dimensions from local slopes
    std::vector<double> CDimSlope;

    /// Vector with information F dimensions from local slopes
    std::vector<double> IDimSlope;

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

    /// Profile of counting boxes
    TProfile* boxProf;

    /// Profile of counting boxes normalized with respect to the number of the smallest boxes
    TProfile* boxProfNorm;

    /// Profile of F dimension
    TProfile* dimProf;

    /// 2D histogram of counting boxes
    TH2D* boxH2;

    /// 2D histogram of counting boxes normalized with respect to the number of the smallest boxes
    TH2D* boxH2Norm;

    /// 2D histogram of F dimension
    TH2D* dimH2;

    /// Profile of correlation integral
    TProfile* corProf;

    /// Profile of correlation integral removing normalization
    TProfile* corProfNorm;

    /// 2D histogram of correlation integral
    TH2D* corH2;

    /// 2D histogram of correlation integral removing normalization
    TH2D* corH2Norm;

    /// Profile of negative entropy
    TProfile* infProf;

    /// 2D histogram of negative entropy
    TH2D* infH2;

    /// Reset output variables to be filled in the tree
    void resetOutTreeVars();

    /// Make bin contents of histogram either 0 (if already 0) or 1 (if != 0)
    void binarize(TH2* h2);

    /// Return number of bins with content bigger than zero
    double countBoxes(TH2* h2);

    /// Return minus Shannon entropy in base 10
    double mentropy(TH2* h2);

    /// Power function to fit scaling
    static double scalingFitPow(double* x, double* par);

    /// Linear function to fit scaling
    static double scalingFitLin(double* x, double* par);

    /// Compute mean F dimension and a cutoff from local slopes and store the results in a vector passed by reference
    void meanFDim(const std::vector<double>& slopes, const std::vector<double>& widths, std::vector<double>& FDim );

    /// Compute F dimension and a cutoff from power scaling and store the results in a vector passed by reference
    void powFDim(TGraph*, const std::vector<double>& widths, std::vector<double>& FDim );

    /// Compute F dimension and a cutoff from linear scaling and store the results in a vector passed by reference
    void linFDim(const double& saturation, TGraph*, const std::vector<double>& widths, std::vector<double>& FDim );

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
