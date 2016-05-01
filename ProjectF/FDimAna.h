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

    /// Bin widths in cm
    std::vector<double> binWidths;

    /// Vector of histograms at different scales
    std::vector<TH2D*> hitHistos;

    /// Profile of counting boxes
    TProfile* boxProf;

    /// Profile of counting boxes normalized with respect to the number of 1 mm x 1mm boxes
    TProfile* boxProfNorm;

    /// Profile of F dimension
    TProfile* dimProf;

    /// Make bin contents of histogram either 0 (if already 0) or 1 (if != 0)
    void binarize(TH2* h2);

    /// Tree to store histograms at different scales
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
