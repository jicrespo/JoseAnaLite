/**
 * \file TPCWaveformAna.h
 *
 * \ingroup TPCWaveformAna
 * 
 * \brief Class def header for a class TPCWaveformAna
 *
 * @author jcrespo
 *
 * Based on the WFViewer class by David Caratelli in SNCompression/Compression
 */

/** \addtogroup TPCWaveformAna

    @{*/

#ifndef LARLITE_TPCWAVEFORMANA_H
#define LARLITE_TPCWAVEFORMANA_H

#include "Analysis/ana_base.h"
#include "LArUtil/GeometryHelper.h"
#include "DataFormat/rawdigit.h"
#include <TH2D.h>

namespace larlite {
  /**
     \class TPCWaveformAna
     Analysis class for basic analysis of TPC waveforms
   */
  class TPCWaveformAna : public ana_base{
  
  public:

    /// Default constructor
    TPCWaveformAna(){ _name="TPCWaveformAna"; _fout=0;}

    /// Default destructor
    virtual ~TPCWaveformAna(){}

    /** IMPLEMENT in TPCWaveformAna.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in TPCWaveformAna.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in TPCWaveformAna.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
    ;
    
  private:

    // Histogram showing the difference between the sample ADC values and the linear
    // interpolation using the previous and the following samples, as a function of channel number
    TH2D* _hDiffToInterpol;

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
