#ifndef LARLITE_TPCWAVEFORMANA_CXX
#define LARLITE_TPCWAVEFORMANA_CXX

#include "TPCWaveformAna.h"

namespace larlite {

  bool TPCWaveformAna::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    _hDiffToInterpol = new TH2D("hDiffToInterpol", 
				"Difference to interpolation; Channel; ADC_{i} - (ADC_{i+1} + ADC_{i-1})/2 (ADC)", 
				8256, 0, 8256, 4096, 0, 4096);

    return true;
  }
  
  bool TPCWaveformAna::analyze(storage_manager* storage) {
  
    //
    // Do your event-by-event analysis here. This function is called for 
    // each event in the loop. You have "storage" pointer which contains 
    // event-wise data. To see what is available, check the "Manual.pdf":
    //
    // http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
    // 
    // Or you can refer to Base/DataFormatConstants.hh for available data type
    // enum values. Here is one example of getting PMT waveform collection.
    //
    // event_fifo* my_pmtfifo_v = (event_fifo*)(storage->get_data(DATA::PMFIFO));
    //
    // if( event_fifo )
    //
    //   std::cout << "Event ID: " << my_pmtfifo_v->event_id() << std::endl;
    //
    
    // Get pointer to vector with waveforms
    event_rawdigit* ev_wf = storage->get_data<event_rawdigit>("daq");

    if( !ev_wf ){
      print(msg::kERROR, __FUNCTION__ , "Data storage did not find associated waveforms!");
      return false;
    }

    // Loop over all waveforms in event
    for(size_t i = 0; i < ev_wf->size(); i++){
      
      // Get a pointer to one of the TPC waveforms
      larlite::rawdigit* tpcwf = &(ev_wf->at(i));      

      UInt_t channel = tpcwf->Channel();
      
      // Check for empty waveforms!
      if( tpcwf->ADCs().size() < 1 ){
	print(msg::kERROR, __FUNCTION__,
	      Form("Found 0-length waveform: Event %d ... Ch. %d", ev_wf->event_id(), channel) );
	continue;
      }
      
      std::vector<short> ADCs = tpcwf->ADCs();
      for(size_t s = 1; s < ADCs.size() - 1; s++){
	_hDiffToInterpol->Fill( channel, ADCs[s] - (ADCs[s+1] + ADCs[s-1])/2. );
      }

    } // End of loop over all waveforms

    return true;
  }

  bool TPCWaveformAna::finalize() {

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    if(_fout){
      _fout->cd(); 
      _hDiffToInterpol->Write();
    }

    return true;
  }

}
#endif
