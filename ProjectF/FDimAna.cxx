#ifndef LARLITE_FDIMANA_CXX
#define LARLITE_FDIMANA_CXX

#include "FDimAna.h"

#include "LArUtil/GeometryHelper.h"
#include "DataFormat/hit.h"

#include <TCanvas.h>

namespace larlite {

  bool FDimAna::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    gROOT->SetBatch();

    // Bin widths in m
    binWidths = {0.001, 0.003, 0.01, 0.03, 0.1, 0.3};

    // Create histograms
    hitHistos.resize(binWidths.size());
    // TODO: triplicate (for the 3 planes)? (directly or using a vector)
    size_t plane = 2; // Collection plane
    // TODO: define histogram limits based on detector (or ROI) dimensions (using Geometry)
    // Detector dimensions in m
    double wireMin = 0;
    double wireMax = 11.;
    double timeMin = 0;
    double timeMax = 11.;
    for(size_t b = 0; b < binWidths.size(); b++){
      hitHistos.at(b) = new TH2D( Form( "hitHisto_plane%zu_binW%imm", plane, (int)(1000*binWidths.at(b)) ), 
				  Form( "Hits on plane %zu / %f m #times %f m; Wire (m); Time (m)", plane, binWidths.at(b), binWidths.at(b) ),
				  int((wireMax - wireMin)/binWidths.at(b)), wireMin, wireMax,
				  int((timeMax - timeMin)/binWidths.at(b)), timeMin, timeMax );
    }

    // std::cout << "hitHistos size: " << hitHistos.size() << std::endl;

    boxProf = new TProfile( Form( "boxProf_plane%zu", plane ), 
			    Form( "Boxes on plane %zu; log_{10}(1/#epsilon (m)); log_{10}(N(#epsilon))", plane ),
			    10, std::log10(1./(2.*binWidths.back())), std::log10(1/(binWidths.front()/2.)), "s");

    boxProfNorm = new TProfile( Form( "boxProfNorm_plane%zu", plane ), 
				Form( "Boxes on plane %zu; #epsilon (m); N(0.001 m)/N(#epsilon)", plane ),
				400, binWidths.front()/2., 2.*binWidths.back(), "s");

    dimProf = new TProfile( Form( "dimProf_plane%zu", plane ), 
			    Form( "Boxes on plane %zu; #epsilon (m); FD", plane ),
			    400, binWidths.front()/2., 2.*binWidths.back() );

    
    // Output tree to store histograms
    // USELESS TO VIEW HISTOGRAMS
    outTree = new TTree("outTree", "Output tree");
    for(size_t b = 0; b < hitHistos.size(); b++){
      outTree->Branch( hitHistos.at(b)->GetName(), hitHistos.at(b)->ClassName(), &(hitHistos.at(b)) );
    }

    return true;
  }

  bool FDimAna::analyze(storage_manager* storage) {

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
  
    const larutil::GeometryHelper* geom = ::larutil::GeometryHelper::GetME();
    // TODO: Loop over 3 planes?
    // const size_t& nplanes = geom->Nplanes();
    // Wire to m conversion factor
    const float& wire2m = geom->WireToCm()/100.;
    // Time to m conversion factor
    const float& time2m = geom->TimeToCm()/100.;

    // TODO: Set producer from python script
    event_hit* ev_hit = storage->get_data<event_hit>("gaushit");
    if(ev_hit == nullptr) throw DataFormatException("Could not locate hit data product!");
    
    // TODO: Add support for ROI

    for(auto const& hit : *ev_hit) {
      // Get plane
      auto const& plane = hit.WireID().Plane;

      // TODO: use also the other planes?
      if( plane != 2 ) continue;
	  
      // Hit coordinates
      double hitWire = hit.WireID().Wire * wire2m;
      double hitTime = hit.PeakTime() * time2m;
      //      std::cout << "Wire: " << hitWire << " Time: " << hitTime << std::endl;

      for(size_t b = 0; b < hitHistos.size(); b++){
	hitHistos.at(b)->Fill(hitWire, hitTime);
      }
    }

    int index = storage->get_index();

    for(size_t b = 0; b < hitHistos.size(); b++){
      TCanvas* canvas = new TCanvas( Form("event%i_%s", index, hitHistos.at(b)->GetName()), hitHistos.at(b)->GetTitle() );
      binarize(hitHistos.at(b));
      hitHistos.at(b)->Draw("col");
      // canvas->SaveAs(".png");
      if( _fout ){ 
	_fout->cd();
	canvas->Write();	
      }
      delete canvas;
      if( hitHistos.at(b)->Integral() != 0 ){
	boxProf->Fill( std::log10(1./binWidths.at(b)), std::log10(hitHistos.at(b)->Integral()) );
	boxProfNorm->Fill( binWidths.at(b), (hitHistos.at(0)->Integral())/(hitHistos.at(b)->Integral()) );
	dimProf->Fill( binWidths.at(b), -std::log10(hitHistos.at(b)->Integral())/(std::log10(binWidths.at(b))) );
      }
    }

    outTree->Fill();

    for(size_t b = 0; b < hitHistos.size(); b++){
      hitHistos.at(b)->Reset();
    }

    return true;
  }

  bool FDimAna::finalize() {

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

    if( _fout ){ 
      _fout->cd();
      boxProf->Write();
      boxProfNorm->Write();
      dimProf->Write();
      outTree->Write();
    }
    return true;
  }

  void FDimAna::binarize(TH2* h2) {
    // Keep original entries 
    double entries = h2->GetEntries();
    for(int i = 0; i < h2->GetNcells(); i++){
      if( h2->GetBinContent(i) > 0 ) h2->SetBinContent(i,1);
    }
    h2->SetEntries(entries);
  }

}
#endif
