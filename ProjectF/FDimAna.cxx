#ifndef LARLITE_FDIMANA_CXX
#define LARLITE_FDIMANA_CXX

#include "FDimAna.h"

#include "LArUtil/GeometryHelper.h"
#include "DataFormat/hit.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"

#include <TCanvas.h>
#include <TF1.h>

namespace larlite {

  // TODO: Move to FDimAna.h
  enum FDimType { kLinFDim = 0, /// F dimension from linear fit
		  kLinCutoff, /// Cutoff scale from linear fit
		  kLinChiSq,  /// Chi square from linear fit
		  kMeanFDim, /// F dimension from average of slopes
		  kMeanCutoff, /// Cutoff scale from average of slopes
		  kPowFDim, /// F dimension from power fit
		  kPowCutoff, /// Cutoff scale from power fit
		  kPowChiSq, /// Chi square from power fit
		  kMaxFDim
  };

  bool FDimAna::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D").
    //

    gROOT->SetBatch();

    // Bin widths in m (must be in increasing order)
    binWidths = {0.004375, 0.00875, 0.0175, 0.035, 0.07, 0.14, 0.28, 0.56};

    hitHistos.resize(binWidths.size());
    boxNumbers.resize(binWidths.size());

    BDim.resize(kMaxFDim);
    BDimSlope.resize(binWidths.size() - 1);
    CDim.resize(kMaxFDim);
    CDimSlope.resize(binWidths.size() - 1);
    IDim.resize(kMaxFDim);
    IDimSlope.resize(binWidths.size() - 1);

    // TODO: set from python script
    _chargeMin = 3;

    // TODO: triplicate (for the 3 planes)? (directly or using a vector)
    size_t plane = 2; // Collection plane
    // TODO: define histogram limits based on detector (or ROI) dimensions (using Geometry)
    // Relative lenghts in m (measured from the minimum wire and minimum time events)
    double wireMin = 0;
    double wireMax = 11.;
    double timeMin = 0;
    double timeMax = 11.;

    // Create histograms
    for(size_t b = 0; b < binWidths.size(); b++){
      hitHistos.at(b) = new TH2D( Form( "hitHisto_plane%zu_binW%imm", plane, (int)(1000*binWidths.at(b)) ),
				  Form( "Hits on plane %zu / %f m #times %f m; Rel. wire (m); Rel. time (m)", plane, binWidths.at(b), binWidths.at(b) ),
				  int((wireMax - wireMin)/binWidths.at(b)), wireMin, wireMax,
				  int((timeMax - timeMin)/binWidths.at(b)), timeMin, timeMax );
    }

    boxProf = new TProfile( Form( "boxProf_plane%zu", plane ),
			    Form( "Boxes on plane %zu; log_{10}(#epsilon (m)); log_{10}(1/N(#epsilon))", plane ),
			    10, std::log10(binWidths.front()/2.), std::log10(2.*binWidths.back()), "s" );

    boxProfNorm = new TProfile( Form( "boxProfNorm_plane%zu", plane ),
				Form( "Boxes on plane %zu; #epsilon (m); N(%f m)/N(#epsilon)", plane, binWidths.front() ),
				200, binWidths.front(), binWidths.back() + 0.1*binWidths.back(), "s" );

    dimProf = new TProfile( Form( "dimProf_plane%zu", plane ),
			    Form( "Boxes on plane %zu; #epsilon (m); FD", plane ),
			    200, binWidths.front(), binWidths.back() + 0.1*binWidths.back(), "s" );

    boxH2 = new TH2D( Form( "boxH2_plane%zu", plane ),
		      boxProf->GetTitle(),
		      boxProf->GetNbinsX(), boxProf->GetXaxis()->GetXmin(), boxProf->GetXaxis()->GetXmax(),
		      std::rint(std::fabs(20*std::log10(1./(hitHistos.at(0)->GetNcells())))), std::log10(1./(hitHistos.at(0)->GetNcells())), 0 );

    boxH2Norm = new TH2D( Form( "boxH2Norm_plane%zu", plane ),
			  boxProfNorm->GetTitle(),
			  boxProfNorm->GetNbinsX(), boxProfNorm->GetXaxis()->GetXmin(), boxProfNorm->GetXaxis()->GetXmax(),
			  std::rint(std::pow(binWidths.back(), 2)/std::pow(binWidths.front(), 2) - 1),
			  1., std::pow(binWidths.back(), 2)/std::pow(binWidths.front(), 2) );

    dimH2 = new TH2D( Form( "dimH2_plane%zu", plane ),
		      dimProf->GetTitle(),
		      dimProf->GetNbinsX(), dimProf->GetXaxis()->GetXmin(), dimProf->GetXaxis()->GetXmax(),
		      300,0., 3. );

    corProf = new TProfile( Form( "corProf_plane%zu", plane ),
			    Form( "Correlation integral on plane %zu; log_{10}(#epsilon (m)); log_{10}(C(#epsilon))", plane ),
			    10, std::log10(binWidths.front()/2.), std::log10(2.*binWidths.back()), "s" );

    corProfNorm = new TProfile( Form( "corProfNorm_plane%zu", plane ),
				Form( "Correlation integral on plane %zu; #epsilon (m); N^{2}C(#epsilon)", plane ),
				200, binWidths.front(), binWidths.back() + 0.1*binWidths.back(), "s" );

    corH2 = new TH2D( Form( "corH2_plane%zu", plane ),
		      corProf->GetTitle(),
		      corProf->GetNbinsX(), corProf->GetXaxis()->GetXmin(), corProf->GetXaxis()->GetXmax(),
		      600, -6, 0 ); // hardcoded

    corH2Norm = new TH2D( Form( "corH2Norm_plane%zu", plane ),
		      corProfNorm->GetTitle(),
		      corProfNorm->GetNbinsX(), corProfNorm->GetXaxis()->GetXmin(), corProfNorm->GetXaxis()->GetXmax(),
		      600, -6, 0 ); // hardcoded

    infProf = new TProfile( Form( "infProf_plane%zu", plane ),
			    Form( "Information on plane %zu; log_{10}(#epsilon (m)); #sum_{i}^{N}P_{i}(#epsilon)log_{10}[P_{i}(#epsilon)]", plane ),
			    10, std::log10(binWidths.front()/2.), std::log10(2.*binWidths.back()), "s" );

    infH2 = new TH2D( Form( "infH2_plane%zu", plane ),
		      infProf->GetTitle(),
		      infProf->GetNbinsX(), infProf->GetXaxis()->GetXmin(), infProf->GetXaxis()->GetXmax(),
		      std::rint(std::fabs(20*std::log10(1./(hitHistos.at(0)->GetNcells())))), std::log10(1./(hitHistos.at(0)->GetNcells())), 0 );

    outTree = new TTree("outTree", "FDimAna output tree");
    outTree->Branch("boxN", &boxNumbers );
    outTree->Branch("BDim", &BDim );
    outTree->Branch("BDimSlope", &BDimSlope );
    outTree->Branch("CDim", &CDim );
    outTree->Branch("CDimSlope", &CDimSlope );
    outTree->Branch("IDim", &IDim );
    outTree->Branch("IDimSlope", &IDimSlope );
    outTree->Branch("hitSumInt", &hitSumInt );
    outTree->Branch("showerEDep", &showerDepE );
    outTree->Branch("showerETh", &showerTruthE );
    outTree->Branch("trackEDep", &trackDepE );
    outTree->Branch("trackETh", &trackTruthE );
    outTree->Branch("processTh", &processTruth );

    std::cout << "Finished initialize" << std::endl;

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

    int index = storage->get_index();

    resetOutTreeVars();

    event_mctrack* ev_mctrack = storage->get_data<event_mctrack>("mcreco");
    for(auto const& mctrack : *ev_mctrack) {
      trackTruthE += mctrack.Start().E();
      // Avoid segmentation violation by ill particles
      //std::cout << "Mass " << mctrack.front().Momentum().M() << std::endl;
      //if( mctrack.front().Momentum().M() < 1500 ){ // hardcoded
      //trackDepE += mctrack.front().E() - mctrack.back().E();
      //} else {
      trackDepE = -1; // Disable: the computation causes segmentation violations
	//break;
	//}
      if( mctrack.Process() != "primary" ) processTruth.push_back( mctrack.Process() );
      // std::cout << "PDG code: " << mctrack.PdgCode() << " Mother PDG code: " << mctrack.MotherPdgCode() << " Process: " << mctrack.Process() << std::endl;
    }

    event_mcshower* ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    for(auto const& mcshower : *ev_mcshower) {
      showerTruthE += mcshower.Start().E();
      showerDepE += mcshower.DetProfile().E();
      if( mcshower.Process() != "primary" ) processTruth.push_back( mcshower.Process() );
      // std::cout << "PDG code: " << mcshower.PdgCode() << " Mother PDG code: " << mcshower.MotherPdgCode() << " Process: " << mcshower.Process() << std::endl;
    }

    const larutil::GeometryHelper* geom = ::larutil::GeometryHelper::GetME();
    // TODO: Loop over 3 planes?
    // const size_t& nplanes = geom->Nplanes();
    // Wire to m conversion factor
    const float& wire2m = geom->WireToCm()/100.;
    // Time to m conversion factor
    const float& time2m = geom->TimeToCm()/100.;

    // TODO: Set producer from python script
    event_hit* ev_hit = storage->get_data<event_hit>("gaushit");
    if(ev_hit == nullptr){
      //throw DataFormatException("Could not locate hit data product!");
      return false;
    }
    //if(ev_hit->size() == 0){
      //throw DataFormatException("event_hit is empty");
      //return false;
    //}
    // TODO: Add support for ROI

    // Require at least 2 hits (needed for the correlation dimension)
    if(ev_hit->size() < 2) return false;

    // Take the minimum wire and minimum time as origin of each coordinate to ensure using the minimum of boxes to cover
    double minWire = 1e12, minTime = 1e12;
    for(auto const& hit : *ev_hit) {
      // TODO: use also the other planes?
      if( hit.WireID().Plane != 2 ) continue;

      if( hit.Integral() < _chargeMin ) continue;

      if( hit.WireID().Wire < minWire ) minWire = hit.WireID().Wire;
      if( hit.PeakTime() < minTime ) minTime = hit.PeakTime();
    }

    double maxHitWire = 0., maxHitTime = 0.;
    //    for(auto const& hit_h : *ev_hit) {
    // TODO: Set the maximum and number of bins using detector dimensions from initialize
    TH1D* corH = new TH1D( Form("corH_ev%i", index), "Correlation histogram; Distance (m); Hit pairs", 16000, 0, 16);

    for(size_t h = 0; h < ev_hit->size(); h++){
      hit hit_h = ev_hit->at(h);

      // Get plane
      auto const& plane = hit_h.WireID().Plane;

      // TODO: use also the other planes?
      if( plane != 2 ) continue;

      double hitInt = hit_h.Integral();
      if( hitInt < _chargeMin ) continue;

      // Hit coordinates
      double hitWire = (hit_h.WireID().Wire - minWire) * wire2m;
      double hitTime = (hit_h.PeakTime() - minTime) * time2m;
      //      double hitTime = hitWire; // Test with straight lines ---
      //      std::cout << "Wire: " << hitWire << " Time: " << hitTime << std::endl;

      // Find maximum hit and maximum wire coordinates to set the range when plotting the histogram
      if( hitWire > maxHitWire ) maxHitWire = hitWire;
      if( hitTime > maxHitTime ) maxHitTime = hitTime;

      for(size_t b = 0; b < hitHistos.size(); b++){
	//	for(double hitTime = hitWire; hitTime >= 0.; hitTime -= 0.003){// Test with triangles <|
	// Information method uses hit weights
	hitHistos.at(b)->Fill(hitWire, hitTime); // weight = 1 (default)
	// hitHistos.at(b)->Fill(hitWire, hitTime, hitInt); // Charge weights - useful to see charge distribution
	//	}// Test with triangles <|
      }
      hitSumInt += hitInt;

      // Correlation method
      for(size_t i = h + 1; i < ev_hit->size(); i++){
	hit hit_i = ev_hit->at(i);
	// TODO: use also the other planes?
	if( hit_i.WireID().Plane != 2 ) continue;
	if( hit_i.Integral() < _chargeMin ) continue;
	double hitWire_i = (hit_i.WireID().Wire - minWire) * wire2m;
	double hitTime_i = (hit_i.PeakTime() - minTime) * time2m;
	//	double hitTime_i = hitWire_i; // Test with straight lines ---
	//	for(double hitTime = hitWire; hitTime >= 0; hitTime -= 0.003){ // Test with triangles <|
	//	for(double hitTime_i = hitWire_i; hitTime_i >= 0.; hitTime_i -= 0.003){// Test with triangles <|
	double corDistance = std::sqrt(std::pow(hitWire - hitWire_i, 2) + std::pow(hitTime - hitTime_i, 2));
	if( corDistance != 0. ) corH->Fill(corDistance);
	//	}} // Test with trianges <|
      }
    }

    // Graphs
    //__________________________________________________________________________

    // Correlation method graphs
    TGraph* corG = new TGraph(binWidths.size()); // Same number of points as others to ease comparison, but not required
    TGraph* corGNorm = new TGraph(binWidths.size()); // Same number of points as others to ease comparison, but not required
    // Nevertheless, it is better to have the x points equally spaced (after logarithm), so in the end, the choice will be similar to binWidths
    double lastLogCorIntegral = 0;
    double lastLogDistance = 0;
    double corNorm = corH->Integral();
    for(size_t b = 0; b < binWidths.size(); b++){
      int bin = corH->FindBin(binWidths.at(b));
      double corIntegralN2 = corH->Integral(1, bin);
      double logCorIntegral = std::log10(corIntegralN2/corNorm);
      double distance = corH->GetBinCenter(bin + 1);
      double logDistance = std::log10(distance);
      corProf->Fill(logDistance, logCorIntegral);
      corH2->Fill(logDistance, logCorIntegral);
      corG->SetPoint(b, logDistance, logCorIntegral);
      corProfNorm->Fill(distance, corIntegralN2);
      corH2Norm->Fill(distance, corIntegralN2);
      corGNorm->SetPoint(b, distance, corIntegralN2);

      if(b > 0 ){
	// Compute local slopes
	double localSlope = (logCorIntegral - lastLogCorIntegral)/(logDistance - lastLogDistance);
	CDimSlope.at(b - 1) = localSlope;
      }
      lastLogCorIntegral = logCorIntegral;
      lastLogDistance = logDistance;
    }

    // int corFirstBin = corH->FindFirstBinAbove(0.);
    // int corLastBin = corH->FindLastBinAbove(0.);
    // int nBins = 0;
    // for(int bin = 2*corFirstBin; bin <= corLastBin; bin *= 2){ nBins++;}
    // TGraph* corG = new TGraph(nBins);
    // int point = 0;
    // for(int bin = 2*corFirstBin; bin <= corLastBin; bin *= 2){
    //   double corIntegral = corH->Integral(corFirstBin, bin); // Not divided by N^2 because it does not affect the slope when taking logs
    //   corG->SetPoint(point, std::log10(corH->GetBinCenter(bin + 1)), std::log10(corIntegral));
    //   point++;
    // }

    TGraph* boxGNorm = new TGraph(hitHistos.size() - 1);
    TGraph* boxG = new TGraph(hitHistos.size());

    TGraph* infG = new TGraph(hitHistos.size());

    double lastNegentropy = 0;

    for(size_t b = 0; b < hitHistos.size(); b++){
      // Write canvas for 1 of 100 events
      if( index%100 == 0 ){
	//if( 1 ){
	TCanvas* canvas = new TCanvas( Form("event%i_%s", index, hitHistos.at(b)->GetName()), hitHistos.at(b)->GetTitle() );
	// binarize(hitHistos.at(b));
	// Zoom to the event
	hitHistos.at(b)->GetXaxis()->SetRange(1, hitHistos.at(b)->GetXaxis()->FindBin(maxHitWire));
	hitHistos.at(b)->GetYaxis()->SetRange(1, hitHistos.at(b)->GetYaxis()->FindBin(maxHitTime));
	hitHistos.at(b)->Draw("colz");
	// canvas->SaveAs(".png");
	if( _fout ){
	  _fout->cd();
	  canvas->Write();
	}
	delete canvas;
      }

      // Box method graphs
      double boxes = countBoxes(hitHistos.at(b));
      boxNumbers.at(b) = boxes;
      double logInvBoxes = -std::log10(boxes); // == std::log10(1./boxes);
      double logBinWidth = std::log10(binWidths.at(b));
      if( boxes != 0 ){
	boxProf->Fill(logBinWidth, logInvBoxes);
	boxH2->Fill(logBinWidth, logInvBoxes);
	boxG->SetPoint(b, logBinWidth, logInvBoxes);

	dimProf->Fill( binWidths.at(b), logInvBoxes/logBinWidth );
	dimH2->Fill( binWidths.at(b), logInvBoxes/logBinWidth );

	if(b > 0 ){ // Do not waste time to compute a 1
	  boxProfNorm->Fill( binWidths.at(b), boxNumbers.at(0)/boxes );
	  boxH2Norm->Fill( binWidths.at(b), boxNumbers.at(0)/boxes );

	  // Per event
	  boxGNorm->SetPoint(b - 1, binWidths.at(b), boxNumbers.at(0)/boxes);

	  // Compute local slopes
	  // For slopes, there is no need to normalize
	  // And anyway, the normalization disappears due to logarithm properties
	  double localSlope = (std::log10(boxNumbers.at(b - 1)) - std::log10(boxes))
	    /(logBinWidth - std::log10(binWidths.at(b - 1)));
	  BDimSlope.at(b - 1) = localSlope;
	}
      }

      // Information method graphs
      double negentropy = mentropy(hitHistos.at(b));
      infProf->Fill(logBinWidth, negentropy);
      infH2->Fill(logBinWidth, negentropy);
      infG->SetPoint(b, logBinWidth, negentropy);

      if(b > 0 ){
	// Compute local slopes
	double localSlope = (negentropy - lastNegentropy)/(logBinWidth - std::log10(binWidths.at(b - 1)));
	IDimSlope.at(b - 1) = localSlope;
      }
      lastNegentropy = negentropy;
    }

    if( index%100 == 0 ){
      //if( 1 ){
      boxGNorm->SetNameTitle( Form( "boxGNorm_ev%i", index ), boxProfNorm->GetTitle() );
      boxG->SetNameTitle( Form( "boxG_ev%i", index ), boxProf->GetTitle() );
      corG->SetNameTitle( Form( "corG_ev%i", index ), corProf->GetTitle() );
      corGNorm->SetNameTitle( Form( "corGNorm_ev%i", index ), corProfNorm->GetTitle() );
      infG->SetNameTitle( Form( "infG_ev%i", index ), infProf->GetTitle() );
    }

    // F dimensions
    //__________________________________________________________________________

    // Mean method
    // Box
    meanFDim(BDimSlope, binWidths, BDim);
    // Correlation
    meanFDim(CDimSlope, binWidths, CDim);
    // Information
    meanFDim(IDimSlope, binWidths, IDim);

    // Power scaling fit
    // Box
    powFDim(boxGNorm, binWidths, BDim);
    // Correlation
    powFDim(corGNorm, binWidths, CDim);

    // Linear scaling fit
    // Box
    linFDim(0., boxG, binWidths, BDim);
    // Correlation
    linFDim(0., corG, binWidths, CDim);
    // Information
    linFDim(0., infG, binWidths, IDim);

    if( index%100 == 0 ){
      //if( 1 ){
      boxGNorm->Write();
      boxG->Write();
      corG->Write();
      corGNorm->Write();
      infG->Write();
      corH->Write();
    }
    delete boxGNorm;
    delete boxG;
    delete corG;
    delete corGNorm;
    delete infG;
    delete corH;

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
      corProf->Write();
      infProf->Write();

      boxH2->Write();
      boxH2Norm->Write();
      dimH2->Write();
      corH2->Write();
      infH2->Write();

      outTree->Write();
    }
    return true;
  }

  void FDimAna::resetOutTreeVars()
  {
    std::fill(BDim.begin(), BDim.end(), 0.);
    std::fill(BDimSlope.begin(), BDimSlope.end(), 0.);
    std::fill(CDim.begin(), CDim.end(), 0.);
    std::fill(CDimSlope.begin(), CDimSlope.end(), 0.);
    std::fill(IDim.begin(), IDim.end(), 0.);
    std::fill(IDimSlope.begin(), IDimSlope.end(), 0.);
    hitSumInt = showerDepE = showerTruthE = trackDepE = trackTruthE = 0.;
    processTruth.clear();
  }

  void FDimAna::binarize(TH2* h2)
  {
    // Keep original entries
    double entries = h2->GetEntries();
    // Includes underflows and overflows
    for(int i = 0; i < h2->GetNcells(); i++){
      if( h2->GetBinContent(i) > 0 ) h2->SetBinContent(i,1);
    }
    h2->SetEntries(entries);
  }

  double FDimAna::countBoxes(TH2* h2)
  {
    double count = 0;
    // Exclude underflows and overflows
    for(int i = 1; i <= h2->GetNbinsX(); i++){
      for(int j = 1; j <= h2->GetNbinsY(); j++){
	if( h2->GetBinContent( h2->GetBin(i,j) ) > 0 ) count++;
      }
    }
    return count;
  }

  double FDimAna::mentropy(TH2* h2)
  {
    double entropy = 0;
    // Exclude underflows and overflows
    double totHits = h2->GetSumOfWeights();
    for(int i = 1; i <= h2->GetNbinsX(); i++){
      for(int j = 1; j <= h2->GetNbinsY(); j++){
	double probability = (h2->GetBinContent( h2->GetBin(i,j) ))/totHits;
	if( probability > 0 ) entropy += probability*std::log10(probability);
      }
    }
    return entropy;
  }

  double FDimAna::scalingFitPow(double *x, double *par)
  {
    double f = par[0];
    if( x[0] < par[2] ) f *= std::pow( x[0], par[1] );
    else f *= std::pow( par[2], par[1] );
    return f;
  }

  double FDimAna::scalingFitLin(double *x, double *par)
  {
    double f;
    if( x[0] < par[0] ) f = par[1] * ( x[0] - par[0] );
    else f = par[2]; // par[2] is fixed
    return f;
  }

  void FDimAna::meanFDim(const std::vector<double>& slopes, const std::vector<double>& widths, std::vector<double>& FDim)
  {
    double sumSlope = 0;
    double countSlope = 0;
    for(size_t s = 0; s < slopes.size(); s++){
      if( slopes[s] != 0 ){
	sumSlope += slopes[s];
	countSlope++;
      } else {
	// If slope is zero, first point of interval is taken as cutoff
	FDim[kMeanCutoff] = widths.at(s);
	// If cutoff is found, drop the last slope since it might be biased <-- Preliminary results indicate this is not true
	// if( s > 0 ){
	//   sumSlope -= slopes[s - 1];
	//   countSlope--;
	// }
	break; // Stop after first cutoff is found
      }
    }
    FDim[kMeanFDim] = sumSlope/countSlope;
    // If no cutoff is found, set it at the limit
    if( FDim[kMeanCutoff] == 0 ) FDim[kMeanCutoff] = widths.back();
  }

  void FDimAna::powFDim(TGraph* graph, const std::vector<double>& widths, std::vector<double>& FDim)
  {
    TF1* f1 = new TF1( Form("fit_%s", graph->GetName()), FDimAna::scalingFitPow, widths.front(), widths.back(), 3);
    f1->SetParameter(1, ((FDim[kMeanFDim] != 0 )? FDim[kMeanFDim] : 1.) ); // Initialize the F dimension using the one computed in the mean method
    f1->SetParLimits(1, 0., 3.);
    f1->SetParameter(0, 2. * TMath::MaxElement(graph->GetN(), graph->GetY()));
    // f1->SetParLimits(0, 0., 1.e4);
    f1->SetParameter(2, ((FDim[kMeanCutoff] != 0)? FDim[kMeanCutoff] : widths.back()) ); // Initialize the cutoff using the one computed in the mean method
    f1->SetParLimits(2, widths.front(), 2.*widths.back());
    graph->Fit(f1, "QR+");
    FDim[kPowFDim] = f1->GetParameter(1);
    FDim[kPowCutoff] = f1->GetParameter(2);
    FDim[kPowChiSq] = f1->GetChisquare();
    delete f1;
  }

  void FDimAna::linFDim(const double& saturation, TGraph* graph, const std::vector<double>& widths, std::vector<double>& FDim)
  {
    // 3 parameters but 1 is fixed
    TF1* f2 = new TF1( Form("fit_%s", graph->GetName()), FDimAna::scalingFitLin, std::log10(widths.front()), std::log10(widths.back()), 3);
    // Initialize the F dimension using the one computed in the mean method
    f2->SetParameter(1, ((FDim[kMeanFDim] != 0 )? FDim[kMeanFDim] : 1.) );
    f2->SetParLimits(1, 0., 3.);
    // Initialize the cutoff using the one computed in the mean method
    f2->SetParameter(0, std::log10( ((FDim[kMeanCutoff] != 0)? FDim[kMeanCutoff] : widths.back())) );
    f2->SetParLimits(0, std::log10(widths.front()), std::log10(100.*widths.back()));
    // Fixed parameter
    f2->FixParameter(2, saturation);
    graph->Fit(f2, "QR+");
    FDim[kLinFDim] = f2->GetParameter(1);
    FDim[kLinCutoff] = std::pow(10., f2->GetParameter(0));
    FDim[kLinChiSq] = f2->GetChisquare();
    delete f2;
  }

}
#endif
