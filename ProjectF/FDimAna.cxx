#ifndef LARLITE_FDIMANA_CXX
#define LARLITE_FDIMANA_CXX

#include "FDimAna.h"

#include "LArUtil/GeometryHelper.h"
#include "DataFormat/hit.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>

namespace larlite {

  // TODO: Move to FDimAna.h?
  enum FDimType { kFDim = 0, /// F dimension from linear fit (logistic for correlation)
		  kErrFDim, /// Error on F dimension from linear fit (logistic for correlation)
		  kScale, /// Cutoff scale from linear fit (transtion point for correlation)
		  kChiSq,  /// Chi square from linear fit (logistic for correlation)
		  kMeanFDim, /// F dimension from average of slopes
		  kMeanCutoff, /// Cutoff scale from average of slopes
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
    // binWidths = {0.004375, 0.00875, 0.0175, 0.035, 0.07, 0.14, 0.28, 0.56};
    binWidths = {0.00875, 0.0175, 0.035, 0.07, 0.14};
    hitHistos.resize(binWidths.size());
    boxNumbers.resize(binWidths.size());

    // Correlation dimension
    for(double d = binWidths.front(); d <= binWidths.back(); d+= 0.003){
	  distWidths.push_back(d);
    }
    CDim.resize(kMaxFDim);
    CDimSlope.resize(distWidths.size() - 1);

    // 10 dimensions
    for(double q = 0; q <= 10; q++){
      QExp.push_back(q);
    }
    QDim.resize(QExp.size());
    QDimSlope.resize(QDim.size());
    for(size_t qi = 0; qi < QDim.size(); qi++){
      QDim.at(qi).resize(kMaxFDim);
      QDimSlope.at(qi).resize(binWidths.size() - 1);
    }

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

    QDimProf = new TProfile( Form( "QDimProf_plane%zu", plane ),
			     Form( "Generalized dimension spectrum on plane %zu; q; D_{q}", plane ),
			     QExp.back() - QExp.front() + 1, QExp.front() - 0.5, QExp.back() + 0.5, "s" );

    CDimProf = new TProfile( Form( "CDimProf_plane%zu", plane ),
			     Form( "Correlation dimension spectrum on plane %zu; q; D_{q}", plane ),
			     QExp.back() - QExp.front() + 1, QExp.front() - 0.5, QExp.back() + 0.5, "s" );

    outTree = new TTree("outTree", "FDimAna output tree");
    outTree->Branch("boxN", &boxNumbers );
    outTree->Branch("CDim", &CDim );
    outTree->Branch("CDimSlope", &CDimSlope );
    outTree->Branch("QDim", &QDim );
    outTree->Branch("QDimSlope", &QDimSlope );
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
    int each = 50; // Save one event every "each" events

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
      //double hitTime = hitWire; // Test with straight lines ---
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
	//double hitTime_i = hitWire_i; // Test with straight lines ---
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
    TGraph* corG = new TGraph(distWidths.size());
    double lastLogCorIntegral = 0;
    double lastLogDistance = 0;
    double corNorm = corH->Integral();

    for(size_t b = 0; b < distWidths.size(); b++){
      int bin = corH->FindBin(distWidths.at(b));
      double corIntegralN2 = corH->Integral(1, bin);
      double logCorIntegral = std::log10(corIntegralN2/corNorm);
      double distance = corH->GetBinCenter(bin + 1);
      double logDistance = std::log10(distance);
      corG->SetPoint(b, logDistance, logCorIntegral);

      if(b > 0 ){
	// Compute local slopes
	double localSlope = (logCorIntegral - lastLogCorIntegral)/(logDistance - lastLogDistance);
	CDimSlope.at(b - 1) = localSlope;
      }
      lastLogCorIntegral = logCorIntegral;
      lastLogDistance = logDistance;
    }

    // Generalized method graphs
    std::vector<TGraph*> genG(QDim.size());
    for(size_t qi = 0; qi < genG.size(); qi++){
      genG.at(qi) = new TGraph(hitHistos.size());
    }
    std::vector<double> lastGeneralizedNum(QDim.size());

    for(size_t b = 0; b < hitHistos.size(); b++){
      // Write canvas for 1 of 100 events
      if( index%each == 0 ){
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

      double logBinWidth = std::log10(binWidths.at(b));

      // Generalized method graphs
      for(size_t qi = 0; qi < QExp.size(); qi++){
	double probSumQ = probabilitySumQ(hitHistos.at(b), QExp.at(qi));
	if(probSumQ != 0){
	  // O-dimension is the box counting
	  if( QExp.at(qi) == 0. ) boxNumbers.at(b) = probSumQ;
	  double generalizedNum = 0;
	  if( QExp.at(qi) != 1. ) generalizedNum = std::log10(probSumQ)/(QExp.at(qi) - 1);
	  // Information dimension using L'Hôpital rule
	  else generalizedNum = probSumQ;
	  genG.at(qi)->SetPoint(b, logBinWidth, generalizedNum);

	  if(b > 0 ){
	    // Compute local slopes
	    double localSlope = (generalizedNum - lastGeneralizedNum.at(qi))/(logBinWidth - std::log10(binWidths.at(b - 1)));
	    QDimSlope.at(qi).at(b - 1) = localSlope;
	  }
	  lastGeneralizedNum.at(qi) = generalizedNum;
	}
      }

    }

    if( index%each == 0 ){
      //if( 1 ){
      corG->SetNameTitle( Form( "corG_ev%i", index ), "Correlation integral; log_{10}(#epsilon (m)); log_{10}(C(#epsilon))" );
      for(size_t qj = 0; qj < QExp.size(); qj++){
	genG.at(qj)->SetNameTitle( Form( "genG%.0lf_ev%i", QExp.at(qj), index ),
				   Form( "D_{%.0lf}; log_{10}(#epsilon (m)); #frac{1}{%.0lf}log_{10}#sum_{i}^{N}P^{%.0lf}_{i}(#epsilon)",
					 QExp.at(qj), QExp.at(qj) - 1, QExp.at(qj) ) );
      }
    }

    // F dimensions
    //__________________________________________________________________________

    // Mean method

    // Correlation
    meanFDim(CDimSlope, distWidths, CDim);
    // Generalized
    for(size_t qi = 0; qi < QDim.size(); qi++){
      meanFDim(QDimSlope.at(qi), binWidths, QDim.at(qi));
    }

    // Linear scaling fit

    // Generalized

    TGraphErrors* QDimGE = new TGraphErrors(QExp.size());

    for(size_t qi = 0; qi < QDim.size(); qi++){
      linFDim(0., genG.at(qi), binWidths, QDim.at(qi));
      if( QDim.at(qi).at(kFDim) > 0. && QDim.at(qi).at(kFDim) < 3. ){
	QDimProf->Fill(QExp.at(qi), QDim.at(qi).at(kFDim));
	QDimGE->SetPoint((int)qi, QExp.at(qi), QDim.at(qi).at(kFDim));
	QDimGE->SetPointError((int)qi, 0., QDim.at(qi).at(kErrFDim));
      }
    }

    // Logistic scaling fit

    // Correlation
    logiFDim(corG, distWidths, CDim);
    // For reference
    if( CDim.at(kFDim) > 0 && CDim.at(kFDim) < 3 ){
      CDimProf->Fill(2, CDim.at(kFDim));
    }

    if( index%each == 0 ){
      //if( 1 ){
      corG->Write();
      for(size_t qi = 0; qi < genG.size(); qi++){
	genG.at(qi)->Write();
      }
      QDimGE->SetNameTitle( Form("QDimGE_ev%i", index), "Generalized dimension spectrum; q; D_{q}" );
      QDimGE->Write();
      corH->Write();
    }
    delete corG;
    delete corH;
    delete QDimGE;
    for(size_t qi = 0; qi < genG.size(); qi++){
      delete genG.at(qi);
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

      QDimProf->Write();
      CDimProf->Write();

      outTree->Write();
    }
    return true;
  }

  void FDimAna::resetOutTreeVars()
  {
    std::fill(CDim.begin(), CDim.end(), 0.);
    std::fill(CDimSlope.begin(), CDimSlope.end(), 0.);
    for(size_t q = 0; q < QDim.size(); q++){
      std::fill(QDim.at(q).begin(), QDim.at(q).end(), 0.);
      std::fill(QDimSlope.at(q).begin(), QDimSlope.at(q).end(), 0.);
    }
    hitSumInt = showerDepE = showerTruthE = trackDepE = trackTruthE = 0.;
    processTruth.clear();
  }

  double FDimAna::probabilitySumQ(TH2* h2, double exp)
  {
    double probSumQ = 0;
    // Exclude underflows and overflows
    double totHits = h2->GetSumOfWeights();
    for(int i = 1; i <= h2->GetNbinsX(); i++){
      for(int j = 1; j <= h2->GetNbinsY(); j++){
	if( h2->GetBinContent( h2->GetBin(i,j) ) > 0 ){
	  // 0-dim
	  if( exp == 0 ) probSumQ++; // Save time for box
	  // 1-dim
	  else if( exp == 1 ){
	    double probability = (h2->GetBinContent( h2->GetBin(i,j) ))/totHits;
	    // Information dimension using L'Hôpital rule
	    probSumQ+= probability*std::log10(probability);
	  }
	  // q-dim (q != 0, 1)
	  else probSumQ += std::pow( (h2->GetBinContent( h2->GetBin(i,j) ))/totHits, exp );
	}
      }
    }
    return probSumQ;
  }

  double FDimAna::scalingFitLin(double *x, double *par)
  {
    double f;
    if( x[0] < par[0] ) f = par[1] * ( x[0] - par[0] );
    else f = par[2]; // par[2] is fixed
    return f;
  }

  double FDimAna::scalingFitLogi(double *x, double *par)
  {
    // Written so that par[1] is the FDim
    double f = par[2]/( 1 + exp(-4./par[2]*par[1]*( x[0] - par[0] )) ) - par[2];
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

  void FDimAna::linFDim(const double& saturation, TGraph* graph, const std::vector<double>& widths, std::vector<double>& FDim)
  {
    // 3 parameters but 1 is fixed
    TF1* f2 = new TF1( Form("fit_lin_%s", graph->GetName()), FDimAna::scalingFitLin, std::log10(widths.front()), std::log10(widths.back()), 3);
    // Initialize the F dimension using the one computed in the mean method
    f2->SetParameter(1, ((FDim[kMeanFDim] != 0 )? FDim[kMeanFDim] : 1.) );
    f2->SetParLimits(1, 0., 3.);
    // Initialize the cutoff using the one computed in the mean method
    f2->SetParameter(0, std::log10( ((FDim[kMeanCutoff] != 0)? FDim[kMeanCutoff] : widths.back())) );
    f2->SetParLimits(0, std::log10(widths.front()), std::log10(100.*widths.back()));
    // Fixed parameter
    f2->FixParameter(2, saturation);
    f2->SetLineColor(kGreen);
    graph->Fit(f2, "QR+");
    FDim[kFDim] = f2->GetParameter(1);
    FDim[kErrFDim] = f2->GetParError(1);
    FDim[kScale] = std::pow(10., f2->GetParameter(0));
    FDim[kChiSq] = f2->GetChisquare();
    delete f2;
  }

  void FDimAna::logiFDim(TGraph* graph, const std::vector<double>& widths, std::vector<double>& FDim)
  {
    // 3 parameters
    TF1* f = new TF1( Form("fit_logi_%s", graph->GetName()), FDimAna::scalingFitLogi, std::log10(widths.front()), std::log10(widths.back()), 3);
    double minx, miny;
    graph->GetPoint(0, minx, miny);
    // Initialize the F scale limit at the smallest scale measured
    f->SetParameter(0, minx);
    f->SetParLimits(0, std::log10(0.01*widths.front()), std::log10(100.*widths.back()));
    // Initialize the normalization with the F scale limit at the smallest scale measured
    double p2 = -2.*miny;
    double maxp2 = 7; // hardcoded: 1e7 hits
    f->SetParameter(2, p2);
    f->SetParLimits(2, 0., maxp2);
    // Initialize the F dimension using the one computed in the mean method
    f->SetParameter(1, ((FDim[kMeanFDim] != 0 )? FDim[kMeanFDim] : 1.) );
    f->SetParLimits(1, 0., 3.);
    f->SetLineColor(kMagenta);
    graph->Fit(f, "QR+");
    FDim[kFDim] = f->GetParameter(1);
    FDim[kErrFDim] = f->GetParError(1);
    FDim[kScale] = std::pow(10., f->GetParameter(0));
    FDim[kChiSq] = f->GetChisquare();
    delete f;
  }

}
#endif
