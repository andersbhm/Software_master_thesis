
// Written by Anders Lindanger. Master thesis in spacephysics 2018


void integral_hist(TH1D *h, double xmin, double xmax, double &integral){
  TAxis *axis = h->GetXaxis();
  int bmin = axis->FindBin(xmin); //in your case xmin=-1.5
  int bmax = axis->FindBin(xmax); //in your case xmax=0.8
  integral = h->Integral(bmin,bmax);
  integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
  integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/axis->GetBinWidth(bmax);
}

// Create tree of unique WWLLN time from photon_WWLLN_match.root
void create_f_WWLLN_unique_time(TString path, TString name_year){


  double t_W{0}, t, distance_footprint{0}, distance_lightning_AGILE{0}, propagation_time{0}, obt{0};
  int  i{0}, usec{0}, usec_p{0}, year_WWLLN{0}, month_WWLLN{0}, day_WWLLN{0}, hour_WWLLN{0}, minute_WWLLN{0}, sec_WWLLN{0};
  float lat_WWLLN{0}, lon_WWLLN{0}, height_A{0}, lat_A{0}, lon_A{0};

  cout.precision(15);

  gROOT->Reset();
  struct WWLLN_timing {
  };

  WWLLN_timing tdata;

  TFile *f_WWLLN_unique_time = new TFile(path + "f_WWLLN_unique_time" + name_year + ".root","RECREATE");
  TTree *tree_WWLLN = new TTree("tdata","");
  tree_WWLLN->Branch("t",&t_W, "t/D");
  tree_WWLLN->Branch("usec",&usec_p, "usec/I");
  tree_WWLLN->Branch("lat_W",&lat_WWLLN, "lat_W/F");
  tree_WWLLN->Branch("lon_W",&lon_WWLLN, "lon_W/F");
  tree_WWLLN->Branch("year_W", &year_WWLLN, "year_WWLLN/I");
  tree_WWLLN->Branch("month_W",&month_WWLLN, "month_W/I");
  tree_WWLLN->Branch("day_W", &day_WWLLN, "day_W/I");
  tree_WWLLN->Branch("hour_W",&hour_WWLLN, "hour_W/I");
  tree_WWLLN->Branch("minute_W", &minute_WWLLN, "minute_W/I");
  tree_WWLLN->Branch("second_W", &sec_WWLLN, "second_W/I");

  //#################################################################################


  //TFile *f_E = new TFile(path + "photon_WWLLN_match_without_duplicates_obt_folder/" + "photon_WWLLN_match_" + name_year +"_1_without_duplicates_obt" + ".root","r");
  //TTree *t_E = (TTree *) f_E->Get("data");

  TChain *t_E = new TChain("data");
  t_E->Add(path + "photon_WWLLN_match_without_duplicates_obt_folder/" + "*.root");

  t_E->SetBranchAddress("t_W", &t_W);
  t_E->SetBranchAddress("usec_W", &usec_p);
  t_E->SetBranchAddress("lat_W", &lat_WWLLN);
  t_E->SetBranchAddress("lon_W", &lon_WWLLN);
  t_E->SetBranchAddress("year_W", &year_WWLLN);
  t_E->SetBranchAddress("month_W",&month_WWLLN);
  t_E->SetBranchAddress("day_W", &day_WWLLN);
  t_E->SetBranchAddress("hour_W",&hour_WWLLN);
  t_E->SetBranchAddress("minute_W", &minute_WWLLN);
  t_E->SetBranchAddress("second_W", &sec_WWLLN);
  t_E->SetBranchAddress("propagation_time", &propagation_time);
  t_E->SetBranchAddress("obt_MCAL", &obt);

  int entries = t_E->GetEntries();
  //#################################################################################
  int counter = 0;
  double t_W_before{0};
  double dt;
  //Iterate through tree and make a list of all unique WWLLN times.
  set<double> event_list;
  while (t_E->GetEntry(i++)) {
    dt = obt - propagation_time - t_W;
    if (counter > entries/100) {
      cout << int((i/(float(entries + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;

    if (-0.1 < dt && dt < 0.1)  {

    if (event_list.count(t_W) == 0) {



      event_list.insert(t_W);
      //if ((t_W_before - 0.1) < t_W && t_W < (t_W_before + 0.1)) {
      //cout << t_W_before << "\t" << t_W << "\n"<< endl;
      //}
      tree_WWLLN->Fill();
      t_W_before = t_W;
    }
    }
  }

  tree_WWLLN->Print();
  //tree_WWLLN->Scan("","","col=20.6f");
  //delete f_E;
  f_WWLLN_unique_time->Write();
}

//#################################################################################
Double_t func_ostgaard_fit(Double_t*x, Double_t*par){
  return par[0]*TMath::Poisson(x[0],par[1]) + par[2]*pow(x[0],-par[3]);
  //return par[0]*TMath::Poisson(x[0],par[1]);
}

Double_t func_poiss(Double_t*x, Double_t*par){
  return par[0]*TMath::Poisson(x[0],par[1]);
}

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Lorentzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
   return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,
   (x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
   return background(x,par) + lorentzianPeak(x,&par[3]);
}
//#################################################################################

int background_hist(TString path, TString title_name, TString name_year){
  double propagation_time{0}, obt ,t_W{0}, unique_WWLLN_time{0}, avg_obt{0}, delta_time{0}, delta_time_check{0}, delta_time_fill{0}, obt_check{0}, dt{0}, dt_background{0};
  float energy{0}, lon_WWLLN{0}, lat_WWLLN{0}, unique_lat_WWLLN{0}, unique_lon_WWLLN{0}, lat_pos_A{0}, lon_pos_A{0}, h_pos_A{0};
  int entry_a{0}, entry_b{0}, minimum{0}, maximum{0}, x{0}, y{0}, z{0}, i{0}, j{0}, k{0}, usec{0}, usec_W{0}, bincontent{0}, xmin{0}, xmax{0}, fill_histogram{0}, contact{0};
  int year_WWLLN{0}, month_WWLLN{0}, day_WWLLN{0}, hour_WWLLN{0}, minute_WWLLN{0}, sec_WWLLN{0};
  cout.precision(16);

  // parameters for light curve histogram h1
  double binsize = 0.000100;

  // parameters for histogram hmain
  int start = 0;
  int stop = 20;
  int nbins2 = (stop-start)/1;
  int counter =0;


  // ############################ BACKGROUND ##################################
  double adt_background1 = 3.2;
  double adt_background2 = 3.7;
  int nbins_background = 0.5/binsize;

  double obt_MCAL_background, time_WWLLN_background, propagation_time_background,unique_background_time;
  //TFile *f_background = new TFile(path + "background_" + name_year + "_without_duplicates_obt.root","r");
  //TTree *tree_background = (TTree *) f_background->Get("bdata");

  TChain *tree_background = new TChain("bdata");
  tree_background->Add(path + "background/" + "*.root");

  tree_background->SetBranchAddress("obt_MCAL",&obt_MCAL_background);
  tree_background->SetBranchAddress("time_WWLLN",&time_WWLLN_background);
  tree_background->SetBranchAddress("propagation_time", &propagation_time_background);
  tree_background->BuildIndex("time_WWLLN", "0");


  TFile *file_histo = new TFile(path + "background_hist" + name_year +".root", "RECREATE");

  TH2D *hmain_background = new TH2D("3D_histo_background", "3D_histo_background", nbins2, start, stop, nbins_background, adt_background1, adt_background2);



  //TFile *f_background_unique_time = new TFile(path + "f_background_unique_time" + name_year + ".root","r");
  //TTree *tree_unique_background = (TTree *) f_background_unique_time->Get("tdata");
  TChain *tree_unique_background = new TChain("tdata");
  tree_unique_background->Add(path + "f_background_unique_time/" + "*.root");
  tree_unique_background->SetBranchAddress("t",&unique_background_time);
  int entries_unique_background = tree_unique_background->GetEntries();



  TH1D *h1_lightcurve_background[10];
  char *histname1_background = new char[10];
  TH2D *h_binheight_background[10];
  char *histname2_background = new char[10];
  TCanvas *c_background[10];
  char *canvasname_background = new char[10];  // parameters for light curve histogram h1
  int histo_name_background = 0;
  int canvas_counter_background = 0;

  counter = 0;

  int time_2016_03 = 385777915 - 55;
  int time_2016_06 = 393812580;
  int time_2017_03 = 417312419 - 59;
  int time_2017_06 = 425275801 - 239;

  int breaker = 0;
  int lightning_counter = 0;
  cout << "Iterating through background file" <<endl;
  for (x = 0 ; x < entries_unique_background; x++) {


    if (counter > entries_unique_background/100) {
      cout << int((x/(float(entries_unique_background + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;

    tree_unique_background->GetEntry(x);


    lightning_counter = lightning_counter +1 ;

    sprintf(histname1_background,"histo_background_%d_%f",histo_name_background, unique_background_time);
    h1_lightcurve_background[histo_name_background] = new TH1D(histname1_background, histname1_background, nbins_background, adt_background1, adt_background2);

    minimum = tree_background->GetEntryNumberWithBestIndex(unique_background_time) ;


    tree_unique_background->GetEntry(x + 1);
    maximum = tree_background->GetEntryNumberWithBestIndex(unique_background_time) -1;

    //tree_unique_background->GetEntry(x - 1);

    // Fill each lightcurve histogram
    for (y = minimum;  y < maximum; y++) {
      tree_background->GetEntry(y);
      h1_lightcurve_background[histo_name_background]->Fill(obt_MCAL_background - time_WWLLN_background - propagation_time_background);
    }

    sprintf(histname2_background,"hist2_background%d",histo_name_background);
    h_binheight_background[histo_name_background] = new TH2D(histname2_background, histname2_background, nbins2, start, stop, nbins_background, adt_background1, adt_background2);


    for (k = 0; k < nbins_background; k++) {
      bincontent = h1_lightcurve_background[histo_name_background]->GetBinContent(k);
      // Get value along x axis. Return the x value to the right of the bin
      h_binheight_background[histo_name_background]->Fill(bincontent, h1_lightcurve_background[histo_name_background]->GetXaxis()->GetBinCenter(k) - binsize/2);
    }
    hmain_background->Add(h_binheight_background[histo_name_background]);


    delete  h1_lightcurve_background[histo_name_background];
    delete h_binheight_background[histo_name_background];
    histo_name_background += 1;

  }

  cout << "Do projection_poisson_background" << endl;
  TH1D * projection_poisson_background = hmain_background->ProjectionX("projection_poisson_background");

  cout << "lightning_counter: " << lightning_counter<< endl;

  file_histo->Write();

  return 0;
}

void draw_histogram_and_find_candidates(TString path, TString title_name, TString name_year, double adt, int threshold, double plussminusinsideTGFcandidate){

  double propagation_time{0}, obt ,t_W{0}, unique_WWLLN_time{0}, avg_obt{0}, delta_time{0}, delta_time_check{0}, delta_time_fill{0}, obt_check{0}, dt{0}, dt_background{0}, distance_footprint{0};
  float energy{0}, lon_WWLLN{0}, lat_WWLLN{0}, unique_lat_WWLLN{0}, unique_lon_WWLLN{0}, lat_pos_A{0}, lon_pos_A{0}, h_pos_A{0};
  int entry_a{0}, entry_b{0}, minimum{0}, maximum{0}, x{0}, y{0}, z{0}, i{0}, j{0}, k{0}, usec{0}, usec_W{0}, bincontent{0}, xmin{0}, xmax{0}, fill_histogram{0}, contact{0};
  int year_WWLLN{0}, month_WWLLN{0}, day_WWLLN{0}, hour_WWLLN{0}, minute_WWLLN{0}, sec_WWLLN{0};
  cout.precision(16);

  gROOT->Reset();

  TFile *f_candidate = new TFile(path + "f_candidate_" + name_year + "_" + std::to_string(plussminusinsideTGFcandidate) + "_" + std::to_string(threshold) + ".root","RECREATE");
  TTree *tree_candidate = new TTree("tdata","TGF candidates time drift correction");
  tree_candidate->Branch("t_WWLLN",&unique_WWLLN_time, "t_WWLLN/D");
  tree_candidate->Branch("year_WWLLN", &year_WWLLN, "year_WWLLN/I");
  tree_candidate->Branch("month_WWLLN",&month_WWLLN, "month_WWLLN/I");
  tree_candidate->Branch("day_WWLLN", &day_WWLLN, "day_WWLLN/I");
  tree_candidate->Branch("hour_WWLLN",&hour_WWLLN, "hour_WWLLN/I");
  tree_candidate->Branch("minute_WWLLN", &minute_WWLLN, "minute_WWLLN/I");
  tree_candidate->Branch("sec_WWLLN", &sec_WWLLN, "sec_WWLLN/I");
  tree_candidate->Branch("usec_WWLLN",&usec, "usec_WWLLN/I");
  tree_candidate->Branch("dt",&dt, "dt/D");
  tree_candidate->Branch("contact", &contact, "contact/I");
  tree_candidate->Branch("lat_W", &lat_WWLLN, "lat_W/F");
  tree_candidate->Branch("lon_W", &lon_WWLLN, "lon_W/F");
  tree_candidate->Branch("lat_A", &lat_pos_A, "lat_A/F");
  tree_candidate->Branch("lon_A", &lon_pos_A, "lon_A/F");
  tree_candidate->Branch("height_A", &h_pos_A, "height_A/F");
  tree_candidate->Branch("bincontent",&bincontent, "bincontent/I");
  tree_candidate->Branch("propagation_time", &propagation_time, "propagation_time/D");
  tree_candidate->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");

  //#################################################################################

  TChain *t_E = new TChain("data");
  t_E->Add(path + "photon_WWLLN_match_without_duplicates_obt_folder/" + "*.root");
  t_E->SetBranchAddress("obt_MCAL", &obt);
  t_E->SetBranchAddress("t_W", &t_W);
  t_E->SetBranchAddress("usec_W", &usec_W);
  t_E->SetBranchAddress("year_W", &year_WWLLN);
  t_E->SetBranchAddress("month_W",&month_WWLLN);
  t_E->SetBranchAddress("day_W", &day_WWLLN);
  t_E->SetBranchAddress("hour_W",&hour_WWLLN);
  t_E->SetBranchAddress("minute_W", &minute_WWLLN);
  t_E->SetBranchAddress("second_W", &sec_WWLLN);
  t_E->SetBranchAddress("propagation_time", &propagation_time);
  //t_E->SetBranchAddress("Etot", &energy);
  //t_E->SetBranchAddress("lat_W", &lat_WWLLN);
  //t_E->SetBranchAddress("lon_W", &lon_WWLLN);
  t_E->SetBranchAddress("distance_footprint", &distance_footprint);

  t_E->SetBranchAddress("contact", &contact);
  t_E->SetBranchAddress("lat_A", &lat_pos_A);
  t_E->SetBranchAddress("lon_A", &lon_pos_A);
  t_E->SetBranchAddress("height_A", &h_pos_A);
  int entries = t_E->GetEntries();
  cout <<entries<< endl;

  cout << "Building index photon_WWLLN_match..." << endl;
  t_E->BuildIndex("t_W", "usec_W");
  cout << "Building index photon_WWLLN_match DONE" << endl;

  //#################################################################################
  cout << "Get f_WWLLN_unique_time" << endl;
  TFile *f_WWLLN_unique_time = new TFile(path + "f_WWLLN_unique_time" + name_year + ".root","r");
  TTree *tree_WWLLN = (TTree *) f_WWLLN_unique_time->Get("tdata");
  tree_WWLLN->SetBranchAddress("t", &unique_WWLLN_time);
  tree_WWLLN->SetBranchAddress("usec", &usec);
  tree_WWLLN->SetBranchAddress("lat_W", &lat_WWLLN);
  tree_WWLLN->SetBranchAddress("lon_W", &lon_WWLLN);
  int entries_unique_WWLLN = tree_WWLLN->GetEntries();
  //################################################################################

  /*
  //################################################
  double time_WWLLN_match;
  // Martinos WWLLN matches
  TFile *g = new TFile(path + "trg_PAPER_with_simultaneous_WWLLN_TIMING_ANGLES.root", "r");
  TTree *t_WWLLN_match = (TTree*)g->Get("triggers");
  t_WWLLN_match->SetBranchAddress("t0", &time_WWLLN_match);
  int entries_WWLLN_match = t_WWLLN_match->GetEntriesFast();
  int break_variable = 0;
  //################################################
  */

  TH1D *h1_lightcurve[10];
  char *histname1 = new char[10];
  TH2D *h_binheight[10];
  TH2D *h_binheight_signal[10];

  char *histname2 = new char[10];
  char *histname3 = new char[10];

  TCanvas *c[10];
  char *canvasname = new char[10];  // parameters for light curve histogram h1

  // parameters for light curve histogram h1
  double binsize = 0.000100;
  int nbins1 = (int) 2*adt/binsize;
  double adt_signal = adt;
  int nbins1_signal = (int) 2*adt_signal/binsize;

  // parameters for histogram hmain
  int start = 0;
  int stop = 20;
  int nbins2 = (stop-start)/1;
  cout << "Create stackplots_histogram" << endl;

  TFile *file_histo = new TFile(path + "stackplots_histogram_" + name_year + "_" + std::to_string(plussminusinsideTGFcandidate) + "_" + std::to_string(threshold) + ".root", "RECREATE");
  TH2D *hmain = new TH2D("3D_histo", "3D_histo_title", nbins2, start, stop, nbins1, -adt, adt);
  TH2D *hsignal = new TH2D("3D_histo_signal", "3D_histo_signal_title", nbins2, start, stop, nbins1_signal, -adt_signal , adt_signal);

  y = 0;
  set<double> event_list;

  int histo_name = 0;
  int canvas_counter = 0;
  int counter =  0;
  cout << "Iterating through photon file" <<endl;
  for (x = 0 ; x < entries_unique_WWLLN; x++) {




    if (counter > entries_unique_WWLLN/100) {
      cout << int((x/(float(entries_unique_WWLLN + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;

    tree_WWLLN->GetEntry(x);
    /*
    // ################################################
    break_variable = 0;
    for (int j = 0; j < entries_WWLLN_match; j++) {
      t_WWLLN_match->GetEntry(j);
      if (TMath::Abs(time_WWLLN_match - unique_WWLLN_time) < 0.005) {
        break_variable = 1;
        cout << time_WWLLN_match << " " << unique_WWLLN_time << endl;
      }
    }

    if (break_variable == 1) {
      continue;
    }

    //################################################
    */
    sprintf(histname1,"histo_%d_%f_%d",histo_name, unique_WWLLN_time, contact);
    h1_lightcurve[histo_name] = new TH1D(histname1, histname1, nbins1, -adt, adt);

    minimum = t_E->GetEntryNumberWithBestIndex(unique_WWLLN_time, usec) ;

    tree_WWLLN->GetEntry(x + 1);
    maximum = t_E->GetEntryNumberWithBestIndex(unique_WWLLN_time, usec) -1;

    //tree_WWLLN->GetEntry(x - 1);

    // Fill each lightcurve histogram
    for (y = minimum;  y < maximum; y++) {
      t_E->GetEntry(y);
      h1_lightcurve[histo_name]->Fill(obt - t_W - propagation_time);
    }

    // Create new histogram to see for bin height
    sprintf(histname2,"hist2%d",histo_name);
    sprintf(histname3,"hist2%d_signal",histo_name);


    h_binheight[histo_name] = new TH2D(histname2, histname2, nbins2, start, stop, nbins1, -adt, adt);
    h_binheight_signal[histo_name] = new TH2D(histname3, histname3, nbins2, start, stop, nbins1_signal, -adt_signal, adt_signal);

    // Iterate through lightcurve and look at dt = obt - t_W - propagation_time and binheight (bincontent)
    for (k = 0; k < nbins1; k++) {
      bincontent = h1_lightcurve[histo_name]->GetBinContent(k);
      // Get value along x axis. Return the x value to the right of the bin
      dt =  h1_lightcurve[histo_name]->GetXaxis()->GetBinCenter(k) - binsize/2 ;

      if (TMath::Abs(dt) < plussminusinsideTGFcandidate) {
        h_binheight[histo_name]->Fill(bincontent, dt);
      }

      if (bincontent >= threshold) {


        // To not draw the same single lightcurve twice assiciated with one WWLLN
        if (event_list.count(unique_WWLLN_time) == 0 && TMath::Abs(dt) < plussminusinsideTGFcandidate) {
          tree_WWLLN->GetEntry(x);
          t_E->GetEntryWithIndex(unique_WWLLN_time, usec);
          event_list.insert(unique_WWLLN_time);
          tree_candidate->Fill();
          h_binheight_signal[histo_name]->Fill(bincontent, dt);


          //if (361727338.581796 - 1 < obt && obt < 361727338.581796 +1 ) { // 3107-3113
          //sprintf(canvasname,"c%d_%d_%f",canvas_counter, histo_name, unique_WWLLN_time);
          //  cout << obt << "   " <<unique_WWLLN_time << "   " << contact << "  " <<  h1_lightcurve[histo_name]->GetBinContent(k) <<endl;
          //c[histo_name] = new TCanvas(canvasname,canvasname,0 + canvas_counter*500, 0, 500, 400);
          //h1_lightcurve[histo_name]->Draw();
          //canvas_counter += 1;
          //}

        }
      }
    }

    hmain->Add(h_binheight[histo_name]);
    hsignal->Add(h_binheight_signal[histo_name]);
    delete  h1_lightcurve[histo_name];
    delete h_binheight[histo_name];
    delete h_binheight_signal[histo_name];
    histo_name += 1;
  }

  TH1D * projection_poisson = hmain->ProjectionX("projection_poisson_3d_histo");
  //TH1D * projection_poisson = hsignal->ProjectionX();

  TH1D * projection_light_curve = hmain->ProjectionY("lightcurve_over_threshold", threshold, stop);
  TH1D * projection_light_curve_signal = hsignal->ProjectionY("lightcurve_over_threshold_signal", threshold, stop);


  // ############################ BACKGROUND ##################################
  double adt_background1 = 3.2;
  double adt_background2 = 3.7;
  int nbins_background = 0.5/binsize;

  TFile *file_hist = new TFile(path + "background_hist" + name_year + ".root", "r");
  TH1D *projection_poisson_background = (TH1D*)file_hist->Get("projection_poisson_background");

  //tree_candidate->Print();
  //tree_candidate->Scan("t_WWLLN:propagation_time","","col=20.6f");
  cout << "Write f_candidate" << endl;
  f_candidate->Write();

  cout << "Write file_histo" << endl;
  file_histo->Write();
  cout << "Write file_histo DONE" << endl;

  cout << "Number of candidates: " <<tree_candidate->GetEntries() << endl;

}

void plot(TString path, TString title_name, TString name_year, double adt, int threshold, double plussminusinsideTGFcandidate){



  TFile *f_candidate = new TFile(path + "f_candidate_" + name_year + "_" + std::to_string(plussminusinsideTGFcandidate) + "_" + std::to_string(threshold) + ".root","r");
  TTree *tree_candidate = (TTree *) f_candidate->Get("tdata");

  int candidates = tree_candidate->GetEntries();

  // parameters for light curve histogram h1
  double binsize = 0.000100;
  int nbins1 = (int) 2*adt/binsize;
  double adt_signal = adt;
  int nbins1_signal = (int) 2*adt_signal/binsize;

  // parameters for histogram hmain
  int start = 0;
  int stop = 20;
  int nbins2 = (stop-start)/1;

  double adt_background1 = 3.2;
  double adt_background2 = 3.7;
  int nbins_background = 0.5/binsize;


  TFile *histo = new TFile(path + "stackplots_histogram_" + name_year + "_" + std::to_string(plussminusinsideTGFcandidate) + "_" + std::to_string(threshold) + ".root", "r");
  TH2D *hmain = (TH2D*) histo->Get("3D_histo");
  TH1D *projection_poisson = hmain->ProjectionX("projection_poisson_3d_histo");
  TH2D *hsignal = (TH2D*) histo->Get("3D_histo_signal");
  TH1D * projection_light_curve_signal = hsignal->ProjectionY("lightcurve_over_threshold_signal");

  TFile *file_hist = new TFile(path + "background_hist" + name_year + ".root", "r");
  TH1D *projection_poisson_background = (TH1D*)file_hist->Get("projection_poisson_background");
  TH2D *histo_background = (TH2D*)file_hist->Get("3D_histo_background");
  TH1D * projection_light_curve_backround = histo_background->ProjectionY("lightcurve_over_threshold_background", threshold, stop);


  TString path_2015 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_disabled_analyse/";

  TFile *histo_ref_period = new TFile(path_2015 + "stackplots_histogram_" + "2015" + "_" + "0.000500" + "_" + "4" + ".root", "r");
  TH2D *hmain_ref_period = (TH2D*) histo_ref_period->Get("3D_histo");

  TH1D *projection_poisson_ref_period = hmain_ref_period->ProjectionX("projection_poisson_3d_histo");

  TFile *file_hist_ref_period = new TFile(path_2015 + "background_hist" + "2015" + ".root", "r");
  TH1D *projection_poisson_background_ref_period = (TH1D*)file_hist_ref_period->Get("projection_poisson_background");


  // Get mean and sigma of Background

  int binheight, k;
  int nbins = projection_light_curve_backround->GetSize()  ; // included Underflow overflow
  double binheight_sum = 0;
  for (k = 1; k < nbins; k++) {
    binheight = projection_light_curve_backround->GetBinContent(k);
    binheight_sum += binheight;
  }

  double average_counts_per_bin = binheight_sum/(nbins-2);
  //cout << "nbins: " << nbins << endl;
  //cout << "binheight_sum: " << binheight_sum << endl;
  cout << "average_counts_per_bin: " << average_counts_per_bin << endl;

  double x_x_mean_square_sum = 0;
  for (k = 1; k < nbins; k++) {
    binheight = projection_light_curve_backround->GetBinContent(k);
    x_x_mean_square_sum += pow(binheight - average_counts_per_bin,2);
  }

  double variance = x_x_mean_square_sum/(nbins-1);
  double sigma = TMath::Sqrt(variance);
  //cout << "x_x_mean_square_sum: " << x_x_mean_square_sum << endl;
  //cout << "variance: " << variance << endl;
  cout << "sigma: " << sigma << endl;


  TF1 *f_mean = new TF1("f_mean","[0]",-2,2);
  f_mean->SetParameter(0,average_counts_per_bin);
  f_mean->SetParName(0,"mean");

  TF1 *f_sigma1 = new TF1("f_sigma1","[0]",-2,2);
  f_sigma1->SetParameter(0,average_counts_per_bin+sigma);
  f_sigma1->SetParName(0,"sigma1");

  TF1 *f_sigma3 = new TF1("f_sigma3","[0]",-2,2);
  f_sigma3->SetParameter(0,average_counts_per_bin+sigma*3);
  f_sigma3->SetParName(0,"sigma2");


  TF1 *f_sigma5 = new TF1("f_sigma5","[0]",-2,2);
  f_sigma5->SetParameter(0,average_counts_per_bin+sigma*5);
  f_sigma5->SetParName(0,"sigma5");




  //TCanvas *can3_background = new TCanvas("can3_background","can3_background",500*1.5,500, 500*1.5, 400*1.5);
  //can->SetLogy();
  gStyle->SetOptStat(1111);
  projection_light_curve_backround->Draw();
  //double integral_background = 0;
  //integral_hist(projection_light_curve_backround, double(3.5), double(3.8), integral_background);
  //cout <<  integral_background<< endl;
  //hmain->SetMaximum(20);
  /*  hmain_background->GetXaxis()->SetTitle("#particles per x sec time bin");
  hmain_background->GetYaxis()->SetTitle("dt");
  hmain_background->GetXaxis()->CenterTitle();
  hmain_background->GetYaxis()->SetTitleOffset(0);
  hmain_background->GetYaxis()->CenterTitle();
  hmain_background->GetXaxis()->SetTitleOffset(1.2);
  hmain_background->Draw("COLZ");
  can3_background->Print(path + "plots/TH3D_background" + title_name + to_string(plussminusinsideTGFcandidate) + ".gif");
  */


  char text2[200];
  sprintf(text2,"#particles per %.2f ms time bin summed for each single lightcurve +- %.2f ms around lightning", binsize*1000, adt_signal*1000);


  TCanvas *can3 = new TCanvas("can3","can3",500*1.5,500, 500*1.5, 400*1.5);
  //can->SetLogy();
  gStyle->SetOptStat(1111);
  //hmain->SetMaximum(20);

  hmain->GetXaxis()->SetTitle(text2);
  hmain->GetYaxis()->SetTitle("dt");
  hmain->GetXaxis()->CenterTitle();
  hmain->GetYaxis()->SetTitleOffset(0);
  hmain->GetYaxis()->CenterTitle();
  hmain->GetXaxis()->SetTitleOffset(1.5);
  hmain->Draw("COLZ");

  can3->Print(path + "plots/TH3D_" + title_name + std::to_string(plussminusinsideTGFcandidate)  + ".eps");

  double canvas_size = 400;
  double labelsize = .05;
  double textsizeX = .05;
  double textsizeY = .05;

  double offsetx = 1.;
  double offsety = 0.8;
  double ratioxy = 2;

  TCanvas *can2 = new TCanvas("can2","can2",500,0, canvas_size*ratioxy, canvas_size);
  can2->SetLogy();
  gStyle->SetOptStat(0000);
  //gStyle->SetOptFit(1111);

  projection_poisson->GetXaxis()->SetTitle("Counts per 100 #mus");
  projection_poisson->GetYaxis()->SetTitle("Normalized number of bins");
  ///  projection_poisson->SetTitle(text2);
  projection_poisson->SetTitle("");

  //projection_poisson->SetTitle("#splitline{}{bbb}");

  can2->SetTicks(1,1);
  projection_poisson->GetXaxis()->SetLabelSize(labelsize);
  projection_poisson->GetYaxis()->SetLabelSize(labelsize);
  projection_poisson->GetXaxis()->SetTitleSize(textsizeX);
  projection_poisson->GetYaxis()->SetTitleSize(textsizeY);


  projection_poisson->GetXaxis()->CenterTitle();
  //projection_poisson->GetYaxis()->SetTitleOffset(offsety);
  projection_poisson->GetYaxis()->CenterTitle();
  projection_poisson->GetYaxis()->SetTitleOffset(offsety);
  projection_poisson->SetLineColor(kBlue);



  //
  int number_of_ligtning_ref_period = 9202;
  int number_of_ligtning_ref_background = 9670;

  int number_of_ligtning_drift_perid = 175129; // 2016-2017

  //
  int number_of_ligtning_drift_perid_background = 64559; // 2016-2017
  //int number_of_ligtning_drift_perid_background = 6837; // march - june 2016 + dt
  //int number_of_ligtning_drift_perid_background = 7031; // march - june 2016 - dt

  //int number_of_ligtning_drift_perid = 174287;

  //NORMALIZE ON ENTRIES IN BACKGROUND FILE
  //int number_of_ligtning_ref_period =  2205695;
  //int number_of_ligtning_drift_perid = 5068675;


  //projection_poisson->Scale(1/(plussminusinsideTGFcandidate*2*number_of_ligtning_drift_perid));
  //projection_poisson_background->Scale(1/((0.5)*number_of_ligtning_drift_perid_background));


  projection_poisson->Draw();
  can2->Update();

  double xpos = 8;
  double zpos = 10000;
  double text_size = 0.05;

  double textsize_label = 0.04;
  TLatex lx1,lx2,lx3,lx4, lx5;
  lx1.SetTextSize(textsize_label);
  lx1.SetTextColor(kRed);
  lx1.DrawLatex(14, pow(10,3),"Background ");

  lx2.SetTextSize(textsize_label);
  lx2.SetTextColor(kBlue);
  lx2.DrawLatex(14, pow(10,3)*3,"Signal  ");

  lx3.SetTextSize(textsize_label);
  lx3.SetTextColor(kBlack);
  lx3.SetTextAngle(90);
  lx3.DrawLatex(4.4, 1*pow(10,2),"threshold");

  projection_poisson_background->SetLineColor(kRed);
  projection_poisson_background->Draw("same");

  TLine *line = new TLine(4,0,4,pow(10,4)*2);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->Draw();
  can2->Modified();
  //####################################################
  /*
  lx4.SetTextSize(textsize_label);
  lx4.SetTextColor(1);
  lx4.DrawLatex(14, pow(10,1),"Background REF ");

  lx5.SetTextSize(textsize_label);
  lx5.SetTextColor(7);
  lx5.DrawLatex(14, pow(10,1)*9,"Signal REF ");

  projection_poisson_ref_period->Scale(1/(TMath::Abs(0.0005*2)*(number_of_ligtning_ref_period)));
  projection_poisson_background_ref_period->Scale(1/((0.5)*(number_of_ligtning_ref_background)));

  projection_poisson_ref_period->SetLineColor(7);
  projection_poisson_ref_period->Draw("same");
  projection_poisson_background_ref_period->SetLineColor(1);
  projection_poisson_background_ref_period->Draw("same");

  */
  //####################################################
  can2->Print(path + "plots/poisson_" + title_name  + std::to_string(plussminusinsideTGFcandidate) + "_" + std::to_string(threshold) +  ".eps");

  canvas_size = 400;
  labelsize = .04;
  textsizeX = .05;
  textsizeY = .05;

  offsetx = 1.;
  offsety = 0.6;
  ratioxy = 2;

  TCanvas *can = new TCanvas("can","can",0,500, canvas_size*ratioxy, canvas_size);
  //gStyle->SetOptStat(0000);
  can->SetTicks(1,1);

  char text3[200];
  sprintf(text3,"#deltat [sec]");
  projection_light_curve_signal->GetXaxis()->SetTitle(text3);
  char text1[200];
  sprintf(text1,"Events");
  projection_light_curve_signal->GetYaxis()->SetTitle(text1);
  projection_light_curve_signal->SetTitle("");

  projection_light_curve_signal->GetXaxis()->SetLabelSize(labelsize);
  projection_light_curve_signal->GetYaxis()->SetLabelSize(labelsize);
  projection_light_curve_signal->GetXaxis()->SetTitleSize(textsizeX);
  projection_light_curve_signal->GetYaxis()->SetTitleSize(textsizeY);


  projection_light_curve_signal->GetXaxis()->CenterTitle();
  projection_light_curve_signal->GetYaxis()->SetTitleOffset(offsety);
  projection_light_curve_signal->GetYaxis()->CenterTitle();
  projection_light_curve_signal->GetXaxis()->SetTitleOffset(offsetx);
  projection_light_curve_signal->Draw();

  can->Print(path + "plots/lightcurve_signal_" + title_name + std::to_string(plussminusinsideTGFcandidate) + "_" + std::to_string(threshold) + ".eps");
  //can->Print("/Users/anderslindanger/github/Master/Masteroppgave/images/" + title_name + ".eps");

  /*
  TF1 *ostgaard_fit = new TF1("ostgaard_fit", func_ostgaard_fit ,0,13,4); // x in [0;30], 4 parameters
  ostgaard_fit->SetParName(0,"A");
  ostgaard_fit->SetParName(1,"alpha");
  ostgaard_fit->SetParName(2,"B");
  ostgaard_fit->SetParName(3,"lambda");
  ostgaard_fit->SetParameters(2.03*pow(10,6) ,0.37, 2*pow(10,5), 4);
  //projection_poisson->Fit("ostgaard_fit", "L:R");


  xmin = 0;
  xmax = 7;
  TF1 *poiss = new TF1("poiss", func_poiss ,xmin,xmax);
  poiss->SetParName(0,"A");
  poiss->SetParName(1,"alpha");

  poiss->SetParameters(2*pow(10,7),0.2);
  //hmain->Draw("poiss");

  projection_poisson->Draw();



  */
  // ##############################################

  // FITTING
  TF1 *gauss_background = new TF1("fitFunction", fitFunction,-0.1,0.1,4); // x in [0;30], 4 parameters





  double integral_b = 0;
  integral_hist(projection_poisson_background, threshold, stop, integral_b);
  //char text4[200];
  //sprintf(text4,"Integral background between %d and %d: ", threshold, stop);
  //cout << text4 << integral_b << endl;

  TH1D * projection_poisson_background_threshold = histo_background->ProjectionY("background", threshold, stop);

  //cout << "projection_poisson_background_threshold Mean: "<< projection_poisson_background_threshold->GetMean()<<endl;
  //cout <<"projection_poisson_background_threshold StdDev: "<< projection_poisson_background_threshold->GetStdDev()<<endl;
  //cout << "projection_poisson_background Mean: "<< projection_poisson_background->GetMean()<<endl;
  //cout <<"projection_poisson_background StdDev: "<< projection_poisson_background->GetStdDev()<<endl;


  double integral = 0;
  integral_hist(projection_poisson, threshold, stop, integral);
  //
  char text5[200];
  //
  sprintf(text5,"Integral signal between %d and %d: ", threshold, stop);
  //
  cout << text5 << integral << endl;


  char text6[200];
  sprintf(text6,"N/S: %f ", (integral_b/integral));
  cout << text6 << endl;

  cout << plussminusinsideTGFcandidate*pow(10,3) << "," << threshold << "," << candidates << "," << integral << "," << integral_b << endl;


  // ##############################################

}


void plot_background_DRIFT_REF(TString path, TString title_name, TString name_year, double adt, int threshold, double plussminusinsideTGFcandidate){



  TFile *f_candidate = new TFile(path + "f_candidate_" + name_year + "_" + std::to_string(plussminusinsideTGFcandidate) + "_" + std::to_string(threshold) + ".root","r");
  TTree *tree_candidate = (TTree *) f_candidate->Get("tdata");

  int candidates = tree_candidate->GetEntries();

  // parameters for light curve histogram h1
  double binsize = 0.000100;
  int nbins1 = (int) 2*adt/binsize;
  double adt_signal = adt;
  int nbins1_signal = (int) 2*adt_signal/binsize;

  // parameters for histogram hmain
  int start = 0;
  int stop = 20;
  int nbins2 = (stop-start)/1;

  double adt_background1 = 3.2;
  double adt_background2 = 3.7;
  int nbins_background = 0.5/binsize;


  //TFile *histo = new TFile(path + "stackplots_histogram_" + name_year + "_" + std::to_string(plussminusinsideTGFcandidate) + "_" + std::to_string(threshold) + ".root", "r");
  //TH2D *hmain = (TH2D*) histo->Get("3D_histo");
  //TH1D *projection_poisson = hmain->ProjectionX("projection_poisson_3d_histo");

  //
  TFile *file_hist = new TFile(path + "background_hist" + name_year + ".root", "r");
  //
  TH1D *projection_poisson_background = (TH1D*)file_hist->Get("projection_poisson_background");


    //TFile *file_hist_2 = new TFile(path + "background_plus/background_hist" + name_year + ".root", "r");
    //TH1D *projection_poisson_background_ref_period = (TH1D*)file_hist_2->Get("projection_poisson_background");
  //TH2D *histo_background = (TH2D*)file_hist->Get("3D_histo_background");
  //TH1D * projection_light_curve_backround = histo_background->ProjectionY("lightcurve_over_threshold_background", threshold, stop);


  TString path_2015 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_disabled_analyse/";

  //TFile *histo_ref_period = new TFile(path_2015 + "stackplots_histogram_" + "2015" + "_" + "0.000500" + "_" + "4" + ".root", "r");
  //TH2D *hmain_ref_period = (TH2D*) histo_ref_period->Get("3D_histo");
  //TH1D *projection_poisson_ref_period = hmain_ref_period->ProjectionX("projection_poisson_3d_histo");

  TFile *file_hist_ref_period = new TFile(path_2015 + "background_hist" + "2015" + ".root", "r");
  TH1D *projection_poisson_background_ref_period = (TH1D*)file_hist_ref_period->Get("projection_poisson_background");

  //TFile *file_hist_ref_period_2 = new TFile(path_2015 + "background_pluss/background_hist" + "2015" + ".root", "r");
  //TH1D *projection_poisson_background = (TH1D*)file_hist_ref_period_2->Get("projection_poisson_background");



  char text2[200];
  sprintf(text2,"#particles per %.2f ms time bin summed for each single lightcurve +- %.2f ms around lightning", binsize*1000, adt_signal*1000);

  double canvas_size = 400;
  double labelsize = .05;
  double textsizeX = .05;
  double textsizeY = .05;

  double offsetx = 1.;
  double offsety = 0.8;
  double ratioxy = 2;

  TCanvas *can2 = new TCanvas("can2","can2",500,0, canvas_size*ratioxy, canvas_size);
  can2->SetLogy();
  gStyle->SetOptStat(0000);
  //gStyle->SetOptFit(1111);

  projection_poisson_background->GetXaxis()->SetTitle("Counts per 100 #mus");
  projection_poisson_background->GetYaxis()->SetTitle("Normalized number of bins");
  ///  projection_poisson->SetTitle(text2);
  projection_poisson_background->SetTitle("");

  //projection_poisson->SetTitle("#splitline{}{bbb}");

  can2->SetTicks(1,1);
  projection_poisson_background->GetXaxis()->SetLabelSize(labelsize);
  projection_poisson_background->GetYaxis()->SetLabelSize(labelsize);
  projection_poisson_background->GetXaxis()->SetTitleSize(textsizeX);
  projection_poisson_background->GetYaxis()->SetTitleSize(textsizeY);


  projection_poisson_background->GetXaxis()->CenterTitle();
  //projection_poisson->GetYaxis()->SetTitleOffset(offsety);
  projection_poisson_background->GetYaxis()->CenterTitle();
  projection_poisson_background->GetYaxis()->SetTitleOffset(offsety);
  projection_poisson_background->SetLineColor(kBlue);

  int number_of_ligtning_ref_background = 9670; // - dt
  //  number_of_ligtning_ref_background = 9601; //+ dt


  //int number_of_ligtning_drift_perid_background = 7031; // march - june 2016 - dt

  //
  int number_of_ligtning_drift_perid_background = 6837; // march - june 2016 + dt

  //int number_of_ligtning_drift_perid = 174287;

  //NORMALIZE ON ENTRIES IN BACKGROUND FILE
  //int number_of_ligtning_ref_period =  2205695;
  //int number_of_ligtning_drift_perid = 5068675;


  //projection_poisson->Scale(1/(plussminusinsideTGFcandidate*2*number_of_ligtning_drift_perid));
  //rojection_poisson->GetYaxis()->SetRangeUser(1, pow(10,9));


  //projection_poisson_background->Scale(1/((0.5)*number_of_ligtning_drift_perid_background));


  //projection_poisson->Draw();

  projection_poisson_background->SetLineColor(kRed);
  projection_poisson_background->Draw();
  can2->Update();

  double xpos = 8;
  double zpos = 10000;
  double text_size = 0.05;

  double textsize_label = 0.04;
  TLatex lx1,lx2,lx3,lx4, lx5;
  lx1.SetTextSize(textsize_label);
  lx1.SetTextColor(kRed);
  lx1.DrawLatex(14, pow(10,1)*3,"Background DRIFT");

  lx2.SetTextSize(textsize_label);
  lx2.SetTextColor(kBlue);
  //lx2.DrawLatex(14, pow(10,1)*3*3*3,"Signal DRIFT ");

  lx3.SetTextSize(textsize_label);
  lx3.SetTextColor(kBlack);
  lx3.SetTextAngle(90);
  //lx3.DrawLatex(7.4, 1*pow(10,3),"threshold");



  TLine *line = new TLine(7,0,7,pow(10,4));
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  //line->Draw();
  can2->Modified();
  //####################################################

  lx4.SetTextSize(textsize_label);
  lx4.SetTextColor(6);
  lx4.DrawLatex(14, pow(10,1),"Background REF ");

  lx5.SetTextSize(textsize_label);
  lx5.SetTextColor(7);
  //lx5.DrawLatex(14, pow(10,1)*9,"Signal REF ");

  //projection_poisson_ref_period->Scale(1/(TMath::Abs(0.0005*2*number_of_ligtning_ref_period)));
  //projection_poisson_background_ref_period->Scale(1/((0.5*number_of_ligtning_ref_background)));

  //projection_poisson_ref_period->SetLineColor(7);
  //projection_poisson_ref_period->Draw();
  projection_poisson_background_ref_period->SetLineColor(6);
  //
  //projection_poisson_background_ref_period->Draw("same");



}


int J3_analyse_histogram(){
  //################ < 2015 ############################# DONT USE
  TString path_less_2015 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_enabled_analyse/";
  //create_f_WWLLN_unique_time(path_less_2015, "2008_2015");
  //background_hist(path_less_2015, "", "2015");

  //draw_histogram_and_find_candidates(path_less_2015, "", "2008_2015", 0.1, 5, 0.1);
  //plot(path_less_2015, "", "2008_2015", 0.1, 5, 0.1);





  //################ 2015 #############################
  TString path_2015 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_disabled_analyse/";
  //create_f_WWLLN_unique_time(path_2015, "2015");
  //background_hist(path_2015, "", "2015");

  //draw_histogram_and_find_candidates(path_2015, "reference_periode", "2015", 0.0005, 3, 0.0005);
  //plot(path_2015, "reference_periode", "2015", 0.0005, 3, 0.0005);

  //draw_histogram_and_find_candidates(path_2015, "reference_periode", "2015", 0.000500, 4, 0.000500);
  //plot(path_2015, "reference_periode", "2015", 0.000500, 4, 0.000500);

  //draw_histogram_and_find_candidates(path_2015, "reference_periode", "2015", 0.0005, 5, 0.0005);
  //plot(path_2015, "reference_periode", "2015", 0.0005, 5, 0.0005);



  //draw_histogram_and_find_candidates(path_2015, "reference_periode", "2015", 0.000500, 6, 0.000500);
  //plot(path_2015, "reference_periode", "2015", 0.000500, 6, 0.000500);

  //draw_histogram_and_find_candidates(path_2015, "reference_periode", "2015", 0.000250, 4, 0.000250);
  //plot(path_2015, "reference_periode", "2015", 0.000250, 4, 0.000250);

  //draw_histogram_and_find_candidates(path_2015, "reference_periode", "2015", 0.000250, 5, 0.000250);
  //plot(path_2015, "reference_periode", "2015", 0.000250, 5, 0.000250);

  //draw_histogram_and_find_candidates(path_2015, "reference_periode", "2015", 0.1, 5, 0.1);

  //################ GRID 2015 #############################
  //TString path_2015_GRID = "/Volumes/TOSHIBA/WWLLN_AGILE/GRID/analyse/2015/3901_corrected_data_analyse/try/";
  //create_f_WWLLN_unique_time(path_2015, "2015");
  //draw_histogram_and_find_candidates(path_2015_GRID, "", "2015", 0.0005, 3);

  //draw_histogram_and_find_candidates(path_2015, "reference_periode", "2015", 0.1, 6, 0.1);
  //plot(path_2015, "reference_periode", "2015", 0.1, 6, 0.1);


  //################ Timedrift 2016_2017 #############################
  TString path_2017 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/";
  //create_f_WWLLN_unique_time(path_2017, "2016_2017");
  //background_hist(path_2017, "", "2016_2017");


  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.005, 5, 0.005);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.005, 5, 0.005);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.005, 6, 0.005);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.005, 6, 0.005);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.005, 7, 0.005);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.005, 7, 0.005);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.005, 8, 0.005);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.005, 8, 0.005);




  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 5, 0.01);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 5, 0.01);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 6, 0.01);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 6, 0.01);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 7, 0.01);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 7, 0.01);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 8, 0.01);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 8, 0.01);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 9, 0.01);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.01, 9, 0.01);


  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 7, 0.02);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 5, 0.02);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 5, 0.02);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 6, 0.02);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 6, 0.02);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 7, 0.02);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 7, 0.02);
  //plot_background_DRIFT_REF(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 7, 0.02);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 8, 0.02);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 8, 0.02);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 9, 0.02);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 9, 0.02);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 10, 0.02);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.02, 10, 0.02);



  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 6, 0.05);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 6, 0.05);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 7, 0.05);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 7, 0.05);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 8, 0.05);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 8, 0.05);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 9, 0.05);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 9, 0.05);

  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 10, 0.05);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.05, 10, 0.05);



  //##########################################################
  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.1, 6, 0.1);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.1, 6, 0.1);
  //draw_histogram_and_find_candidates(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.1, 7, 0.1);
  //plot(path_2017, "Timedrift_period_2016_2017", "2016_2017", 0.1, 7, 0.1);
  //##########################################################




  //################ 2018 extra #############################
  TString path_2018 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2018_NEWDRIFT/";
  //create_f_WWLLN_unique_time(path_2018, "2018");
  //background_hist(path_2018, "", "2018");

  //draw_histogram_and_find_candidates(path_2018, "", "2018", 0.5, 6, 0.5);
  //plot(path_2018, "", "2018", 0.5, 6, 0.5);

  //draw_histogram_and_find_candidates(path_2018, "", "2018", 0.02, 6, 0.02);
  //plot(path_2018, "", "2018", 0.02, 6, 0.02);



  //draw_histogram_and_find_candidates(path_2018, "", "2018", 0.02, 4, 0.02);
  //plot(path_2018, "", "2018", 0.02, 4, 0.02);

  //################ 2018 3dfix #############################
  TString path_2018_3dfix = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/3DFIX_4ms/";
  //create_f_WWLLN_unique_time(path_2018_3dfix, "2018");
  //background_hist(path_2018_3dfix, "", "2018");

  //draw_histogram_and_find_candidates(path_2018_3dfix, "", "2018", 0.5, 6, 0.5);
  //plot(path_2018_3dfix, "", "2018", 0.5, 6, 0.5);

  //draw_histogram_and_find_candidates(path_2018_3dfix, "", "2018", 0.02, 6, 0.02);
  //plot(path_2018_3dfix, "", "2018", 0.02, 6, 0.02);

  //################ 2018 3dfix0 #############################
  TString path_2018_3dfix0 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/3DFIX_0ms/";
  //create_f_WWLLN_unique_time(path_2018_3dfix0, "2018");
  //background_hist(path_2018_3dfix0, "", "2018");

  //draw_histogram_and_find_candidates(path_2018_3dfix0, "", "2018", 0.0005, 5, 0.0005);
  //plot(path_2018_3dfix0, "", "2018", 0.0005, 5, 0.0005);
  //draw_histogram_and_find_candidates(path_2018_3dfix0, "", "2018", 0.02, 6, 0.02);

  //draw_histogram_and_find_candidates(path_2018_3dfix0, "", "2018", 0.5, 5, 0.5);




  return 0;
}
