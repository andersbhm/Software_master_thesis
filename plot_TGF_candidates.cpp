int plot(TString save_path, TString inputfile, double t0, int contact, int usec, double propagation_time, double &obt_MCAL_mean)
{
  // MM 14/11/2016 simplified version for Anders: read data file, plot light curve and energy vs. time scatter plot, convert AGILE internal time to UTC

  // reference date for AGILE internal time (need to convert from internal time to UTC)

  TTimeStamp epoch(2004, 1, 1, 0, 0, 0, 0, 1, 0);
  TTimeStamp ts;
  ts.Set((int)  t0, 0, epoch.GetSec(),0);
  ts.SetNanoSec((int) 1.e9*fmod(t0, 1.));

  // open input root data file

  TFile *f = new TFile(inputfile,"r");
  TTree *t = (TTree *) f->Get("tdata");
  double time{0};
  float energy{0};
  t->SetBranchAddress("time", &time);
  t->SetBranchAddress("Etot", &energy);
  double dt = 0.1; // time interval for analysis (and for writing the photon list)
  double adt = 0.1; // time interval for zoom window, default 0.002
  double atbin1 = 0.00005; // short time bin for event histo
  double atbin2 = 0.0001; // long time bin for event histo
  int nbins1 = (int) 2*adt/atbin1;
  int nbins2 = (int) 2*adt/atbin2;
  TH1F *hlc1 = new TH1F("hlc1", "Light curve", nbins1, -adt, adt);
  TH1F *hlc2 = new TH1F("hlc2", "Light curve", nbins2, -adt, adt);
  TH2F *h2e;
  h2e = new TH2F("h2e", "Energy vs. time scatter plot", nbins1, -adt, adt, 2000, 0.1, 2000.);
  int npoint = 0;
  TGraph *ge = new TGraph();

  double obt_MCAL;
  TTree *t_MCAL_mean = new TTree("data","");
  t_MCAL_mean->Branch("obt_MCAL",&time, "obt_MCAL/D");


  // scan the tree, fill histo

  int i=0;
  while (t->GetEntry(i++)) {
    if (time - t0 - propagation_time>=-dt && time - t0 - propagation_time <=dt) {
      //printf("%.6f \t%.2f \t %20.6f \n", time-t0, energy, time);
      if (time - t0 - propagation_time >=-adt && time - t0 - propagation_time <= adt) {
        hlc1->Fill(time - t0 - propagation_time);
        hlc2->Fill(time - t0 - propagation_time);
        ge->SetPoint(npoint++, time - t0 - propagation_time, energy);
        t_MCAL_mean->Fill();
      }
    }
  }


  t_MCAL_mean->SetBranchAddress("obt_MCAL", &obt_MCAL);
  t_MCAL_mean->GetEntry(0);
  double start = obt_MCAL;
  t_MCAL_mean->GetEntry(t_MCAL_mean->GetEntries() -1 );
  double stop = obt_MCAL;

  int nbins3 = TMath::Abs(stop-start)/atbin2;
  TH1F *h_t_MCAL_mean = new TH1F("h_t_MCAL_mean", "time_MCAL", nbins3, start, stop);
  for (int x = 0; x < t_MCAL_mean->GetEntries(); x++) {
    t_MCAL_mean->GetEntry(x);
    h_t_MCAL_mean->Fill(obt_MCAL);
  }

  obt_MCAL_mean = h_t_MCAL_mean->GetXaxis()->GetBinCenter(h_t_MCAL_mean->GetMaximumBin()) - atbin2/2 ;

  //TCanvas *ch1 = new TCanvas("ch1", "", 800, 600);
  //h_t_MCAL_mean->Draw();
  //ch1->Print("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/obt_MCAL.jpg");
  //delete ch1;
  delete t_MCAL_mean;



  // graphics section

  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->Divide(1,2);

  // plot light curve

  c1->cd(1);
  hlc2->SetStats(0);
  hlc2->SetLineColor(kBlue);
  hlc2->GetXaxis()->SetTitle("Time - t0 (s)");
  hlc2->GetXaxis()->SetTitleSize(0.05);
  hlc2->GetXaxis()->CenterTitle();
  hlc2->GetXaxis()->SetLabelSize(0.05);
  hlc2->GetYaxis()->SetTitle("Counts / bin");
  hlc2->GetYaxis()->SetTitleSize(0.05);
  hlc2->GetYaxis()->SetTitleOffset(1.);
  hlc2->GetYaxis()->CenterTitle();
  hlc2->GetYaxis()->SetLabelSize(0.05);
  gPad->SetTicks(1,1);
  hlc2->Draw();

  hlc1->SetStats(0);
  hlc1->SetFillColor(kRed);
  hlc1->SetLineColor(kRed);
  hlc1->Draw("same");

  double xpos = 0.15;
  double zpos = 0.8;
  double text_size = 0.05;

  TLatex t1, t2, t3;   // write TGF time info
  t3.SetTextSize(text_size);
  char str3[100];
  t1.SetTextSize(text_size);
  TString str1("");	// T0 =
  str1 += ts.AsString("s");
  str1 += " UT";
  t1.DrawTextNDC(xpos, zpos-0.1, str1);
  t2.SetTextSize(text_size);
  char str2[200];	// T0 =
  sprintf(str2,"Time_WWLLN %.6f, proptime %.6f contact %d", t0, propagation_time, contact);
  t1.DrawTextNDC(xpos, zpos-0.2, str2);


  // plot energy vs. time plot
  c1->cd(2);
  h2e->Draw();
  h2e->SetStats(0);
  h2e->GetXaxis()->SetTitle("Time - t0 (s)");
  h2e->GetXaxis()->SetTitleSize(0.05);
  h2e->GetXaxis()->CenterTitle();
  h2e->GetXaxis()->SetLabelSize(0.05);
  h2e->GetYaxis()->SetTitle("Energy (MeV)");
  h2e->GetYaxis()->SetTitleSize(0.05);
  h2e->GetYaxis()->SetTitleOffset(1.);
  h2e->GetYaxis()->CenterTitle();
  h2e->GetYaxis()->SetLabelSize(0.05);
  gPad->SetTicks(1,1);
  gPad->SetLogy();
  ge->SetMarkerStyle(5);
  ge->SetMarkerColor(kBlue);
  ge->Draw("Psame");
  // save picture
  TString filename = std::to_string(contact) + "_" + std::to_string(usec);
  c1->Print(save_path + "/plots/Lightcurves/" + filename +".png");
  delete c1;





  delete f;
  return 0;
}

int plot_light_curve_single_TGF_8(TString path_candidate, TString path_MCAL, TString name_year){

  double t_WWLLN{0}, dt{0}, t_W{0}, propagation_time{0}, obt_MCAL_mean{0}, obt_MCAL, distance_footprint{0};
  int contact{0}, bincontent{0}, usec{0}, usec_W{0};
  float  lat_WWLLN{0}, lon_WWLLN{0}, lat_pos_A{0}, lon_pos_A{0}, h_pos_A{0};
  int year_WWLLN{0}, month_WWLLN{0}, day_WWLLN{0}, hour_WWLLN{0}, minute_WWLLN{0}, sec_WWLLN{0};

  TFile *f_candidate = new TFile(path_candidate + "f_candidate_" + name_year + ".root","r");
  TTree *tree_candidate = (TTree *) f_candidate->Get("tdata");
  tree_candidate->SetBranchAddress("t_WWLLN",&t_WWLLN);
  tree_candidate->SetBranchAddress("usec_WWLLN",&usec);
  tree_candidate->SetBranchAddress("year_WWLLN", &year_WWLLN);
  tree_candidate->SetBranchAddress("month_WWLLN",&month_WWLLN);
  tree_candidate->SetBranchAddress("day_WWLLN", &day_WWLLN);
  tree_candidate->SetBranchAddress("hour_WWLLN",&hour_WWLLN);
  tree_candidate->SetBranchAddress("minute_WWLLN", &minute_WWLLN);
  tree_candidate->SetBranchAddress("sec_WWLLN", &sec_WWLLN);
  tree_candidate->SetBranchAddress("dt",&dt);
  tree_candidate->SetBranchAddress("contact", &contact);
  tree_candidate->SetBranchAddress("lat_W", &lat_WWLLN);
  tree_candidate->SetBranchAddress("lon_W", &lon_WWLLN);
  tree_candidate->SetBranchAddress("lat_A", &lat_pos_A);
  tree_candidate->SetBranchAddress("lon_A", &lon_pos_A);
  tree_candidate->SetBranchAddress("height_A", &h_pos_A);
  tree_candidate->SetBranchAddress("bincontent",&bincontent);
  tree_candidate->SetBranchAddress("propagation_time", &propagation_time);
  tree_candidate->SetBranchAddress("distance_footprint", &distance_footprint);

  int entries_candidate = tree_candidate->GetEntries();
  cout << "entries_candidate: " << entries_candidate << endl;
  /*
  //######### ########################################################################
  //TFile *f_E = new TFile(path_candidate + "photon_WWLLN_match_" + name_year + ".root","r");
  TFile *f_E = new TFile(path_candidate + "photon_WWLLN_match_without_duplicates_obt_folder/" + "photon_WWLLN_match_2016_without_duplicates_obt.root","r");
  TTree *t_E = (TTree *) f_E->Get("data");


  //TChain *t_E = new TChain("data");
  //t_E->Add(path_candidate + "photon_WWLLN_match_without_duplicates_obt_folder/" + "*.root");


  //t_E->SetBranchAddress("obt_MCAL", &obt);
  t_E->SetBranchAddress("t_W", &t_W);
  t_E->SetBranchAddress("usec_W", &usec_W);
  t_E->SetBranchAddress("propagation_time", &propagation_time);
  //t_E->SetBranchAddress("Etot", &energy);
  //t_E->SetBranchAddress("lat_W", &lat_WWLLN);
  //t_E->SetBranchAddress("lon_W", &lon_WWLLN);
  //t_E->SetBranchAddress("contact", &contact);
  //t_E->SetBranchAddress("lat_A", &lat_pos_A);
  //t_E->SetBranchAddress("lon_A", &lon_pos_A);
  //t_E->SetBranchAddress("height_A", &h_pos_A);
  int entries = t_E->GetEntries();
  t_E->BuildIndex("t_W", "usec_W");

  //#################################################################################
  */


  TFile *TGF_candidate = new TFile(path_candidate + "TGF_candidates" + name_year + ".root","RECREATE");
  TTree *tree_TGF = new TTree("data","TGF candidates time drift correction period 2016");
  tree_TGF->Branch("t_WWLLN",&t_WWLLN, "t_WWLLN/D");
  tree_TGF->Branch("usec_WWLLN",&usec, "usec_WWLLN/I");
  tree_TGF->Branch("year_WWLLN", &year_WWLLN, "year_WWLLN/I");
  tree_TGF->Branch("month_WWLLN",&month_WWLLN, "month_WWLLN/I");
  tree_TGF->Branch("day_WWLLN", &day_WWLLN, "day_WWLLN/I");
  tree_TGF->Branch("hour_WWLLN",&hour_WWLLN, "hour_WWLLN/I");
  tree_TGF->Branch("minute_WWLLN", &minute_WWLLN, "minute_WWLLN/I");
  tree_TGF->Branch("sec_WWLLN", &sec_WWLLN, "sec_WWLLN/I");
  tree_TGF->Branch("lat_W", &lat_WWLLN, "lat_W/F");
  tree_TGF->Branch("lon_W", &lon_WWLLN, "lon_W/F");
  tree_TGF->Branch("lat_A", &lat_pos_A, "lat_A/F");
  tree_TGF->Branch("lon_A", &lon_pos_A, "lon_A/F");
  tree_TGF->Branch("height_A", &h_pos_A, "height_A/F");
  tree_TGF->Branch("bincontent",&bincontent, "bincontent/I");
  tree_TGF->Branch("obt_MCAL_mean",&obt_MCAL_mean, "obt_MCAL_mean/D");
  tree_TGF->Branch("propagation_time",&propagation_time, "propagation_time/D");
  tree_TGF->Branch("contact", &contact, "contact/I");
  tree_TGF->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");

  // comes from python software with same name
  //int array_candidates[] = {0, 1, 2, 4, 5, 7, 8, 9, 11, 12, 17, 20, 25, 26, 29, 30, 31, 34, 36, 37, 38, 39, 40, 41, 44, 45, 46, 47, 49, 50, 51,
  // 52, 54, 57, 58, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 78, 82, 83, 84};
  /*
  t_WWLLN = 393460618.775256;
  usec = 775256;
  contact = 47394;
  TString contactname = path_MCAL + "RT0" + std::to_string(contact) + "_3908.root";

  t_E->GetEntryWithIndex(t_WWLLN);
  plot(contactname, t_WWLLN, contact, usec, propagation_time);
  */
  for (int i = 0; i < entries_candidate; i++) {
    //for (int i = 0; i < 1; i++) {
    cout << i << " / " << entries_candidate << endl;
    //tree_candidate->GetEntry(array_candidates[i]);
    tree_candidate->GetEntry(i);

    //t_E->GetEntryWithIndex(t_WWLLN, usec);
    //RT045622_3908.root

    TString contactname = path_MCAL + "RT0" + std::to_string(contact) + "_3908.root";
    plot(path_candidate, contactname, t_WWLLN, contact, usec, propagation_time, obt_MCAL_mean);
    tree_TGF->Fill();

    //cout << contactname << endl;
    //cout << "lon_W: " << lon_WWLLN << "\t lat_W: " << lat_WWLLN << "\t propagation_time: " << propagation_time << "\n" << endl;
  }
  tree_TGF->Print();
  tree_TGF->Show(0);
  tree_TGF->Show(tree_TGF->GetEntries() - 1);
  //tree_TGF->Scan("obt_MCAL_mean","","col=20.6f");
  //tree_TGF->Scan("t_WWLLN:propagation_time","","col=20.6f");

  TGF_candidate->Write();

  return 0 ;
}


int K_plot_TGF_candidates(){
  //plot_light_curve_single_TGF_8("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_disabled_analyse/", "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2015/", "2015_0.000500_4");


  //plot_light_curve_single_TGF_8("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/", "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", "2016");
  //plot_light_curve_single_TGF_8("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_disabled_analyse/", "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2015/", "2015");

  //plot_light_curve_single_TGF_8("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/", "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", "2016_2017_0.100000_5");


  //plot_light_curve_single_TGF_8("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/", "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", "2016_2017_0.010000_6");


  //
  plot_light_curve_single_TGF_8("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017_2018_extra/", "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/MCAL_55655_55935/", "2018_0.100000_5");


  return 0;
}
