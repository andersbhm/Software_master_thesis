// Written by Anders Lindanger. Master thesis in spacephysics 2018

void plot_two_hist(TString path_candidate){

  double t_WWLLN{0}, dt{0}, t_W{0}, propagation_time{0}, obt_MCAL_mean{0}, obt_MCAL, distance_footprint{0};
  int contact{0}, bincontent{0}, usec{0}, usec_W{0};
  float  lat_WWLLN{0}, lon_WWLLN{0}, lat_pos_A{0}, lon_pos_A{0}, h_pos_A{0};
  int year_WWLLN{0}, month_WWLLN{0}, day_WWLLN{0}, hour_WWLLN{0}, minute_WWLLN{0}, sec_WWLLN{0};

  TChain *tree_candidate = new TChain("data");
  tree_candidate->Add(path_candidate + "*.root");
  tree_candidate->SetBranchAddress("distance_footprint", &distance_footprint);
  tree_candidate->SetBranchAddress("t_W", &t_WWLLN);

  int entries_candidate = tree_candidate->GetEntries();
  cout << "entries_candidate: " << entries_candidate << endl;


  double distance_footprint2{0};
  TFile *f_candidate2 = new TFile("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/number_of_lightning/total_lightning_along_path_2008-23.06.2015.root" ,"r");
  TTree *tree_candidate2 = (TTree *) f_candidate2->Get("tdata");
  //tree_candidate2->SetBranchAddress("t_WWLLN",&t_WWLLN);
  tree_candidate2->SetBranchAddress("distance_footprint", &distance_footprint2);
  int entries_candidate2 = tree_candidate2->GetEntries();
  cout << "entries_candidate2: " << entries_candidate2 << endl;



  double r0{0}, r_nminus1{0}, r_a{0}, A{0};
  int counter{0}, n{0}, x{0};

  r0 = 200;


  double bins[100];
  counter = 0;
  int stop = 1000;
  for (n = 0; n < stop; n++) {
    r_nminus1 = TMath::Sqrt(n)*r0;

    A = 3.14*(pow(r_nminus1,2) - pow(r_a,2));
    r_a = r_nminus1;
    //cout << A << endl;
    if (r_nminus1 >= stop-100 ) {
      break;
    }
    bins[counter] = r_nminus1;
    //cout << r_nminus1 << endl;
    counter += 1;
  }


//  for (x = 0; x < counter; x++) {
    //cout << bins[x] << endl;
  //}
  TH1F *hist = new TH1F("Data", "", counter - 1, bins);
  TH1F *hist2 = new TH1F("Data2", "", counter - 1, bins);

  set<double> event_list;

  for (x = 0; x < entries_candidate; x++) {
    tree_candidate->GetEntry(x);
    //cout << distance_footprint << endl;
    if (event_list.count(t_WWLLN) == 0) {
      event_list.insert(t_WWLLN);
      hist->Fill(distance_footprint/1000);
    }
  }

  for (x = 0; x < entries_candidate2; x++) {
    tree_candidate2->GetEntry(x);
    //cout << distance_footprint << endl;
    hist2->Fill(distance_footprint2/1000);
  }

  double canvas_size = 400;
  double labelsize = .05;
  double textsize = .045;
  double offsetx = 1.;
  double offsety =0.7;
  double ratioxy = 2;
  // #########################################################


  Double_t scale = 1/hist->Integral("width");
  Double_t scale2 = 1/hist2->Integral("width");

  //h->Scale(scale);

  hist->Scale(scale);
  hist2->Scale(scale2);

  cout << "integral  : " << hist->Integral("width") << " " << hist->Integral() << endl;
  cout << "integral 2: " << hist2->Integral("width") << " " << hist2->Integral() << endl;

  TCanvas *c2 = new TCanvas("c2","c2",100,100,canvas_size*ratioxy, canvas_size);
  //gPad->SetLogy();

  //projecton_x->SetFillColor(kBlue+1);
  c2->SetTicks(1,1);

  gStyle->SetOptStat(0000);
  hist->SetTitle("");

  hist->GetXaxis()->SetTitle("Distance between sub-satellite point and WWLLN detection [km]");
  //
  hist->GetYaxis()->SetTitle("Normalized WWLLN detections ");
  //hist->GetYaxis()->SetTitle("TGF candidates");

  hist->GetXaxis()->SetLabelSize(labelsize);
  hist->GetYaxis()->SetLabelSize(labelsize);
  hist->GetXaxis()->SetTitleSize(textsize);
  hist->GetYaxis()->SetTitleSize(textsize);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleOffset(offsetx);
  hist->GetYaxis()->SetTitleOffset(offsety);

  //projecton_x->SetMaximum(300);
  hist->Draw("e1");
  //hist->Draw("same");
  hist2->SetLineColor(kRed);
  hist2->Draw("same");

  hist->SetTitle("WWLLN with #deltat #leq #pm500");
  hist2->SetTitle("WWLLN detections in TGF field of view");

  c2->BuildLegend();

  //c2->Print("plots/plot1_" + title_name + "_TGF.jpg");
  // #########################################################

}

void plot_hist(TString path_candidate){

  double t_WWLLN{0}, dt{0}, t_W{0}, propagation_time{0}, obt_MCAL_mean{0}, obt_MCAL, distance_footprint{0};
  int contact{0}, bincontent{0}, usec{0}, usec_W{0};
  float  lat_WWLLN{0}, lon_WWLLN{0}, lat_pos_A{0}, lon_pos_A{0}, h_pos_A{0};
  int year_WWLLN{0}, month_WWLLN{0}, day_WWLLN{0}, hour_WWLLN{0}, minute_WWLLN{0}, sec_WWLLN{0};

  TChain *tree_candidate = new TChain("data");
  tree_candidate->Add(path_candidate + "*.root");
  tree_candidate->SetBranchAddress("distance_footprint", &distance_footprint);
  tree_candidate->SetBranchAddress("t_W", &t_WWLLN);

  int entries_candidate = tree_candidate->GetEntries();
  cout << "entries_candidate: " << entries_candidate << endl;


  double r0{0}, r_nminus1{0}, r_a{0}, A{0};
  int counter{0}, n{0}, x{0};

  r0 = 200;


  double bins[100];
  counter = 0;
  int stop = 1000;
  for (n = 0; n < stop; n++) {
    r_nminus1 = TMath::Sqrt(n)*r0;

    A = 3.14*(pow(r_nminus1,2) - pow(r_a,2));
    r_a = r_nminus1;
    //cout << A << endl;
    if (r_nminus1 >= stop-100 ) {
      break;
    }
    bins[counter] = r_nminus1;
    //cout << r_nminus1 << endl;
    counter += 1;
  }


//  for (x = 0; x < counter; x++) {
    //cout << bins[x] << endl;
  //}
  TH1F *hist = new TH1F("Data", "", counter - 1, bins);
  set<double> event_list;

  for (x = 0; x < entries_candidate; x++) {
    tree_candidate->GetEntry(x);
    //cout << distance_footprint << endl;

    if (event_list.count(t_WWLLN) == 0) {
      event_list.insert(t_WWLLN);

      hist->Fill(distance_footprint/1000);

    }
  }
  double canvas_size = 400;
  double labelsize = .05;
  double textsize = .045;
  double offsetx = 1.;
  double offsety =0.7;
  double ratioxy = 2;
  // #########################################################

  TCanvas *c2 = new TCanvas("c2","c2",100,100,canvas_size*ratioxy, canvas_size);
  //gPad->SetLogy();

  //projecton_x->SetFillColor(kBlue+1);
  c2->SetTicks(1,1);

  gStyle->SetOptStat(0000);
  hist->SetTitle("");

  hist->GetXaxis()->SetTitle("Distance between sub-satellite point and WWLLN detection [km]");
  //hist->GetYaxis()->SetTitle("WWLLN detections in TGF field of view");
  hist->GetYaxis()->SetTitle("TGF candidates");

  hist->GetXaxis()->SetLabelSize(labelsize);
  hist->GetYaxis()->SetLabelSize(labelsize);
  hist->GetXaxis()->SetTitleSize(textsize);
  hist->GetYaxis()->SetTitleSize(textsize);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleOffset(offsetx);
  hist->GetYaxis()->SetTitleOffset(offsety);

  //projecton_x->SetMaximum(300);
  hist->Draw("e1");
  hist->Draw("same");

  //c2->Print("plots/plot1_" + title_name + "_TGF.jpg");
  // #########################################################

}




int plot_distance_footprint(){

  //plot_hist("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/GRID/GRID_total/");
  //plot_two_hist("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/GRID/GRID_total/");

  //plot_hist("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_enabled_analyse/photon_WWLLN_match_without_duplicates_obt_folder/");
  //
  plot_two_hist("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_enabled_analyse/photon_WWLLN_match_without_duplicates_obt_folder/");

  return 0;
}
