


double golden_ratio = 1.618;



// This software makes plots.
// Input: .root you want to plot
// Output: three histograms


void plot_histo_2D(TString path, TString inputfile, TString title, TString title_name){

  //#################################################################################
  // Open histo tree

  double propagation_time, obt ,t_W ;
  float energy;
  int nbins_x_energy = 10;


  //######################   PARAMETRE   #########################################

  // Bin length and number of bins of time and energy
  double adt=0.1;
  double binsize = 0.000100;
  int nbins1 = (int) 2*adt/binsize;

  double adt_background1 = 1.5 - adt;
  double adt_background2 = 1.5 + adt;
  int nbins_background = (adt_background2-adt_background1)/binsize;

  double bins[] = {0.1, 0.3, 0.55, 1, 1.7, 3, 5.5, 10, 17, 30,55, 100}; // This is a nice one

  int binnum = sizeof(bins)/sizeof(double) - 1;
  double limit = 0.0001;

  //#################################################################################

  //TFile *f_E = new TFile(path + inputfile,"r");
  //TTree *t_E = (TTree *) f_E->Get("hdata");

  TChain *t_E = new TChain("hdata");
  t_E->Add(path + inputfile);

  t_E->SetBranchAddress("obt_MCAL", &obt);
  t_E->SetBranchAddress("time_WWLLN", &t_W);
  t_E->SetBranchAddress("propagation_time", &propagation_time);
  t_E->SetBranchAddress("Etot_MCAL", &energy);
  int entries = t_E->GetEntries();
  cout << "Entries TChain: " << entries << endl;

  //#################################################################################
  /*
  TFile *histo = new TFile(path + "h_plot_stack.root", "RECREATE");

  TH2F *hcol1 = new TH2F("Data", "Data", nbins1, -adt, adt, binnum, bins);
  TH1F *hcol1_background = new TH1F("Data5", "Data5", nbins_background, adt_background1, adt_background2);



  TH1F *h_zero = new TH1F("Data2", "Data2", nbins_x_energy, 0, 100);
  TH1F *h_background = new TH1F("Data3", "Data3", nbins_x_energy, 0, 100);

  int i=0;

  int counter = 0;
  // scan the tree, fill histo
  while (t_E->GetEntry(i++)) {
    if (counter > entries/100) {
      cout << int((i/(float(entries + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;

    if (-adt <= obt - propagation_time - t_W && obt - propagation_time - t_W <= adt) {
        hcol1->Fill(obt - propagation_time - t_W, energy);
    }

    if (- limit <= obt - propagation_time - t_W && obt - propagation_time - t_W <= limit) {
      h_zero->Fill(energy);
    }

    if (adt_background1 <= obt - propagation_time - t_W && obt - propagation_time - t_W <= adt_background2) {
      hcol1_background->Fill(obt - propagation_time - t_W);
    }


    if (adt_background1 <= obt - propagation_time - t_W && obt - propagation_time - t_W <= adt_background2) {
      h_background->Fill(energy);
    }
  }
  */
  cout << "Opening histograms" << endl;
  //#################################################################################
  TFile *histo = new TFile(path + "h_plot_stack.root", "r");
  TH2F *hcol1 = (TH2F*) histo->Get("Data");

  TH1F *hcol1_background = (TH1F*) histo->Get("Data5");

  TH1F *h_zero = (TH1F*) histo->Get("Data2");
  TH1F *h_background = (TH1F*) histo->Get("Data3");
  // ######################################## GRAPHICS ########################################
  cout << "Opening histograms DONE" << endl;



  // Get mean and sigma of Background

  int binheight, k;
  int nbins = hcol1_background->GetSize()  ; // included Underflow overflow
  double binheight_sum = 0;
  for (k = 1; k < nbins; k++) {
    binheight = hcol1_background->GetBinContent(k);
    binheight_sum += binheight;
  }

  double average_counts_per_bin = binheight_sum/(nbins-2);
  //cout << "nbins: " << nbins << endl;
  //cout << "binheight_sum: " << binheight_sum << endl;
  cout << "average_counts_per_bin: " << average_counts_per_bin << endl;

  double x_x_mean_square_sum = 0;
  for (k = 1; k < nbins; k++) {
    binheight = hcol1_background->GetBinContent(k);
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



  double canvas_size = 400;
  double labelsize = .04;
  double textsize = .05;
  double offsetx = 1.;
  double offsety =0.7;
  double ratioxy = 2;
  // #########################################################
  TCanvas *c1 = new TCanvas("c1","c1",0,0, canvas_size*ratioxy, canvas_size);

  gPad->SetLogy();
  //gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetPalette(55);
  // hcol1->SetMaximum(165);
  hcol1->SetTitle("");

  hcol1->GetXaxis()->SetTitle("Time centered around WWLLN time [sec]");
  hcol1->GetYaxis()->SetTitle("Energy [MeV]");

  hcol1->GetXaxis()->SetLabelSize(labelsize);
  hcol1->GetYaxis()->SetLabelSize(labelsize);
  hcol1->GetXaxis()->SetTitleSize(textsize);
  hcol1->GetYaxis()->SetTitleSize(textsize);

  hcol1->GetXaxis()->CenterTitle();
  hcol1->GetYaxis()->SetTitleOffset(offsety);
  hcol1->GetXaxis()->SetTitleOffset(offsetx);

  hcol1->GetYaxis()->CenterTitle();
  hcol1->Draw("COLZ");
  c1->Print(path + "plots/plot2_" + title_name + "_TGF.jpg");


  // #########################################################




  // #########################################################

  TCanvas *c22 = new TCanvas("c22","c22",0,500, 500*1.5, 400*1.5);
  //gPad->SetLogy();

  //projecton_x->SetFillColor(kBlue+1);
  //Double_t scale3 = 1/(adt_background2-adt_background1);
  //hcol1_background->Scale(scale3);
  //gStyle->SetOptStat(1111);

  hcol1_background->GetXaxis()->SetTitle("Background time [sec]");
  hcol1_background->GetYaxis()->SetTitle("Counts per 100 usec bin");

  hcol1_background->GetXaxis()->CenterTitle();
  hcol1_background->GetYaxis()->SetTitleOffset(1.);
  hcol1_background->GetYaxis()->CenterTitle();
  //projecton_x->SetMaximum(300);
  hcol1_background->Draw();
  //cout << hcol1_background->GetMean(3) << endl;
  // #########################################################

  // #########################################################
  canvas_size = 400;
   labelsize = .04;
   textsize = .05;
   offsetx = 1.;
   offsety = 0.7;
   ratioxy = 2;
  TCanvas *c2 = new TCanvas("c2","c2",0,500, canvas_size*ratioxy, canvas_size);
  TH1D * projecton_x = hcol1 -> ProjectionX();

  //projecton_x->SetFillColor(kBlue+1);
  projecton_x->SetTitle("");
  projecton_x->GetXaxis()->SetTitle("Time centered around WWLLN time [sec]");
  projecton_x->GetYaxis()->SetTitle("Triggers per 100 usec bin");

  projecton_x->GetXaxis()->SetLabelSize(labelsize);
  projecton_x->GetYaxis()->SetLabelSize(labelsize);
  projecton_x->GetXaxis()->SetTitleSize(textsize);
  projecton_x->GetYaxis()->SetTitleSize(textsize);



  projecton_x->GetXaxis()->CenterTitle();
  projecton_x->GetYaxis()->SetTitleOffset(offsety);
  projecton_x->GetXaxis()->SetTitleOffset(offsetx);

  projecton_x->GetYaxis()->CenterTitle();
  //projecton_x->SetMaximum(300);
  projecton_x->Draw();
  //gStyle->SetOptStat(111);

  f_mean->SetLineColor(kRed);
  f_mean->Draw("same");
  f_sigma1->SetLineStyle(7);
  f_sigma1->Draw("same");
  f_sigma3->SetLineStyle(2);
  f_sigma3->Draw("same");
  f_sigma5->SetLineStyle(3);
  f_sigma5->Draw("same");


  c2->Print(path + "plots/plot1_" + title_name + "_TGF.eps");
  //c2->Print("/Users/anderslindanger/github/Master/Masteroppgave/images/" + title_name + ".eps");

  // #########################################################
  // #########################################################
  /*

  // #########################################################
  TCanvas *c3 = new TCanvas("c3","c3",510,0, 500*1.5, 400*1.5);

  //Double_t scale = (0.0002/h_zero->Integral());
  Double_t scale = (1/(limit*2));
  h_zero->Scale(scale);


  //Double_t scale2 = (1/h_background->Integral());
  Double_t scale2 = 1/(adt_background2-adt_background1);
  h_background->Scale(scale2);


  h_zero->SetStats(0);
  h_zero->SetLineColor(kBlack);
  h_zero->GetXaxis()->SetTitle("Energy [MeV]");
  h_zero->GetXaxis()->CenterTitle();
  h_zero->GetYaxis()->SetTitle("Counts per 100 usec bin");
  h_zero->GetYaxis()->SetTitleOffset(1.3);
  h_zero->GetYaxis()->CenterTitle();

  gPad->SetTicks(1,1);
  //gPad->SetLogy();
  //gPad->SetGrid();
  //h_zero->SetMinimum(1);
  //h_zero->SetMaximum(1*pow(10,7));
  h_zero->Draw("h_zero" "E1");
  h_background->SetStats(0);
  h_background->SetLineColor(kRed);
  //
  h_background->Draw("same");


  TLatex latex;
  latex.SetTextSize(0.03);
  latex.DrawLatex(50, 3 * pow(10,6),"Inside +- 500 usec around lightning");
  latex.SetTextColor(kRed);
  latex.DrawLatex(50, 2*pow(10,6),"Background");

  c3->Print("plots/plot3_" + title_name + "_TGF.jpg");
  // #########################################################


  // #########################################################
  TCanvas *c4 = new TCanvas("c4","c4",1000,0, 500*1.5, 400*1.5);

  h_zero->Add(h_background, -1);

  TH1F *h = new TH1F("Data4", title, nbins_x_energy, 0, 100);

  h = h_zero;

  h->SetStats(0);
  h->SetLineColor(kBlack);
  h->GetXaxis()->SetTitle("Energy [MeV]");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("Counts");
  h->GetYaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->CenterTitle();

  gPad->SetTicks(1,1);
  //gPad->SetLogy();
  //gPad->SetGrid();
  //
  h->SetMinimum(0);
  //
  h->SetMaximum(3.3*pow(10,6));
  h->Draw("h_zero" "E1");
  h->SetStats(0);



  latex.SetTextSize(0.03);
  latex.SetTextColor(kBlack);


  latex.DrawLatex(50, 2.8 * pow(10,6),"Background is subtracted");
  //
  c4->Print("plots/plot4_" + title_name + "_TGF.jpg");
  */
  // #########################################################
  //histo->Write();

}

int H3_plot(){

  TString path_less_2015 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_enabled_analyse/";
  //plot_histo_2D(path_less_2015,  "histo_tree/*.root",  "Histo title",  "title_histo");



  TString path_2015 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_disabled_analyse/";
  //plot_histo_2D( path_2015,  "histo_tree/*.root",  "",  "stackplot_histo_2015");


  TString path_2016 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/";
  //plot_histo_2D( path_2016,  "histo_tree/*.root",  "",  "stackplot_histo_2016");

  TString path_2018 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017_2018_extra/";
  //plot_histo_2D( path_2018,  "histo_tree/*.root",  "",  "stackplot_histo_2016");

  TString path_missing_months = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_missing_months/";
  //
  plot_histo_2D( path_missing_months,  "histo_tree/*.root",  "",  "stackplot_histo_2016_2017");

  return 0;
}
