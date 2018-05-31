// Written by Anders Lindanger. Master thesis in spacephysics 2018

// This software open MCAL_time_Etot.root, reduced_super_WWLLN.root and reduced_AGILE_pos.root.
//It creates a file called histo_tree.root which is photons inside 0.5 second from production point,
// and a file calles photon_WWLLN_match which is photon inside 500 us. These new files will be used later to plot histograms.

// Input: MCAL_time_Etot.root, reduced WWLLN.root, reduced_AGILE_pos.root
// Output: histo_tree.root, photon_WWLLN_match


static const double pi = M_PI; // pi
static const int R =6378137; // Radius of Earth at equator
static const int c = 299792458; // speed of light
static const int max_distance_footprint = 1000*pow(10,3); // maximum distance in meter between footprint of AGILE and WWLLN lightning
static const int production_altitude_lightning = 15000;

//#################################################################################
// Input: latitude, longitude and height in degrees and meters, and height of lightning
// Returns distance_lightning_AGILE [m] and distance_footprint [m], which is the distance from production point of lightning to AGILE, and WWLLN to AGILE footprint. This is on a sphere ( Earth). See "Haversine formula" for more information.
//calculate_distance_lightning_AGILE(h_pos_A, lat_pos_A, lon_pos_A, lat_WWLLN, lon_WWLLN, production_altitude_lightning, distance_footprint, distance_lightning_AGILE);
void calculate_distance_lightning_AGILE(double h_A, double lat_A_deg, double lon_A_deg, double lat_w_deg, double lon_w_deg, int h_L, double &distance_footprint, double &distance_lightning_AGILE)
{
  // Degress to radians for WWLLN and AGILE
  double lat_A_rad = (lat_A_deg*pi/180);
  double lon_A_rad = (lon_A_deg*pi/180);
  double lat_w_rad = (lat_w_deg*pi/180);
  double lon_w_rad = (lon_w_deg*pi/180);

  // Calculate distance from lightning production to AGILE using haversine and law of cosine
  double x1 = (sin((lat_A_rad - lat_w_rad)/2));
  double x2 = (sin((lon_A_rad - lon_w_rad)/2));

  double theta = 2*asin(sqrt(pow(x1,2) + cos(lat_w_rad)*cos(lat_A_rad)*pow(x2,2)));
  distance_footprint = R*theta;
  distance_lightning_AGILE = sqrt(pow(R+h_A,2) + pow(R+h_L,2) - 2*(R+h_A)*(R+h_L)*cos(theta));
}
//#################################################################################


int create_TGF_histogram(TString path, TString name_year, TString input_MCAL_file, TString number_outputfile){
  clock_t begin = clock();

  // define types
  int x,y,z;
  int year_WWLLN, month_WWLLN, day_WWLLN, hour_WWLLN, minute_WWLLN, sec_WWLLN, usec_WWLLN, year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, minimum_MCAL, maximum_MCAL, minimum_pos, maximum_pos, contact ;
  float lat_WWLLN, lon_WWLLN, lon_pos_A, lat_pos_A, h_pos_A, Etot_MCAL ;
  double time_pos_A,  time_WWLLN, distance_footprint, distance_lightning_AGILE, obt_MCAL, delta_time_corr, delta_time, h_pos_A_meter, propagation_time, dt, usec, dt_background ;


  //#################################################################################
  // Create new tree to store data in
  gROOT->Reset();
  // the structure to hold the variables for the branch.
  struct TGF {
  };
  TGF data;
  // create a new ROOT file
  TFile *f = new TFile(path + "photon_WWLLN_match_" + name_year + "_" +number_outputfile + ".root","RECREATE");
  TTree *tree = new TTree("data","Photons inside plus minus x microseconds around trigger data MCAL");
  // create one branch with all information from the stucture
  tree->Branch("t_W",&time_WWLLN, "time_WWLLN/D");
  tree->Branch("year_W", &year_WWLLN, "year_WWLLN/I");
  tree->Branch("month_W",&month_WWLLN, "month_WWLLN/I");
  tree->Branch("day_W", &day_WWLLN, "day_WWLLN/I");
  tree->Branch("hour_W",&hour_WWLLN, "hour_WWLLN/I");
  tree->Branch("minute_W", &minute_WWLLN, "minute_WWLLN/I");
  tree->Branch("second_W", &sec_WWLLN, "sec_WWLLN/I");
  tree->Branch("usec_W", &usec_WWLLN, "usec_WWLLN/I");
  tree->Branch("lat_W", &lat_WWLLN,"lat_WWLLN/F");
  tree->Branch("lon_W", &lon_WWLLN,"lon_WWLLN/F");
  tree->Branch("obt_MCAL",&obt_MCAL, "obt_MCAL/D");
  tree->Branch("Etot",&Etot_MCAL, "Etot_MCAL/F");
  tree->Branch("lat_A", &lat_pos_A,"lat_A/F");
  tree->Branch("lon_A", &lon_pos_A,"lon_A/F");
  tree->Branch("height_A", &h_pos_A,"height/F");
  tree->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");
  tree->Branch("distance_lightning_AGILE", &distance_lightning_AGILE,"distance_lightning_AGILE/D");
  tree->Branch("propagation_time", &propagation_time,"propagation_time/D");
  tree->Branch("contact", &contact, "contact/I");

  //#################################################################################

  //#################################################################################
  // Create new tree to store data in
  gROOT->Reset();
  // the structure to hold the variables for the branch.
  struct histo {
  };
  histo hdata;

  // create a new ROOT file
  TFile *f_histo = new TFile(path + "histo_tree_" + name_year + "_" + number_outputfile + ".root","RECREATE");
  TTree *tree_histo = new TTree("hdata","Photons inside plus minus 3 seconds from production time of lightning");

  // create one branch with all information from the stucture
  tree_histo->Branch("obt_MCAL",&obt_MCAL, "obt_MCAL/D");
  tree_histo->Branch("time_WWLLN",&time_WWLLN, "time_WWLLN/D");
  tree_histo->Branch("usec_W", &usec_WWLLN, "usec_WWLLN/I");
  tree_histo->Branch("lat_W", &lat_WWLLN,"lat_WWLLN/F");
  tree_histo->Branch("lon_W", &lon_WWLLN,"lon_WWLLN/F");
  tree_histo->Branch("propagation_time", &propagation_time,"propagation_time/D");
  tree_histo->Branch("Etot_MCAL",&Etot_MCAL, "Etot_MCAL/F");

  //#################################################################################


  //#################################################################################
  // Create new tree to store data in

  // create a new ROOT file
  TFile *f_background = new TFile(path + "background_" + name_year + "_" + number_outputfile + ".root","RECREATE");
  TTree *tree_background = new TTree("bdata","background");
  tree_background->Branch("obt_MCAL",&obt_MCAL, "obt_MCAL/D");
  tree_background->Branch("time_WWLLN",&time_WWLLN, "time_WWLLN/D");
  tree_background->Branch("propagation_time", &propagation_time,"propagation_time/D");


  //#################################################################################

  //################################################################################
  // Open reduced WWLLN  (R  <= 1000 km) and extract year, month, day, hour, min, sec, usec, lon, lat
  //TFile *WWLLN_file = new TFile(path + "reduced_double_super_WWLLN_2016_2017_sorted_folder/" + "*.root", "r");
  //TTree *t_WWLLN = (TTree*)WWLLN_file->Get("data");

  TChain *t_WWLLN = new TChain("data");
  t_WWLLN->Add(path + "reduced_double_super_WWLLN_sorted_folder/" + "*.root");

  t_WWLLN->SetBranchAddress("year", &year_WWLLN);
  t_WWLLN->SetBranchAddress("month",&month_WWLLN);
  t_WWLLN->SetBranchAddress("day", &day_WWLLN);
  t_WWLLN->SetBranchAddress("hour",&hour_WWLLN);
  t_WWLLN->SetBranchAddress("minute", &minute_WWLLN);
  t_WWLLN->SetBranchAddress("second", &sec_WWLLN);
  t_WWLLN->SetBranchAddress("usec", &usec_WWLLN);
  t_WWLLN->SetBranchAddress("latitude", &lat_WWLLN);
  t_WWLLN->SetBranchAddress("longitude", &lon_WWLLN);
  t_WWLLN->SetBranchAddress("time_WWLLN", &time_WWLLN);
  int entries_WWLLN = t_WWLLN->GetEntries();
  //#################################################################################

  //################################################################################
  // Open MCAL time and Etot
  TFile *MCAL_file = new TFile(path + input_MCAL_file, "r");
  TTree *t_MCAL = (TTree*)MCAL_file->Get("data");
  t_MCAL->SetBranchAddress("time", &obt_MCAL);
  t_MCAL->SetBranchAddress("usec", &usec);
  t_MCAL->SetBranchAddress("Etot",&Etot_MCAL);
  t_MCAL->SetBranchAddress("contact", &contact);
  int entries_MCAL = t_MCAL->GetEntries();
  cout << "Building index MCAL time" << endl;
  t_MCAL->BuildIndex("time", "usec");
  cout << "DONE Building index MCAL time" << endl;

  //#################################################################################

  //################################################################################
  // Open reduced_AGILE_pos and extract year, month, day, hour, min, sec, obt(onboard time), lon, lat, height
  TFile *file_position = new TFile(path + "reduced_AGILE_pos_" + name_year + "_sorted.root", "r");
  TTree *t_pos_A = (TTree*)file_position->Get("data");
  //t_pos_A->SetBranchAddress("year", &year_pos_A);
  //t_pos_A->SetBranchAddress("month",&month_pos_A);
  //t_pos_A->SetBranchAddress("day", &day_pos_A);
  //t_pos_A->SetBranchAddress("hour",&hour_pos_A);
  //t_pos_A->SetBranchAddress("minute", &minute_pos_A);
  //t_pos_A->SetBranchAddress("second", &sec_pos_A);
  t_pos_A->SetBranchAddress("obt_pos", &time_pos_A);
  t_pos_A->SetBranchAddress("longitude", &lon_pos_A);
  t_pos_A->SetBranchAddress("latitude", &lat_pos_A);
  t_pos_A->SetBranchAddress("height", &h_pos_A);
  int entries_pos_A = t_pos_A->GetEntries();
  cout << "Building index position time" << endl;
  t_pos_A->BuildIndex("obt_pos","0");
  //#################################################################################

  double limit = 0.1;
  double limit_histo = 3;
  int counter = 0;

  //###############################
  // For each entry in WWLLN file
  for (x = 0; x < entries_WWLLN; x++) {
    t_WWLLN->GetEntry(x);

    if (counter > entries_WWLLN/100) {
      cout << int((x/(float(entries_WWLLN + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;

    minimum_MCAL = t_MCAL->GetEntryNumberWithBestIndex(time_WWLLN - limit_histo*2);
    maximum_MCAL = t_MCAL->GetEntryNumberWithBestIndex(time_WWLLN + limit_histo*2);

    //###############################
    // For each entry in MCAL data
    for (y = minimum_MCAL; y< maximum_MCAL; y++) {
      t_MCAL->GetEntry(y);

      // get position
      t_pos_A->GetEntryWithIndex(obt_MCAL);

      // calculate distance between lightning and AGILE -> Propagation time of photons
      h_pos_A_meter = h_pos_A * 1000;
      calculate_distance_lightning_AGILE(h_pos_A_meter, lat_pos_A, lon_pos_A, lat_WWLLN, lon_WWLLN, production_altitude_lightning, distance_footprint, distance_lightning_AGILE);
      propagation_time = distance_lightning_AGILE/c;
      dt = obt_MCAL - time_WWLLN - propagation_time;

      if (dt >= -limit_histo && dt <= limit_histo && distance_footprint <= max_distance_footprint) {
        tree_histo->Fill();
      }
      // if time_WWLLN + propagation time of photons = time MCAL +- limit seconds -> fill tree
      if (dt >= -limit && dt <= limit && distance_footprint <= max_distance_footprint) {
        tree->Fill();
      }
      if (3.2 <= dt && dt <= 3.8) {
        tree_background->Fill();
      }


    }
  }

  //delete WWLLN_file;
  delete file_position;
  delete MCAL_file;

  // check what the tree looks like
  //f->ls();
  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries()-1);

  //tree_histo->Print();
  //tree_histo->Show(0);
  //tree_background->Print();
  //tree_background->Show(0);
  //tree->Scan("obt_MCAL", "obt_MCAL > 362170532 - 500", "col=20.6f");
  f->Write();
  f_histo->Write();
  f_background->Write();
  clock_t end = clock();
  double elapsed_secs = double (end-begin)/ CLOCKS_PER_SEC;
  cout << " " << endl;
  cout << "Program took " <<elapsed_secs << " seconds to complete."<< endl;


  cout << "Entries MCAL: " <<entries_MCAL << endl;
  return 0;
}

int add_photon_WWLLN_together(TString path, TString name_year){

  // define types
  int x,y,z;
  int year_WWLLN, month_WWLLN, day_WWLLN, hour_WWLLN, minute_WWLLN, sec_WWLLN, usec_WWLLN, year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, minimum_MCAL, maximum_MCAL, minimum_pos, maximum_pos, contact ;
  float lat_WWLLN, lon_WWLLN, lon_pos_A, lat_pos_A, h_pos_A, Etot_MCAL ;
  double time_pos_A,  time_WWLLN, distance_footprint, distance_lightning_AGILE, obt_MCAL, delta_time_corr, delta_time, h_pos_A_meter, propagation_time, dt, usec, dt_background ;


  //#################################################################################

  // create a new ROOT file
  TFile *f_mock = new TFile(path + "photon_mock" + name_year + "_joined_together" + ".root","RECREATE");
  TTree *tree0 = new TTree("data","Photons inside plus minus x microseconds around trigger data MCAL");
  // create one branch with all information from the stucture
  tree0->Branch("t_W",&time_WWLLN, "time_WWLLN/D");
  tree0->Branch("year_W", &year_WWLLN, "year_WWLLN/I");
  tree0->Branch("month_W",&month_WWLLN, "month_WWLLN/I");
  tree0->Branch("day_W", &day_WWLLN, "day_WWLLN/I");
  tree0->Branch("hour_W",&hour_WWLLN, "hour_WWLLN/I");
  tree0->Branch("minute_W", &minute_WWLLN, "minute_WWLLN/I");
  tree0->Branch("second_W", &sec_WWLLN, "sec_WWLLN/I");
  tree0->Branch("usec_W", &usec_WWLLN, "usec_WWLLN/I");
  tree0->Branch("lat_W", &lat_WWLLN,"lat_WWLLN/F");
  tree0->Branch("lon_W", &lon_WWLLN,"lon_WWLLN/F");
  tree0->Branch("obt_MCAL",&obt_MCAL, "obt_MCAL/D");
  tree0->Branch("Etot",&Etot_MCAL, "Etot_MCAL/F");
  tree0->Branch("lat_A", &lat_pos_A,"lat_A/F");
  tree0->Branch("lon_A", &lon_pos_A,"lon_A/F");
  tree0->Branch("height_A", &h_pos_A,"height/F");
  tree0->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");
  tree0->Branch("distance_lightning_AGILE", &distance_lightning_AGILE,"distance_lightning_AGILE/D");
  tree0->Branch("propagation_time", &propagation_time,"propagation_time/D");
  tree0->Branch("contact", &contact, "contact/I");

  //#################################################################################


  //#################################################################################
  //TFile *f_E = new TFile(path + "photon_WWLLN_match_" + name_year + ".root","r");
  //TTree *t_E = (TTree *) f_E->Get("data");
  TChain *t_E = new TChain("data");
  t_E->Add(path + "photon_WWLLN_match_without_duplicates_obt_folder/*.root");
  // create one branch with all information from the stucture
  t_E->SetBranchAddress("t_W",&time_WWLLN);
  t_E->SetBranchAddress("year_W", &year_WWLLN);
  t_E->SetBranchAddress("month_W",&month_WWLLN);
  t_E->SetBranchAddress("day_W", &day_WWLLN);
  t_E->SetBranchAddress("hour_W",&hour_WWLLN);
  t_E->SetBranchAddress("minute_W", &minute_WWLLN);
  t_E->SetBranchAddress("second_W", &sec_WWLLN);
  t_E->SetBranchAddress("usec_W", &usec_WWLLN);
  t_E->SetBranchAddress("lat_W", &lat_WWLLN);
  t_E->SetBranchAddress("lon_W", &lon_WWLLN);
  t_E->SetBranchAddress("obt_MCAL",&obt_MCAL);
  t_E->SetBranchAddress("Etot",&Etot_MCAL);
  t_E->SetBranchAddress("lat_A", &lat_pos_A);
  t_E->SetBranchAddress("lon_A", &lon_pos_A);
  t_E->SetBranchAddress("height_A", &h_pos_A);
  t_E->SetBranchAddress("distance_footprint", &distance_footprint);
  t_E->SetBranchAddress("distance_lightning_AGILE", &distance_lightning_AGILE);
  t_E->SetBranchAddress("propagation_time", &propagation_time);
  t_E->SetBranchAddress("contact", &contact);
  int entries = t_E->GetEntries();
  cout << entries << endl;
  //#################################################################################
  cout << "fill tree0" << endl;
  for (int x = 0; x < entries; x++) {
    t_E->GetEntry(x);
    tree0->Fill();
  }
  cout << "fill tree0 done" << endl;
  cout << "Buildindex" << endl;
  tree0->BuildIndex("t_W","usec_W");
  cout << "Buildindex done" << endl;

  TFile *f = new TFile(path + "photon_WWLLN_match_without_duplicates_" + name_year + "_joined_together" + ".root","RECREATE");

  // create a TTree
  TTree *tree = new TTree("data","Photons inside plus minus x microseconds around trigger data MCAL");

  // create one branch with all information from the stucture
  tree->Branch("t_W",&time_WWLLN, "time_WWLLN/D");
  tree->Branch("year_W", &year_WWLLN, "year_WWLLN/I");
  tree->Branch("month_W",&month_WWLLN, "month_WWLLN/I");
  tree->Branch("day_W", &day_WWLLN, "day_WWLLN/I");
  tree->Branch("hour_W",&hour_WWLLN, "hour_WWLLN/I");
  tree->Branch("minute_W", &minute_WWLLN, "minute_WWLLN/I");
  tree->Branch("second_W", &sec_WWLLN, "sec_WWLLN/I");
  tree->Branch("usec_W", &usec_WWLLN, "usec_WWLLN/I");
  tree->Branch("lat_W", &lat_WWLLN,"lat_WWLLN/F");
  tree->Branch("lon_W", &lon_WWLLN,"lon_WWLLN/F");
  tree->Branch("obt_MCAL",&obt_MCAL, "obt_MCAL/D");
  tree->Branch("Etot",&Etot_MCAL, "Etot_MCAL/F");
  tree->Branch("lat_A", &lat_pos_A,"lat_A/F");
  tree->Branch("lon_A", &lon_pos_A,"lon_A/F");
  tree->Branch("height_A", &h_pos_A,"height/F");
  tree->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");
  tree->Branch("distance_lightning_AGILE", &distance_lightning_AGILE,"distance_lightning_AGILE/D");
  tree->Branch("propagation_time", &propagation_time,"propagation_time/D");
  tree->Branch("contact", &contact, "contact/I");

  cout << "Sorting" << endl;
  TTreeIndex *index = (TTreeIndex*) tree0->GetTreeIndex();
  Long64_t local;
  for (int i = 0; i < index->GetN(); i++){
    local = tree0->LoadTree(index->GetIndex()[i]);
    tree0->GetEntry(local);

    tree->Fill();

  }
  cout << "Sorting done" << endl;

  int counter = 0; double a = 0; double b = 0; tree->SetBranchAddress("t_W", &time_WWLLN);
  for (int x = 0; x < tree->GetEntries(); x++) {
    tree->GetEntry(x);
    b =time_WWLLN;

    if (a > b) {
      counter += 1;
    }
    a = time_WWLLN;
  }



  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries() -1);
  f->Write();
  cout << "COUNTER: " << counter << endl;

  f->Close();
  //t_E->Close();
  return 0;
}

int remove_duplicates_obt(TString path, TString name_year){

  int x,y,z;
  int year_WWLLN, month_WWLLN, day_WWLLN, hour_WWLLN, minute_WWLLN, sec_WWLLN, usec_WWLLN, year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, minimum_MCAL, maximum_MCAL, minimum_pos, maximum_pos, contact ;
  float lat_WWLLN, lon_WWLLN, lon_pos_A, lat_pos_A, h_pos_A, Etot_MCAL ;
  double time_pos_A,  time_WWLLN, distance_footprint, distance_lightning_AGILE, obt_MCAL, delta_time_corr, delta_time, h_pos_A_meter, propagation_time, usec  ;

  //#################################################################################
  // Create new tree to store data in
  gROOT->Reset();

  TFile *f_mock = new TFile(path + "photon_WWLLN_match_" + name_year + "mock.root","RECREATE");

  // create a TTree
  TTree *tree0 = new TTree("data","Photons inside plus minus x microseconds around trigger data MCAL");

  // create one branch with all information from the stucture
  tree0->Branch("t_W",&time_WWLLN, "time_WWLLN/D");
  tree0->Branch("year_W", &year_WWLLN, "year_WWLLN/I");
  tree0->Branch("month_W",&month_WWLLN, "month_WWLLN/I");
  tree0->Branch("day_W", &day_WWLLN, "day_WWLLN/I");
  tree0->Branch("hour_W",&hour_WWLLN, "hour_WWLLN/I");
  tree0->Branch("minute_W", &minute_WWLLN, "minute_WWLLN/I");
  tree0->Branch("second_W", &sec_WWLLN, "sec_WWLLN/I");
  tree0->Branch("usec_W", &usec_WWLLN, "usec_WWLLN/I");
  tree0->Branch("lat_W", &lat_WWLLN,"lat_WWLLN/F");
  tree0->Branch("lon_W", &lon_WWLLN,"lon_WWLLN/F");
  tree0->Branch("obt_MCAL",&obt_MCAL, "obt_MCAL/D");
  tree0->Branch("Etot",&Etot_MCAL, "Etot_MCAL/F");
  tree0->Branch("lat_A", &lat_pos_A,"lat_A/F");
  tree0->Branch("lon_A", &lon_pos_A,"lon_A/F");
  tree0->Branch("height_A", &h_pos_A,"height/F");
  tree0->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");
  tree0->Branch("distance_lightning_AGILE", &distance_lightning_AGILE,"distance_lightning_AGILE/D");
  tree0->Branch("propagation_time", &propagation_time,"propagation_time/D");
  tree0->Branch("contact", &contact, "contact/I");



  //#################################################################################

  //#################################################################################
  TFile *f_E = new TFile(path + "photon_WWLLN_match_" + name_year + ".root","r");
  TTree *t_E = (TTree *) f_E->Get("data");
  // create one branch with all information from the stucture
  t_E->SetBranchAddress("t_W",&time_WWLLN);
  t_E->SetBranchAddress("year_W", &year_WWLLN);
  t_E->SetBranchAddress("month_W",&month_WWLLN);
  t_E->SetBranchAddress("day_W", &day_WWLLN);
  t_E->SetBranchAddress("hour_W",&hour_WWLLN);
  t_E->SetBranchAddress("minute_W", &minute_WWLLN);
  t_E->SetBranchAddress("second_W", &sec_WWLLN);
  t_E->SetBranchAddress("usec_W", &usec_WWLLN);
  t_E->SetBranchAddress("lat_W", &lat_WWLLN);
  t_E->SetBranchAddress("lon_W", &lon_WWLLN);
  t_E->SetBranchAddress("obt_MCAL",&obt_MCAL);
  t_E->SetBranchAddress("Etot",&Etot_MCAL);
  t_E->SetBranchAddress("lat_A", &lat_pos_A);
  t_E->SetBranchAddress("lon_A", &lon_pos_A);
  t_E->SetBranchAddress("height_A", &h_pos_A);
  t_E->SetBranchAddress("distance_footprint", &distance_footprint);
  t_E->SetBranchAddress("distance_lightning_AGILE", &distance_lightning_AGILE);
  t_E->SetBranchAddress("propagation_time", &propagation_time);
  t_E->SetBranchAddress("contact", &contact);
  int entries = t_E->GetEntries();
  cout << "Entries: " << entries << endl;
  //#################################################################################

  double time_WWLLN2, obt_MCAL2, propagation_time2;
  int usec_WWLLN2;
  //#################################################################################
  TFile *f_E2 = new TFile(path + "photon_WWLLN_match_" + name_year + ".root","r");
  TTree *t_E2 = (TTree *) f_E2->Get("data");
  // create one branch with all information from the stucture
  t_E2->SetBranchAddress("t_W",&time_WWLLN2);
  t_E2->SetBranchAddress("usec_W", &usec_WWLLN2);
  t_E2->SetBranchAddress("obt_MCAL",&obt_MCAL2);
  t_E2->SetBranchAddress("propagation_time", &propagation_time2);
  //t_E->->SetBranchAddress("contact", &contact);
  int entries2 = t_E->GetEntries();
  t_E2->BuildIndex("t_W", "usec_W");

  int counter = 0;
  //#################################################################################
  int minimum, maximum, fill_histogram;
  for ( x = 0; x < entries; x++) {

    if (counter > entries/100) {
      cout << int((x/(float(entries + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;

    t_E->GetEntry(x);

    minimum = t_E2->GetEntryNumberWithBestIndex(time_WWLLN - 2);
    maximum = t_E2->GetEntryNumberWithBestIndex(time_WWLLN + 2);
    fill_histogram = 1;
    for ( y = minimum; y < maximum; y++) {
      t_E2->GetEntry(y);

      if (obt_MCAL2 == obt_MCAL && std::fabs(obt_MCAL - propagation_time - time_WWLLN) > std::fabs(obt_MCAL2 - propagation_time2 - time_WWLLN2)  ) {
        fill_histogram = 0;
      }
    }

    if (fill_histogram == 1) {
      tree0->Fill();

    }

  }
  tree0->BuildIndex("t_W", "usec_W");
  // create a new ROOT file
  TFile *f = new TFile(path + "photon_WWLLN_match_" + name_year + "_without_duplicates_obt.root","RECREATE");

  // create a TTree
  TTree *tree = new TTree("data","Photons inside plus minus x microseconds around trigger data MCAL");

  // create one branch with all information from the stucture
  tree->Branch("t_W",&time_WWLLN, "time_WWLLN/D");
  tree->Branch("year_W", &year_WWLLN, "year_WWLLN/I");
  tree->Branch("month_W",&month_WWLLN, "month_WWLLN/I");
  tree->Branch("day_W", &day_WWLLN, "day_WWLLN/I");
  tree->Branch("hour_W",&hour_WWLLN, "hour_WWLLN/I");
  tree->Branch("minute_W", &minute_WWLLN, "minute_WWLLN/I");
  tree->Branch("second_W", &sec_WWLLN, "sec_WWLLN/I");
  tree->Branch("usec_W", &usec_WWLLN, "usec_WWLLN/I");
  tree->Branch("lat_W", &lat_WWLLN,"lat_WWLLN/F");
  tree->Branch("lon_W", &lon_WWLLN,"lon_WWLLN/F");
  tree->Branch("obt_MCAL",&obt_MCAL, "obt_MCAL/D");
  tree->Branch("Etot",&Etot_MCAL, "Etot_MCAL/F");
  tree->Branch("lat_A", &lat_pos_A,"lat_A/F");
  tree->Branch("lon_A", &lon_pos_A,"lon_A/F");
  tree->Branch("height_A", &h_pos_A,"height/F");
  tree->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");
  tree->Branch("distance_lightning_AGILE", &distance_lightning_AGILE,"distance_lightning_AGILE/D");
  tree->Branch("propagation_time", &propagation_time,"propagation_time/D");
  tree->Branch("contact", &contact, "contact/I");


  TTreeIndex *index = (TTreeIndex*) tree0->GetTreeIndex();
  Long64_t local;
  for (int i = 0; i < index->GetN(); i++){
    local = tree0->LoadTree(index->GetIndex()[i]);
    tree0->GetEntry(local);

    tree->Fill();

  }


  //tree->BuildIndex("t_W", "usec_W");
  //tree->Write();
  counter = 0; double a = 0; double b = 0; tree->SetBranchAddress("obt_MCAL", &obt_MCAL);
  for (int x = 0; x < tree->GetEntries(); x++) {
    tree->GetEntry(x);
    b =obt_MCAL;

    if (a > b) {
      counter += 1;
    }
    a = obt_MCAL;
  }

  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries() -1);
  f->Write();
  f_mock->Write();
  cout << "COUNTER: " << counter << endl;

  return 0;
}

int remove_duplicates_obt_background(TString path, TString name_year, TString number_outputfile){

  int x,y,z;
  int year_WWLLN, month_WWLLN, day_WWLLN, hour_WWLLN, minute_WWLLN, sec_WWLLN, usec_WWLLN, year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, minimum_MCAL, maximum_MCAL, minimum_pos, maximum_pos, contact ;
  float lat_WWLLN, lon_WWLLN, lon_pos_A, lat_pos_A, h_pos_A, Etot_MCAL ;
  double time_pos_A,  time_WWLLN, distance_footprint, distance_lightning_AGILE, obt_MCAL, delta_time_corr, delta_time, h_pos_A_meter, propagation_time, usec  ;

  //#################################################################################
  TFile *f = new TFile(path + "background_" + name_year + "_" + number_outputfile + "without_duplicates_obt.root","RECREATE");
  TTree *tree = new TTree("bdata","background");
  tree->Branch("obt_MCAL",&obt_MCAL, "obt_MCAL/D");
  tree->Branch("time_WWLLN",&time_WWLLN, "time_WWLLN/D");
  tree->Branch("propagation_time", &propagation_time,"propagation_time/D");


  //#################################################################################

  //#################################################################################

  TChain *t_E = new TChain("bdata");
  t_E->Add(path + "background/*.root");
  //TFile *f_background = new TFile(path + "background_" + name_year + ".root","r");
  //TTree *t_E = (TTree *) f_background->Get("bdata");
  t_E->SetBranchAddress("obt_MCAL",&obt_MCAL);
  t_E->SetBranchAddress("time_WWLLN",&time_WWLLN);
  t_E->SetBranchAddress("propagation_time", &propagation_time);
  int entries = t_E->GetEntries();
  cout << "Entries: " << entries << endl;
  //#################################################################################

  double time_WWLLN2, obt_MCAL2, propagation_time2;
  int usec_WWLLN2;
  //#################################################################################

  TChain *t_E2 = new TChain("bdata");
  t_E2->Add(path + "background/*.root");
  //TFile *f_background_2 = new TFile(path + "background_" + name_year + ".root","r");
  //TTree *t_E2 = (TTree *) f_background_2->Get("bdata");
  t_E2->SetBranchAddress("obt_MCAL",&obt_MCAL2);
  t_E2->SetBranchAddress("time_WWLLN",&time_WWLLN2);
  t_E2->SetBranchAddress("propagation_time", &propagation_time2);
  int entries2 = t_E2->GetEntries();
  t_E2->BuildIndex("time_WWLLN", "0");

  int counter = 0;
  int i{0};

  set<double> event_list_obt;
  while (t_E->GetEntry(i++)) {

    if (counter > entries/100) {
      cout << int((i/(float(entries + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;


    if (event_list_obt.count(obt_MCAL) == 0) {
      event_list_obt.insert(obt_MCAL);

      tree->Fill();

    }
  }
  //#################################################################################
  /*
  int minimum, maximum, fill_histogram;
  for ( x = 0; x < entries; x++) {

  if (counter > entries/100) {
  cout << int((x/(float(entries + 0.0)))*100) << endl;
  counter = 0;
}
counter += 1;

t_E->GetEntry(x);

minimum = t_E2->GetEntryNumberWithBestIndex(time_WWLLN - 2);
maximum = t_E2->GetEntryNumberWithBestIndex(time_WWLLN + 2);
fill_histogram = 1;
for ( y = minimum; y < maximum; y++) {
t_E2->GetEntry(y);

if (obt_MCAL2 == obt_MCAL && std::fabs(obt_MCAL - propagation_time - time_WWLLN) > std::fabs(obt_MCAL2 - propagation_time2 - time_WWLLN2)  ) {
fill_histogram = 0;
}
}

if (fill_histogram == 1) {
tree->Fill();

}

}
*/

//t_E->Print();
//t_E->Show(0);

//tree->BuildIndex("t_W", "usec_W");
//tree->Write();

tree->Print();
tree->Show(0);
tree->Show(tree->GetEntries() -1);
f->Write();



TFile *f_background_unique_time = new TFile(path + "f_background_unique_time" + name_year + ".root","RECREATE");
TTree *tree_WWLLN = new TTree("tdata","");
tree_WWLLN->Branch("t",&time_WWLLN, "t/D");


i = 0;
counter = 0;
double t_W_before{0};
//Iterate through tree and make a list of all unique WWLLN times.
set<double> event_list;
while (t_E->GetEntry(i++)) {

  if (counter > entries/100) {
    cout << int((i/(float(entries + 0.0)))*100) << endl;
    counter = 0;
  }
  counter += 1;


  if (event_list.count(time_WWLLN) == 0) {
    event_list.insert(time_WWLLN);
    //if ((t_W_before - 0.1) < t_W && t_W < (t_W_before + 0.1)) {
    //cout << t_W_before << "\t" << t_W << "\n"<< endl;
    //}
    tree_WWLLN->Fill();
    t_W_before = time_WWLLN;

  }
}
tree_WWLLN->Print();
//tree_WWLLN->Scan("","","col=20.6f");
//delete f_E;
f_background_unique_time->Write();


return 0;
}



int G3_create_TGF_histogram(){


  //################ < 2015 #############################
  TString path_less_2015 = "/scratch/Master/fresh_WWLLN_AGILE/AC_enabled_analyse/";

  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_2981-6050_usec_sorted.root", "1");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_6051-7999_usec_sorted.root", "2");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_8000-9000_usec_sorted.root", "3");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_9001-10000_usec_sorted.root", "4");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_10001-10700_usec_sorted.root", "5");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_10701-11300_usec_sorted.root", "6");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_11301-12000_1_usec_sorted.root", "7");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_11301-12000_2_usec_sorted.root", "8");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_12001-12500_usec_sorted.root", "9");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_12501-13300_1_usec_sorted.root", "10");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_12501-13300_2_usec_sorted.root", "11"); // zero entries
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_12501-13300_3_usec_sorted.root", "12");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_13301-14000_usec_sorted.root", "13");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_14001-15000_usec_sorted.root", "14");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_15001-16000_usec_sorted.root", "15");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_16001-17000_usec_sorted.root", "16");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_17001-18000_usec_sorted.root", "17");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_18001-19000_usec_sorted.root", "18");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_19001-20000_usec_sorted.root", "19");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_20001-20700_usec_sorted.root", "20");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_20701-21430_usec_sorted.root", "21");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_21431-22200_usec_sorted.root", "22");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_22201-22700_usec_sorted.root", "23");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_22701-23000_usec_sorted.root", "24");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_23001-23400_usec_sorted.root", "25");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_23401-23800_usec_sorted.root", "26");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_23801-24200_usec_sorted.root", "27");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_24201-25400_usec_sorted.root", "28");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_25401-30800_usec_sorted.root", "29");
  //create_TGF_histogram(path_less_2015, "2008_2015", "MCAL_folder/MCAL_time_Etot_2_30801-40900_usec_sorted.root", "30");

  //add_photon_WWLLN_together(path_less_2015, "2008_2015");

  //remove_duplicates_obt(path_less_2015, "2008_2015_joined_together");
  //remove_duplicates_obt_background(path_less_2015, "2008_2015", "");


  //################ 2015 #############################
  TString path_2015 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/AC_disabled_analyse/";
  //create_TGF_histogram(path_2015, "2015", "MCAL_time_Etot_40900_42298_sorted.root", "" );
  //remove_duplicates_obt(path_2015, "2015");
  //remove_duplicates_obt_background(path_2015, "2015", "");

  //################ 2016 #############################
  TString path_2016 = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/";
  //create_TGF_histogram(path_2016, "2016", "MCAL_time_Etot_45604_47597_sorted.root" , "");
  //remove_duplicates_obt(path_2016, "2016");
  //remove_duplicates_obt_background(path_2016, "2016", "");


  //################ Timedrift 2016-2017 #############################
  TString path_2016_2017 = "/scratch/Master/fresh_WWLLN_AGILE/Time_drift_2017/";

  //create_TGF_histogram(path_2016_2017, "2016_2017", "/MCAL_folder/MCAL_time_Etot_45000_49700_sorted.root", "1");
  //remove_duplicates_obt(path_2016_2017, "2016_2017_1");
  //remove_duplicates_obt_background(path_2016_2017, "2016_2017_1", "");


  //create_TGF_histogram(path_2016_2017, "2016_2017", "/MCAL_folder/MCAL_time_Etot_49701_51000_sorted.root", "2");
  //remove_duplicates_obt(path_2016_2017, "2016_2017_2");
  //remove_duplicates_obt_background(path_2016_2017, "2016_2017_2", "");

  /*
  create_TGF_histogram(path_2016_2017, "2016_2017", "/MCAL_folder/MCAL_time_Etot_51001_52600_sorted.root", "3");
  remove_duplicates_obt(path_2016_2017, "2016_2017_3");
  remove_duplicates_obt_background(path_2016_2017, "2016_2017_3", "");

  create_TGF_histogram(path_2016_2017, "2016_2017", "/MCAL_folder/MCAL_time_Etot_52601_53800_sorted.root", "4");
  remove_duplicates_obt(path_2016_2017, "2016_2017_4");
  remove_duplicates_obt_background(path_2016_2017, "2016_2017_4", "");

  create_TGF_histogram(path_2016_2017, "2016_2017", "/MCAL_folder/MCAL_time_Etot_53801_54583_sorted.root", "5");
  remove_duplicates_obt(path_2016_2017, "2016_2017_5");
  remove_duplicates_obt_background(path_2016_2017, "2016_2017_5", "");

  */
  //################ Timedrift 2016-2017 missing months #############################
  TString path_2016_2017_missing_months = "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_missing_months/";

  //create_TGF_histogram(path_2016_2017_missing_months, "2016_2017", "MCAL_folder/MCAL_time_Etot_45000_49700_sorted.root", "1");
  //remove_duplicates_obt(path_2016_2017_missing_months, "2016_2017_1");
  //remove_duplicates_obt_background(path_2016_2017_missing_months, "2016_2017_1", "");


  //create_TGF_histogram(path_2016_2017_missing_months, "2016_2017", "/MCAL_folder/MCAL_time_Etot_49701_51000_sorted.root", "2");
  //remove_duplicates_obt(path_2016_2017_missing_months, "2016_2017_2");
  //remove_duplicates_obt_background(path_2016_2017_missing_months, "2016_2017_2", "");


  //create_TGF_histogram(path_2016_2017_missing_months, "2016_2017", "/MCAL_folder/MCAL_time_Etot_51001_52600_sorted.root", "3");
  //remove_duplicates_obt(path_2016_2017_missing_months, "2016_2017_3");
  //remove_duplicates_obt_background(path_2016_2017_missing_months, "2016_2017_3", "");

  //create_TGF_histogram(path_2016_2017_missing_months, "2016_2017", "/MCAL_folder/MCAL_time_Etot_52601_53800_sorted.root", "4");
  //remove_duplicates_obt(path_2016_2017_missing_months, "2016_2017_4");
  //remove_duplicates_obt_background(path_2016_2017_missing_months, "2016_2017_4", "");

  //create_TGF_histogram(path_2016_2017_missing_months, "2016_2017", "/MCAL_folder/MCAL_time_Etot_53801_54583_sorted.root", "5");
  //remove_duplicates_obt(path_2016_2017_missing_months, "2016_2017_5");
  //remove_duplicates_obt_background(path_2016_2017_missing_months, "2016_2017_5", "");

  add_photon_WWLLN_together(path_2016_2017_missing_months,"2016_2017");

  //################ Timedrift 2018 #############################
  TString path_2018 = "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017_2018_extra/";
  //create_TGF_histogram(path_2018, "2018", "MCAL_folder/MCAL_time_Etot_55655_55935_sorted.root" , "");
  //remove_duplicates_obt(path_2018, "2018");
  //remove_duplicates_obt_background(path_2018, "2018", "");


  return 0;
}
