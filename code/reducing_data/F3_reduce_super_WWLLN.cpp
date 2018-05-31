// Written by Anders Lindanger. Master thesis in spacephysics 2018

//Software to find all WWLLN lightning in radius of 1000 km around footprint of reduced AGILE trajectory. WWLLN lightning will be stored in one .root ttree file

// Input: reduced_AGILE_pos.root (Position of satellite only at time of triggered MCAL), super_WWLLN.root (All WWLLN lightning between 23.03-23.06 2015)
// Output: reduced_WWLLN.root which is WWLLN lightning inside +- 0.5 second and less than 1000 km away from reduced_AGILE_pos


static const double pi = M_PI; // pi
static const int R =6378137; // Radius of Earth at equator
static const int c = 299792458; // speed of light
static const int max_distance_footprint = 1000*pow(10,3); // maximum distance in meter between footprint of AGILE and WWLLN lightning
static const int production_altitude_lightning = 15000;

//#################################################################################
// Input: latitude, longitude and height in degrees and meters, and height of lightning
// Returns distance_lightning_AGILE [m] and distance_footprint [m], which is the distance from production point of lightning to AGILE, and WWLLN to AGILE footprint. This is on a sphere (Earth). See "Haversine formula" for more information.
//calculate_distance_lightning_AGILE(h_pos_A, lat_pos_A, lon_pos_A, lat_WWLLN, lon_WWLLN, production_altitude_lightning, distance_footprint, distance_lightning_AGILE);

void calculate_distance_lightning_AGILE(double h_A, double lat_A_deg, double lon_A_deg, double lat_w_deg, double lon_w_deg, double h_L, double &distance_footprint, double &distance_lightning_AGILE)
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


// Inputfiles: reduced_AGILE_pos_2016.root, super_WWLLN_2016_sorted.root
// Outputfile: reduced_super_WWLLN_2016
int reduce_super_WWLLN(TString path, TString path_timedrift, TString path_WWLLN, TString filename, TString name_year, TString name_month){
  clock_t begin = clock();

  // define types
  int year_WWLLN, month_WWLLN, day_WWLLN, hour_WWLLN, minute_WWLLN, sec_WWLLN, usec_WWLLN, year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, entry_pos_AGILE, check, maximum, minimum ;
  float lat_WWLLN, lon_WWLLN, lon_pos_A, lat_pos_A, h_pos_A ;
  double obt_pos_A, sec_in_day_AGILE, sec_in_day_WWLLN, distance_footprint, distance_lightning_AGILE, delta_t, time_WWLLN  ;

  //#################################################################################
  // Create new tree to store data in
  gROOT->Reset();
  // the structure to hold the variables for the branch.
  struct WWLLN {
  };

  WWLLN data;

  // create a new ROOT file
  TFile *f = new TFile(path + path_timedrift + "reduced_super_WWLLN_" + name_year + ".root","RECREATE");
  TTree *tree = new TTree("data","WWLLN along AGILE trajectory");

  // create one branch with all information from the stucture
  tree->Branch("time_WWLLN", &time_WWLLN,"time_WWLLN/D");
  tree->Branch("year", &year_WWLLN,"year/I");
  tree->Branch("month", &month_WWLLN,"month/I");
  tree->Branch("day", &day_WWLLN,"day/I");
  tree->Branch("hour", &hour_WWLLN,"hour/I");
  tree->Branch("minute", &minute_WWLLN,"minute/I");
  tree->Branch("second", &sec_WWLLN,"second/I");
  tree->Branch("usec", &usec_WWLLN,"usec/I");
  tree->Branch("latitude", &lat_WWLLN,"latitude/F");
  tree->Branch("longitude", &lon_WWLLN,"longitude/F");
  //tree->Branch("Entry pos A", &entry_pos_AGILE,"entry_pos_AGILE/I");
  //tree->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");
  //tree->Branch("delta t", &delta_t,"delta t/D");

  //#################################################################################

  //################################################################################
  // Open reduced_AGILE_pos and get year, month, day, hour, min, sec, obt(onboard time), lon, lat, height
  TFile *file_position = new TFile(path + path_timedrift + "reduced_AGILE_pos_" + filename + "_sorted.root", "r");
  TTree *t_pos_A = (TTree*) file_position->Get("data");


  //t_pos_A->SetBranchAddress("year", &year_pos_A);
  //t_pos_A->SetBranchAddress("month",&month_pos_A);
  //t_pos_A->SetBranchAddress("day", &day_pos_A);
  //t_pos_A->SetBranchAddress("hour",&hour_pos_A);
  //t_pos_A->SetBranchAddress("minute", &minute_pos_A);
  //t_pos_A->SetBranchAddress("second", &sec_pos_A);
  t_pos_A->SetBranchAddress("longitude", &lon_pos_A);
  t_pos_A->SetBranchAddress("latitude", &lat_pos_A);
  t_pos_A->SetBranchAddress("height", &h_pos_A);
  t_pos_A->SetBranchAddress("obt_pos", &obt_pos_A);

  int entries_pos_A = t_pos_A->GetEntries();
  cout << entries_pos_A << endl;
  //#################################################################################

  //################################################################################
  // Open super WWLLN and get year, month, day, hour, min, sec, usec, lon, lat
  TFile *WWLLN_file = new TFile(path + path_WWLLN + "super_WWLLN_" + name_year + name_month + "_sorted.root", "r");
  TTree *t_WWLLN = (TTree*)WWLLN_file->Get("WWLLN data");

  //TChain *t_WWLLN = new TChain("WWLLN data");
  //t_WWLLN->Add(path + path_WWLLN + name_year);

  t_WWLLN->SetBranchAddress("time_WWLLN", &time_WWLLN);
  t_WWLLN->SetBranchAddress("year", &year_WWLLN);
  t_WWLLN->SetBranchAddress("month",&month_WWLLN);
  t_WWLLN->SetBranchAddress("day", &day_WWLLN);
  t_WWLLN->SetBranchAddress("hour",&hour_WWLLN);
  t_WWLLN->SetBranchAddress("minute", &minute_WWLLN);
  t_WWLLN->SetBranchAddress("second", &sec_WWLLN);
  t_WWLLN->SetBranchAddress("usec", &usec_WWLLN);
  t_WWLLN->SetBranchAddress("latitude", &lat_WWLLN);
  t_WWLLN->SetBranchAddress("longitude", &lon_WWLLN);
  int entries_WWLLN = t_WWLLN->GetEntries();
  cout << "Building index..." << endl;
  t_WWLLN->BuildIndex("time_WWLLN","usec");

  //TChainIndex *index = new TChainIndex(t_WWLLN,"time_WWLLN", "usec");
  //t_WWLLN->SetTreeIndex(index);
  //#################################################################################

  double limit = 5;
  set<int> event_list;

  cout.precision(15);
  int counter = 0;
  // for every entry in reduced_AGILE_pos
  for (int x= 0 ; x < entries_pos_A ; x++) {
    t_pos_A->GetEntry(x);

    if (counter > entries_pos_A/100) {
      cout << int((x/(float(entries_pos_A + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;

    minimum = t_WWLLN->GetEntryNumberWithBestIndex(obt_pos_A - limit*4);
    maximum = t_WWLLN->GetEntryNumberWithBestIndex(obt_pos_A + limit*4);

    // for every entry in WWLLN files
    for (int z = minimum ; z <= maximum; z++) {
      t_WWLLN->GetEntry(z);

      // If time_AGILE +- limit second = time_WWLLN
      if (obt_pos_A - limit <= time_WWLLN && time_WWLLN <= obt_pos_A + limit) {

        // Calculate distance between lightning and AGILEs footprint
        h_pos_A = h_pos_A*1000;
        calculate_distance_lightning_AGILE(h_pos_A, lat_pos_A, lon_pos_A, lat_WWLLN, lon_WWLLN, production_altitude_lightning, distance_footprint, distance_lightning_AGILE);

        // if distance between lightning and AGILE footpring <= max_distance_footprint -> fill tree
        if (distance_footprint <= max_distance_footprint) {

          // list of entries that keep track of already found WWLLN entries. This is for not getting dublicates.
          if (event_list.count(z) == 0) {
            event_list.insert(z);
            tree->Fill();
          }
        }
      }
    }
  }

  //delete WWLLN_file;
  delete file_position;

  // check what the tree looks like
  //f->ls();
  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries() - 1);

  //tree->Scan("day:hour:minute:second:usec:latitude:longitude");

  //tree->Scan("time_WWLLN","time_WWLLN>357233390-2000","col=20.6f");

  f->Write();

  clock_t end = clock();
  double elapsed_secs = double (end-begin)/ CLOCKS_PER_SEC;
  cout << " " << endl;
  cout << "Program took " << elapsed_secs << " seconds to complete."<< endl;

  return 0;
}

//Inputfile: reduced_super_WWLLN_2016.root
//Outputfile: reduced_super_WWLLN_2016_sorted.root
void sort_WWLLN(TString path, TString name_year){

  // define types
  int year_WWLLN, month_WWLLN, day_WWLLN, hour_WWLLN, minute_WWLLN, sec_WWLLN, usec_WWLLN, year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, entry_pos_AGILE, check, maximum, minimum ;
  float lat_WWLLN, lon_WWLLN, lon_pos_A, lat_pos_A, h_pos_A ;
  double obt_pos_A, sec_in_day_AGILE, sec_in_day_WWLLN, distance_footprint, distance_lightning_AGILE, delta_t, time_WWLLN, propagation_time  ;

  //#################################################################################
  // Create new tree to store data in
  gROOT->Reset();
  // the structure to hold the variables for the branch.
  struct WWLLN {
  };

  WWLLN data;

  // create a new ROOT file
  TFile *f = new TFile(path + "reduced_super_WWLLN_" + name_year + "_sorted.root","RECREATE");

  // create a TTree
  TTree *tree = new TTree("data","WWLLN along AGILE trajectory");

  // create one branch with all information from the stucture
  tree->Branch("time_WWLLN", &time_WWLLN,"time_WWLLN/D");
  tree->Branch("year", &year_WWLLN,"year/I");
  tree->Branch("month", &month_WWLLN,"month/I");
  tree->Branch("day", &day_WWLLN,"day/I");
  tree->Branch("hour", &hour_WWLLN,"hour/I");
  tree->Branch("minute", &minute_WWLLN,"minute/I");
  tree->Branch("second", &sec_WWLLN,"second/I");
  tree->Branch("usec", &usec_WWLLN,"usec/I");
  tree->Branch("latitude", &lat_WWLLN,"latitude/F");
  tree->Branch("longitude", &lon_WWLLN,"longitude/F");
  //tree->Branch("Entry pos A", &entry_pos_AGILE,"entry_pos_AGILE/I");
  //tree->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");
  //tree->Branch("delta t", &delta_t,"delta t/D");

  //#################################################################################

  //################################################################################
  // Open reduced WWLLN  (R  <= 1000 km) and extract year, month, day, hour, min, sec, usec, lon, lat

  TFile *WWLLN_file = new TFile(path + "reduced_super_WWLLN_" + name_year + ".root", "r");
  TTree *t_WWLLN = (TTree*)WWLLN_file->Get("data");

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

  t_WWLLN->BuildIndex("time_WWLLN", "usec");
  int entries_WWLLN = t_WWLLN->GetEntries();
  //#################################################################################

  TTreeIndex *index = (TTreeIndex*) t_WWLLN->GetTreeIndex();

  for (int i = 0; i < index->GetN(); i++){

    Long64_t local = t_WWLLN->LoadTree(index->GetIndex()[i]);

    t_WWLLN->GetEntry(local);
    tree->Fill();

  }

  int counter = 0; double a = 0; double b = 0; tree->SetBranchAddress("time_WWLLN", &time_WWLLN);
  for (int x = 0; x < tree->GetEntries(); x++) {
    tree->GetEntry(x);
    b =time_WWLLN;

    if (a > b) {
      counter += 1;
    }
    a = time_WWLLN;
  }

  //delete WWLLN_file;
  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries() - 1);
  f->Write();
  f->Close();
  cout << "COUNTER: " << counter << endl;

}

//Inputfiles: reduced_super_WWLLN_2016_sorted.root, reduced_AGILE_pos_2016_sorted.root
// Outputfile: reduced_double_super_WWLLN_2016_sorted.root
int remove_double_WWLLN(TString path, TString name_year, TString name_position){

  // define types
  int year_WWLLN, month_WWLLN, day_WWLLN, hour_WWLLN, minute_WWLLN, sec_WWLLN, usec_WWLLN, year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, entry_pos_AGILE, check, maximum, minimum ;
  float lat_WWLLN, lon_WWLLN, lon_pos_A, lat_pos_A, h_pos_A ;
  double obt_pos_A, sec_in_day_AGILE, sec_in_day_WWLLN, distance_footprint, distance_lightning_AGILE, delta_t, time_WWLLN, propagation_time;
  //#################################################################################
  // Create new tree to store data in
  gROOT->Reset();
  // the structure to hold the variables for the branch.
  struct WWLLN {
  };

  WWLLN data;

  // create a new ROOT file
  TFile *f = new TFile(path + "reduced_double_super_WWLLN_" + name_year + "_sorted.root","RECREATE");
  TTree *tree = new TTree("data","WWLLN along AGILE trajectory");
  tree->Branch("time_WWLLN", &time_WWLLN,"time_WWLLN/D");
  tree->Branch("year", &year_WWLLN,"year/I");
  tree->Branch("month", &month_WWLLN,"month/I");
  tree->Branch("day", &day_WWLLN,"day/I");
  tree->Branch("hour", &hour_WWLLN,"hour/I");
  tree->Branch("minute", &minute_WWLLN,"minute/I");
  tree->Branch("second", &sec_WWLLN,"second/I");
  tree->Branch("usec", &usec_WWLLN,"usec/I");
  tree->Branch("latitude", &lat_WWLLN,"latitude/F");
  tree->Branch("longitude", &lon_WWLLN,"longitude/F");
  //tree->Branch("Entry pos A", &entry_pos_AGILE,"entry_pos_AGILE/I");
  //tree->Branch("distance_footprint", &distance_footprint,"distance_footprint/D");
  //tree->Branch("delta t", &delta_t,"delta t/D");

  //#################################################################################

  //################################################################################
  // Open reduced WWLLN  (R  <= 1000 km) and extract year, month, day, hour, min, sec, usec, lon, lat
  TFile *WWLLN_file = new TFile(path + "reduced_super_WWLLN_" + name_year + "_sorted.root", "r");
  TTree *t_WWLLN = (TTree*)WWLLN_file->Get("data");
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
  int entries_WWLLN = t_WWLLN->GetEntriesFast();
  //#################################################################################

  //################################################################################
  // Open reduced_AGILE_pos and get year, month, day, hour, min, sec, obt(onboard time), lon, lat, height
  TFile *file_position = new TFile(path + "reduced_AGILE_pos_" + name_position + "_sorted.root", "r");
  TTree *t_pos_A = (TTree*) file_position->Get("data");

  //t_pos_A->SetBranchAddress("year", &year_pos_A);
  //t_pos_A->SetBranchAddress("month",&month_pos_A);
  //t_pos_A->SetBranchAddress("day", &day_pos_A);
  //t_pos_A->SetBranchAddress("hour",&hour_pos_A);
  //t_pos_A->SetBranchAddress("minute", &minute_pos_A);
  //t_pos_A->SetBranchAddress("second", &sec_pos_A);
  t_pos_A->SetBranchAddress("longitude", &lon_pos_A);
  t_pos_A->SetBranchAddress("latitude", &lat_pos_A);
  t_pos_A->SetBranchAddress("height", &h_pos_A);
  t_pos_A->SetBranchAddress("obt_pos", &obt_pos_A);
  t_pos_A->BuildIndex("obt_pos","0");
  int entries_pos_A = t_pos_A->GetEntriesFast();
  //#################################################################################

  double dt0 = 0, time_WWLLN0 = 0, propagation_time0 = 0;

  for (int i = 0; i < entries_WWLLN; i++) {
    t_WWLLN->GetEntry(i);
    t_pos_A->GetEntryWithIndex(time_WWLLN);


    // Calculate distance between lightning and AGILEs footprint
    h_pos_A = h_pos_A*1000;
    calculate_distance_lightning_AGILE(h_pos_A, lat_pos_A, lon_pos_A, lat_WWLLN, lon_WWLLN, production_altitude_lightning, distance_footprint, distance_lightning_AGILE);
    propagation_time = distance_lightning_AGILE/c;

    if (TMath::Abs((time_WWLLN + propagation_time) - (time_WWLLN0 + propagation_time0)) >= 0.000100 && TMath::Abs((time_WWLLN + propagation_time) - (time_WWLLN0 + propagation_time0)) !=0) { // binsize
      tree->Fill();
    }
    time_WWLLN0 = time_WWLLN;
    propagation_time0 = propagation_time;


  }
  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries() - 1);

  f->Write();

  return 0;
}


int F3_reduce_super_WWLLN(){

  //################ < 2015 #############################
  TString path_less_2015 = "/scratch/Master/fresh_WWLLN_AGILE/";
  //reduce_super_WWLLN("path_less_2015", "AC_enabled_analyse/", "super_WWLLN/", "2008_2015", "2008", "");
  //reduce_super_WWLLN(path_less_2015, "AC_enabled_analyse/", "super_WWLLN/", "2008_2015", "2009", "");
  //reduce_super_WWLLN(path_less_2015, "AC_enabled_analyse/", "super_WWLLN/", "2008_2015", "2010", "");
  //reduce_super_WWLLN(path_less_2015, "AC_enabled_analyse/", "super_WWLLN/", "2008_2015", "2011", "");
  //reduce_super_WWLLN(path_less_2015, "AC_enabled_analyse/", "super_WWLLN/", "2008_2015", "2012", "");
  //reduce_super_WWLLN(path_less_2015, "AC_enabled_analyse/", "super_WWLLN/", "2008_2015", "2013", "");
  //reduce_super_WWLLN(path_less_2015, "AC_enabled_analyse/", "super_WWLLN/", "2008_2015", "2014", "");
  //reduce_super_WWLLN(path_less_2015, "AC_enabled_analyse/", "super_WWLLN/", "2008_2015", "2015", ".01.01");

  //sort_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2008");
  //sort_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2009");
  //sort_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2010");
  //sort_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2011");
  //sort_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2012");
  ///sort_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2013");
  //sort_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2014");
  //sort_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2015");

  //remove_double_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2008", "2008_2015");
  //remove_double_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2009", "2008_2015");
  //remove_double_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2010", "2008_2015");
  //remove_double_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2011", "2008_2015");
  //remove_double_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2012", "2008_2015");
  //remove_double_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2013", "2008_2015");
  //remove_double_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2014", "2008_2015");
  //remove_double_WWLLN(path_less_2015 + "AC_enabled_analyse/", "2015", "2008_2015");

  //################ 2015 #############################
  //reduce_super_WWLLN("/scratch/Master/fresh_WWLLN_AGILE/", "AC_disabled_analyse/", "super_WWLLN/", "2015", "2015", ".03.23_06.23");
  //sort_WWLLN("/scratch/Master/fresh_WWLLN_AGILE/AC_disabled_analyse/", "2015");
  //remove_double_WWLLN("/scratch/Master/fresh_WWLLN_AGILE/AC_disabled_analyse/", "2015", "2015");



  //################ 2016 #############################
  //reduce_super_WWLLN("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift/", "super_WWLLN/", "2016", "2016", "");
  //sort_WWLLN("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/", "2016");
  //remove_double_WWLLN("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/", "2016");


  //################ Timedrift 2017 #############################
  //reduce_super_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift_2017/", "super_WWLLN/", "2016_2017", "2016", "");
  //reduce_super_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift_2017/", "super_WWLLN/", "2016_2017", "2017", ".01.01_12.23");

  //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/", "2016");
  //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/", "2017");
  //remove_double_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/", "2016", "2016_2017");
  //remove_double_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/", "2017", "2016_2017");

  //################ Timedrift 2017 missing months #############################
  //reduce_super_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift_missing_months/", "super_WWLLN/", "2016_2017", "2016", "");
  //reduce_super_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift_missing_months/", "super_WWLLN/", "2016_2017", "2017", "");

  //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_missing_months/", "2016");
  //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_missing_months/", "2017");
  //remove_double_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_missing_months/", "2016", "2016_2017");
  //remove_double_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_missing_months/", "2017", "2016_2017");

  //################ Timedrift 2018 #############################
  //reduce_super_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift_2017_2018_extra/", "super_WWLLN/", "2018", "2018_3.2", "");
  //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017_2018_extra/", "2018");
  //remove_double_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017_2018_extra/", "2018", "2018");


  return 0;
}
