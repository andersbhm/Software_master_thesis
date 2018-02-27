// This software merge all WWLLN.root files (A20150323.root) to one large super_WWLLN.root file

// Input: Example: A20150323.root to A20150623.root
// Output: super_WWLLN_year.root

int create_super_tree_WWLLN(TString path, TString name_W, TString filename){
  int year_W, month_W, day_W, hour_W, minute_W, sec_W, usec_W ;
  float lat_W_deg, lon_W_deg;
  double time_WWLLN;

  gROOT->Reset();

  // the structure to hold the variables for the branch
  struct WWLLN {
  };
  WWLLN data;

  // create a new ROOT file
  TFile *f = new TFile(path + filename+ ".root","RECREATE");
  TTree *tree = new TTree("WWLLN data","" );

  // create one branch with all information from the stucture
  tree->Branch("time_WWLLN", &time_WWLLN,"time_WWLLN/D");
  tree->Branch("year", &year_W,"year/I");
  tree->Branch("month", &month_W,"month/I");
  tree->Branch("day", &day_W,"day/I");
  tree->Branch("hour", &hour_W,"hour/I");
  tree->Branch("minute", &minute_W,"minute/I");
  tree->Branch("second", &sec_W,"second/I");
  tree->Branch("usec", &usec_W,"usec/I");
  tree->Branch("latitude", &lat_W_deg,"latitude/F");
  tree->Branch("longitude", &lon_W_deg,"longitude/F");

  // Timestamp used by AGILE on board time clock
  TTimeStamp *epoch = new TTimeStamp(2004, 1, 1, 0, 0, 0, 0, 1, 0);

  // For every day in month (WWLLN files, last three digits)
  for (int i{0}; i <= 1300; i++){

    // Finding day in month for WWLLN file. Example: Month = 03. i = 01,02,03,04, ... ,31
    TString day_month_str = Form("%04d", i);
    TString WWLLNtree = name_W + day_month_str + ".root";

    // if file does not exist -> continue
    std::ifstream find(WWLLNtree);
    if (!find) {
      continue;
    }

    cout << WWLLNtree << endl;

    //#################################################################################
    // Open and specify parameters in WWLLN file ...
    TFile *file_W = new TFile(WWLLNtree, "r");
    TTree *WWLLN_data = (TTree*) file_W->Get("WWLLN data");
    WWLLN_data->SetBranchAddress("year", &year_W);
    WWLLN_data->SetBranchAddress("month", &month_W);
    WWLLN_data->SetBranchAddress("day", &day_W);
    WWLLN_data->SetBranchAddress("hour", &hour_W);
    WWLLN_data->SetBranchAddress("minute", &minute_W);
    WWLLN_data->SetBranchAddress("second", &sec_W);
    WWLLN_data->SetBranchAddress("usec", &usec_W);
    WWLLN_data->SetBranchAddress("latitude", &lat_W_deg);
    WWLLN_data->SetBranchAddress("longitude", &lon_W_deg);
    int entries_W = WWLLN_data->GetEntriesFast();
    //#################################################################################

    // For every entry in day
    for (int x{0}; x < entries_W ; x++) {
      WWLLN_data->GetEntry(x);

      // Convert WWLLN time to second since epoch. This is to compare with AGILE onboard time
      TTimeStamp *t = new TTimeStamp(year_W, month_W, day_W, hour_W, minute_W, sec_W, 0,1,0);
      long time_WWLLN_long =  (t->GetSec()  - epoch->GetSec());
      time_WWLLN = time_WWLLN_long + usec_W*pow(10,-6);
      tree->Fill();
    }
    delete file_W;
  }

  // check what the tree looks like
  //tree->Print();
  //tree->Show(0);
  //tree->Scan("year:month:day:hour:minute:second:usec:latitude:longitude");
  f->Write();
  return 0;
}

// This software sorts super WWLLN files into increasing time. Creates a new file named originalname_sorted.root
void sort_WWLLN(TString path_sort, TString filename){
  gROOT->Reset();

  int year_W, month_W, day_W, hour_W, minute_W, sec_W, usec_W ;
  float lat_W_deg, lon_W_deg;
  double time_WWLLN;



  //################################################################################
  TFile *file = new TFile(path_sort + filename + ".root", "r");
  TTree *tdata = (TTree*) file->Get("WWLLN data");
  tdata->SetBranchAddress("time_WWLLN", &time_WWLLN);
  tdata->SetBranchAddress("year", &year_W);
  tdata->SetBranchAddress("month", &month_W);
  tdata->SetBranchAddress("day", &day_W);
  tdata->SetBranchAddress("hour", &hour_W);
  tdata->SetBranchAddress("minute", &minute_W);
  tdata->SetBranchAddress("second", &sec_W);
  tdata->SetBranchAddress("usec", &usec_W);
  tdata->SetBranchAddress("latitude", &lat_W_deg);
  tdata->SetBranchAddress("longitude", &lon_W_deg);
  int entries_tdata = tdata->GetEntriesFast();
  cout << "Building index" << endl;
  tdata->BuildIndex("time_WWLLN");
  cout << "Building index done" << endl;
  //#################################################################################


  //#################################################################################
  // Create new tree to store data in
  TFile *f = new TFile(path_sort + filename +"_sorted.root","RECREATE"); // create a new ROOT file
  TTree *tree = new TTree("WWLLN data","sorted" );   // create a TTree
  tree->Branch("time_WWLLN", &time_WWLLN,"time_WWLLN/D");   // create one branch with all information from the stucture
  tree->Branch("year", &year_W,"year/I");
  tree->Branch("month", &month_W,"month/I");
  tree->Branch("day", &day_W,"day/I");
  tree->Branch("hour", &hour_W,"hour/I");
  tree->Branch("minute", &minute_W,"minute/I");
  tree->Branch("second", &sec_W,"second/I");
  tree->Branch("usec", &usec_W,"usec/I");
  tree->Branch("latitude", &lat_W_deg,"latitude/F");
  tree->Branch("longitude", &lon_W_deg,"longitude/F");
  //#################################################################################


  TTreeIndex *index = (TTreeIndex*) tdata->GetTreeIndex();

  int counter = 0;

  for (int i = 0; i < index->GetN(); i++){
    Long64_t local = tdata->LoadTree(index->GetIndex()[i]);
    tdata->GetEntry(local);

    if (counter > entries_tdata/100) {
      cout << int((i/(float(entries_tdata + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;
    tree->Fill();
  }

  //tree->BuildIndex("time_WWLLN". "usec");
  //tree->Write();
  delete file;
  f->Write();
  f->Close();
}


int B3_create_super_tree_WWLLN(){
  //################ before 2015 #############################
  TString path_2008_15 = "/scratch/Master/fresh_WWLLN_AGILE/super_WWLLN/";
  //sort_WWLLN(path_2008_15, "super_WWLLN_2008");
  //sort_WWLLN(path_2008_15, "super_WWLLN_2009");
  //sort_WWLLN(path_2008_15, "super_WWLLN_2010");
  //sort_WWLLN(path_2008_15, "super_WWLLN_2011");
  //sort_WWLLN(path_2008_15, "super_WWLLN_2012");
  //sort_WWLLN(path_2008_15, "super_WWLLN_2013");
  //sort_WWLLN(path_2008_15, "super_WWLLN_2014");
  //sort_WWLLN(path_2008_15, "super_WWLLN_2015.01.01");





  //################ 2015 #############################
  //sort_WWLLN("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/super_WWLLN/", "super_WWLLN_2015.03.23_06.23");



    //################  2016 #############################

    //create_super_tree_WWLLN("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/", "/scratch/WWLLN_AGILE/WWLLN_data/2016/A2016", "super_WWLLN_2016")
    //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/super_WWLLN/", "super_WWLLN_2016");


  //################ Timedrift 2017 #############################
  //create_super_tree_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/super_WWLLN/", "/media/fer003/TOSHIBA/WWLLN_AGILE/WWLLN_data/2017/A2017", "super_WWLLN_2017.01.01_12.23");
  //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/super_WWLLN/", "super_WWLLN_2017.01.01_12.23");

  //################ Timedrift 2017_2018_extra #############################
  //create_super_tree_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/super_WWLLN/", "/media/fer003/TOSHIBA/WWLLN_AGILE/WWLLN_data/2017/A2017", "super_WWLLN_2017");
  //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/super_WWLLN/", "super_WWLLN_2017");
  //create_super_tree_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/super_WWLLN/", "/media/fer003/TOSHIBA/WWLLN_AGILE/WWLLN_data/2018/A2018", "super_WWLLN_2018_3.2.2018");
  //sort_WWLLN("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/super_WWLLN/", "super_WWLLN_2018_3.2.2018");

  return 0;
}
