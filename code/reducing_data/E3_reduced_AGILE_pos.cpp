// Written by Anders Lindanger. Master thesis in spacephysics 2018

#include <math.h>
#include <stdio.h>      /* printf */
#include <iostream>
#include <ctime>
#include <fstream>

using namespace std;





// This software exlude all data in AGILE_position file which is not +- x second around trigger data (MCAL) from AGILE. AGILE_position.root is AGILEs position for every second for ca 02.2015 - 11.2015.

// Input: reduce_RT_MCAL.root and AGILE_position.root
// Output: reduced_AGILE_pos.root
int reduced_AGILE_pos(TString path, TString path_working, TString path_AGILE, TString output_name_year, TString input_AGILE_pos_file, TString input_MCAL_file ){
  clock_t begin = clock();

  // define types
  int year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, name_number, x, minimum, maximum;
  float lon_pos_A, lat_pos_A, h_pos_A, Etot_MCAL;
  double obt_pos_A, obt_MCAL, t_start, t_end, last_t, first_t, check, MCAL_time;

  //#################################################################################
  // Create new tree to store data in
  gROOT->Reset();

  // the structure to hold the variables for the branch.
  struct MCAL {
  };

  MCAL data;

  // create a new ROOT file
  TFile *f = new TFile(path + path_working + "reduced_AGILE_pos_" + output_name_year + ".root","RECREATE");
  TTree *tree = new TTree("data","AGILE position and time with accuracy one second. Plus minus one second around trigger data MCAL");
  tree->Branch("year", &year_pos_A,"year/I");
  tree->Branch("month", &month_pos_A,"month/I");
  tree->Branch("day", &day_pos_A,"day/I");
  tree->Branch("hour", &hour_pos_A,"hour/I");
  tree->Branch("minute", &minute_pos_A,"minute/I");
  tree->Branch("second", &sec_pos_A,"second/I");
  tree->Branch("latitude", &lat_pos_A,"latitude/F");
  tree->Branch("longitude", &lon_pos_A,"longitude/F");
  tree->Branch("height", &h_pos_A,"height/F");
  tree->Branch("obt_pos", &obt_pos_A,"obt_pos/D");

  //#################################################################################

  //################################################################################
  // Open AGILE_position.root and get year, month, day, hour, min, sec, obt(onboard time), lon, lat, h
  //TFile *file_position = new TFile(path + path_AGILE + input_AGILE_pos_file, "r");
  //TTree *position = (TTree*)file_position->Get("position");

  TChain *position = new TChain("position");
  position->Add(path + path_AGILE + input_AGILE_pos_file);

  position->SetBranchAddress("year", &year_pos_A);
  position->SetBranchAddress("month",&month_pos_A);
  position->SetBranchAddress("day", &day_pos_A);
  position->SetBranchAddress("hour",&hour_pos_A);
  position->SetBranchAddress("min", &minute_pos_A);
  position->SetBranchAddress("sec", &sec_pos_A);
  position->SetBranchAddress("obt", &obt_pos_A);
  position->SetBranchAddress("lon", &lon_pos_A);
  position->SetBranchAddress("lat", &lat_pos_A);
  position->SetBranchAddress("h", &h_pos_A);
  int entries_pos_A = position->GetEntries();
  position->BuildIndex("obt","0");
  //#################################################################################

  //################################################################################
  // Open MCAL data ( RT04....._3908.root) and extract time (onboard time)
  //TFile *file_MCAL = new TFile(path + path_working + input_MCAL_file, "r");
  //TTree *tdata = (TTree*) file_MCAL->Get("data");

  TChain *tdata = new TChain("data");
  tdata->Add(path + path_working + input_MCAL_file);

  tdata->SetBranchAddress("time", &MCAL_time);
  int entries_data = tdata->GetEntries();
  //#################################################################################

  double limit = 6; // Limit for +- one seconds around trigger data.
  set <int> event_list;
  int c = 0;

  // For every line in MCAL reduced file
  for (int z = 0; z < entries_data; z++) {
    tdata->GetEntry(z);

    if (c > entries_data/150) {
      cout << int((z/(entries_data + 0.0))*100) << endl;
      c = 0;
    }
    c += 1;

    minimum = position->GetEntryNumberWithBestIndex(MCAL_time - limit*2);
    maximum = position->GetEntryNumberWithBestIndex(MCAL_time + limit*2);

    // For every second in AGILE_position
    for (int x = minimum ; x <= maximum; x++) {
      position->GetEntry(x);

      // Check if obt_pos_A - limit <= MCAL_time < obt_post_A + limit .
      if (MCAL_time - limit <= obt_pos_A && obt_pos_A <= MCAL_time + limit){

        // list of entries that keep track of already found WWLLN entries. This is for not getting dublicates.
        if (event_list.count(x) == 0) {
          event_list.insert(x);
          tree->Fill();
        }
      }
    }
  }


  //delete file_MCAL;
  //delete file_position;

  // check what the tree looks like
  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries() - 1);

  //tree->Scan("obt_pos","obt_pos>362170532-2000","col=20.6f");
  //tree->Scan("obt_pos:year:month:day:hour:minute:second:obt_pos");
  f->Write();

  clock_t end = clock();
  double elapsed_secs = double (end-begin)/ CLOCKS_PER_SEC;
  cout << "Program took " <<elapsed_secs << " seconds to complete."<< endl;

  return 0;
}



int Add_AGILE_positions_together(){

  int year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, name_number, x, minimum, maximum;
  float lon_pos_A, lat_pos_A, h_pos_A, Etot_MCAL;
  double obt_pos_A, obt_MCAL, t_start, t_end, last_t, first_t, check, MCAL_time;
  //#################################################################################
  // Create new tree to store data in

  TString path = "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/AGILE_position/";
  gROOT->Reset();

  // the structure to hold the variables for the branch.
  struct MCAL {
  };

  MCAL data;

  // create a new ROOT file
  TFile *f = new TFile(path + "AGILE_position_missing_months.root","RECREATE");
  TTree *tree = new TTree("position","AGILE position and time with accuracy one second. Plus minus one second around trigger data MCAL");
  tree->Branch("year", &year_pos_A,"year/I");
  tree->Branch("month", &month_pos_A,"month/I");
  tree->Branch("day", &day_pos_A,"day/I");
  tree->Branch("hour", &hour_pos_A,"hour/I");
  tree->Branch("minute", &minute_pos_A,"minute/I");
  tree->Branch("second", &sec_pos_A,"second/I");
  tree->Branch("latitude", &lat_pos_A,"latitude/F");
  tree->Branch("longitude", &lon_pos_A,"longitude/F");
  tree->Branch("height", &h_pos_A,"height/F");
  tree->Branch("obt_pos", &obt_pos_A,"obt_pos/D");

  //#################################################################################



  //################################################################################
  // Open AGILE_position.root and get year, month, day, hour, min, sec, obt(onboard time), lon, lat, h


  TChain *position = new TChain("position");
  position->Add(path + "working_folder/*.root");
  //position->Add(path + "AGILE_position_2017.01.01_11.09.root");


  position->SetBranchAddress("year", &year_pos_A);
  position->SetBranchAddress("month",&month_pos_A);
  position->SetBranchAddress("day", &day_pos_A);
  position->SetBranchAddress("hour",&hour_pos_A);
  position->SetBranchAddress("min", &minute_pos_A);
  position->SetBranchAddress("sec", &sec_pos_A);
  position->SetBranchAddress("obt", &obt_pos_A);
  position->SetBranchAddress("lon", &lon_pos_A);
  position->SetBranchAddress("lat", &lat_pos_A);
  position->SetBranchAddress("h", &h_pos_A);
  int entries_pos_A = position->GetEntries();
  //#################################################################################

  for (int i = 0; i < entries_pos_A; i++) {
    position->GetEntry(i);
    tree->Fill();
  }

  f->Write();

  tree->Print();
  tree->Show(0);
  tree->Show(entries_pos_A-1);

  return 0;
}





// This software sorts AGILE files into increasing time.
void sort_AGILE(TString path, TString name_year){


  int year_pos_A, month_pos_A, day_pos_A, hour_pos_A, minute_pos_A, sec_pos_A, name_number, x, minimum, maximum;
  float lon_pos_A, lat_pos_A, h_pos_A, Etot_MCAL;
  double obt_pos_A, obt_MCAL, t_start, t_end, last_t, first_t, check, MCAL_time;


  //#################################################################################
  // Create new tree to store data in
  gROOT->Reset();

  // the structure to hold the variables for the branch.
  struct MCAL {
  };

  MCAL data;

  // create a new ROOT file

  TFile *f = new TFile(path + "reduced_AGILE_pos_" + name_year + "_sorted.root","RECREATE");

  // create a TTree
  TTree *tree = new TTree("data","AGILE position and time with accuracy one second. Plus minus one second around trigger data MCAL");

  // create one branch with all information from the stucture
  tree->Branch("year", &year_pos_A,"year/I");
  tree->Branch("month", &month_pos_A,"month/I");
  tree->Branch("day", &day_pos_A,"day/I");
  tree->Branch("hour", &hour_pos_A,"hour/I");
  tree->Branch("minute", &minute_pos_A,"minute/I");
  tree->Branch("second", &sec_pos_A,"second/I");
  tree->Branch("latitude", &lat_pos_A,"latitude/F");
  tree->Branch("longitude", &lon_pos_A,"longitude/F");
  tree->Branch("height", &h_pos_A,"height/F");
  tree->Branch("obt_pos", &obt_pos_A,"obt_pos/D");

  //#################################################################################

  //################################################################################
  // Open reduced_AGILE_pos and get year, month, day, hour, min, sec, obt(onboard time), lon, lat, height

  TFile *file_position = new TFile(path + "reduced_AGILE_pos_" + name_year + ".root", "r");
  TTree *t_pos_A = (TTree*) file_position->Get("data");
  //TChain *t_pos_A = new TChain("data");
  //t_pos_A->Add(path + "working_folder/*.root");

  t_pos_A->SetBranchAddress("year", &year_pos_A);
  t_pos_A->SetBranchAddress("month",&month_pos_A);
  t_pos_A->SetBranchAddress("day", &day_pos_A);
  t_pos_A->SetBranchAddress("hour",&hour_pos_A);
  t_pos_A->SetBranchAddress("minute", &minute_pos_A);
  t_pos_A->SetBranchAddress("second", &sec_pos_A);
  t_pos_A->SetBranchAddress("longitude", &lon_pos_A);
  t_pos_A->SetBranchAddress("latitude", &lat_pos_A);
  t_pos_A->SetBranchAddress("height", &h_pos_A);
  t_pos_A->SetBranchAddress("obt_pos", &obt_pos_A);
  t_pos_A->BuildIndex("obt_pos", "0");
  int entries_pos_A = t_pos_A->GetEntries();
  cout << entries_pos_A << endl;
  //#################################################################################
  TTreeIndex *index = (TTreeIndex*) t_pos_A->GetTreeIndex();
  Long64_t local;
  for (int i = 0; i < index->GetN(); i++){
    local = t_pos_A->LoadTree(index->GetIndex()[i]);
    t_pos_A->GetEntry(local);

    tree->Fill();

  }

  int counter = 0; double a = 0; double b = 0; tree->SetBranchAddress("obt_pos", &obt_pos_A);
  for (int x = 0; x < tree->GetEntries(); x++) {
    tree->GetEntry(x);
    b =obt_pos_A;

    if (a > b) {
      counter += 1;
    }
    a = obt_pos_A;
  }

  //tree->BuildIndex("obt_pos", "0");
  //tree->Write();


  //delete file_position;
  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries() - 1);
  f->Write();
  f->Close();

  cout << "COUNTER: " << counter << endl;
}



int E3_reduced_AGILE_pos(){


  //################ before 2015 #############################
  //reduced_AGILE_pos("/scratch/Master/fresh_WWLLN_AGILE/", "AC_enabled_analyse/", "AGILE_position/", "2015", "AGILE_position_2015.03.22_11.05.root", "MCAL_time_Etot_40900_42298_sorted.root" );
  //Add_AGILE_positions_together();

  //sort_AGILE("/scratch/Master/fresh_WWLLN_AGILE/AC_enabled_analyse/reduced_AGILE_pos/", "2008_2015");



  //################ 2015 #############################
  //reduced_AGILE_pos("/scratch/Master/fresh_WWLLN_AGILE/", "AC_disabled_analyse/", "AGILE_position/", "2015", "AGILE_position_2015.03.22_11.05.root", "MCAL_time_Etot_40900_42298_sorted.root" );
  //sort_AGILE("/scratch/Master/fresh_WWLLN_AGILE/AC_disabled_analyse/", "2015");

  //################ 2016 #############################
  //reduced_AGILE_pos("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift/", "AGILE_position/", "2016", "AGILE_position_2016.01.01_12.16.root", "MCAL_time_Etot_45604_47597_sorted.root" );
  //sort_AGILE("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/", "2016");

  //################ Timedrift 2016-2017 #############################
  //Add_AGILE_positions_together();
  //reduced_AGILE_pos("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift_2017/", "AGILE_position/working_folder/", "2016_2017", "*.root", "MCAL_folder/*.root" );
  //sort_AGILE("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/", "2016_2017");


    //################ Timedrift 2016-2017 missing months #############################
    //Add_AGILE_positions_together();
    //reduced_AGILE_pos("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift_missing_months/", "AGILE_position/working_folder/", "2016_2017", "*.root", "MCAL_folder/*.root" );
    //
    sort_AGILE("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_missing_months/", "2016_2017");

    //################ Timedrift 2018 #############################
    //Add_AGILE_positions_together();
    //reduced_AGILE_pos("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/", "Time_drift_2017_2018_extra/", "AGILE_position/working_folder/", "2017_2018", "*.root", "MCAL_folder/*.root" );
    //sort_AGILE("/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017_2018_extra/", "2017_2018_extra");




  return 0;
}
