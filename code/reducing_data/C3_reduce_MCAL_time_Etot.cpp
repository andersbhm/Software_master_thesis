// Written by Anders Lindanger. Master thesis in spacephysics 2018

using namespace std;
static const double pi = M_PI; // pi
//#include <math.h>       /* modf */

// This software join MCAL_trigger data to one file, keeps only onboard time and energy (MeV) of photons.

// Input: example: RT040900_3908.root to RT042298_3908.root
// Output: MCAL_time_Etot.root
int reduce_MCAL_time_Etot(TString path, TString path_MCAL, int file_MCAL_name_start_int, int file_MCAL_name_stop_int){

  gROOT->Reset();

  double obt_MCAL;
  float Etot_MCAL;
  int entries_tdata, file_MCAL_name_int, contact;
  double usec;
  //#################################################################################
  // Create new tree to store data in

  // the structure to hold the variables for the branch
  struct MCAL {
  };

  MCAL data;

  // create a new ROOT file
  TFile *f = new TFile(path + "MCAL_time_Etot_" + std::to_string(file_MCAL_name_start_int) + "_" + std::to_string(file_MCAL_name_stop_int) + ".root", "RECREATE");
  // create a TTree
  TTree *tree = new TTree("data","");

  // create one branch with all information from the stucture
  tree->Branch("time", &obt_MCAL,"time/D");
  tree->Branch("usec", &usec,"usec/D");
  tree->Branch("Etot", &Etot_MCAL,"Etot/F");
  tree->Branch("contact", &contact, "contact/I");
  //#################################################################################

  // For every file
  for (int y{file_MCAL_name_start_int}; y <= file_MCAL_name_stop_int; y++) {
    //################################################################################
    // Open MCAL data ( RT04....._3908.root) and get time (onboard time)
    TString name = "RT";
    TString file_MCAL_name_str = Form("%06d", y); // Format y = 1 to y = 000001
    TString inputfile = name + file_MCAL_name_str + "_3908.root";
    contact = y;
    // if file does not exist -> continue
    std::ifstream find(path_MCAL + inputfile);
    if (!find) {
      continue;
    }

    //file_MCAL_name_int = y;

    cout << inputfile<< endl;

    TFile *file_MCAL = new TFile(path_MCAL + inputfile, "r");

    TTree *tdata = (TTree*) file_MCAL->Get("tdata");
    tdata->SetBranchAddress("time", &obt_MCAL);
    tdata->SetBranchAddress("Etot", &Etot_MCAL);

    entries_tdata = tdata->GetEntries();
    //#################################################################################

    // For every entry in file -> fill tree
    for (int x{0}; x < entries_tdata; x++) {
      tdata->GetEntry(x);
      usec = round((obt_MCAL - int(obt_MCAL))*1000000);
      tree->Fill();
    }

    delete file_MCAL;
  }

  // check what the tree looks like
  //f->ls();
  tree->Print();
  tree->Show(0);
  tree->Show(tree->GetEntries() - 1);
  //tree->Scan("time:usec", "", "col=20.6f");
  //tree->Scan("time");
  f->Write();
  f->Close();

  return 0;
}


// This software sorts 2 GB MCAL files into increasing time.
void sort_MCAL(TString path_1, TString filename, int minimum, int maximum){
  gROOT->Reset();
  double obt_MCAL, usec;
  float Etot_MCAL;
  int entries_tdata, contact, y, x;

  //#################################################################################
  // Create new tree to store data in
  // the structure to hold the variables for the branch
  struct MCAL {
  };

  MCAL data;

  // create a new ROOT file
  TFile *f = new TFile(path_1 + filename +"_sorted.root","RECREATE");
  TTree *tree = new TTree("data","");
  tree->Branch("time", &obt_MCAL,"time/D");
  tree->Branch("usec", &usec,"usec/D");
  tree->Branch("Etot", &Etot_MCAL,"Etot/F");
  tree->Branch("contact", &contact,"contact/I");
  //#################################################################################

  //################################################################################
  TFile *file_MCAL = new TFile(path_1 + filename + ".root", "r");
  TTree *tdata = (TTree*) file_MCAL->Get("data");
  tdata->SetBranchAddress("time", &obt_MCAL);
  tdata->SetBranchAddress("usec", &usec);
  tdata->SetBranchAddress("Etot", &Etot_MCAL);
  tdata->SetBranchAddress("contact", &contact);

  entries_tdata = tdata->GetEntriesFast();

  tdata->BuildIndex("time", "usec");
  //#################################################################################

  TTreeIndex *index = (TTreeIndex*) tdata->GetTreeIndex();

  int c = 0;

  //int minimum = 385776000 - 60; //23.06.2016 23.59.59 is 385776000
  //int maximum = 393811199 + 60; //23.06.2016 23.59.59 is 393811199

  for (int i = 0; i < index->GetN() - 1; i++){
    Long64_t local = tdata->LoadTree(index->GetIndex()[i]);
    tdata->GetEntry(local);

    c += 1;
    if (c > entries_tdata/100) {
      cout << int((i/(entries_tdata + 0.0))*100) << endl;
      c = 0;
    }

    if (minimum <= obt_MCAL && obt_MCAL <= maximum) {
      tree->Fill();
    }
  }

  int counter = 0; double a = 0; double b = 0; tree->SetBranchAddress("time", &obt_MCAL);
  for (int x = 0; x < tree->GetEntries(); x++) {
    tree->GetEntry(x);
    b = obt_MCAL;

    if (a > b) {
      counter += 1;
    }
    a = obt_MCAL;
  }

  cout << "COUNTER:" << counter << endl;
  delete file_MCAL;
  tree->Show(0);
  tree->Show(tree->GetEntries()-1);
  f->Write();
  f->Close();
}

int save_buildindex(TString path_1, TString filename){

  gROOT->Reset();
  double obt_MCAL, usec;
  float Etot_MCAL;
  int entries_tdata, contact, y, x;



  //################################################################################
  TFile *file_MCAL = new TFile(path_1 + filename + "_sorted.root", "r");
  TTree *tdata = (TTree*) file_MCAL->Get("data");
  tdata->SetBranchAddress("time", &obt_MCAL);
  tdata->SetBranchAddress("usec", &usec);
  tdata->SetBranchAddress("Etot", &Etot_MCAL);
  tdata->SetBranchAddress("contact", &contact);

  entries_tdata = tdata->GetEntriesFast();
  //#################################################################################


  // create a new ROOT file
  TFile *f = new TFile(path_1 + filename +"_sorted_with_index.root","UPDATE");
  TTree *tree = new TTree("data","index");


  tree->Branch("time", &obt_MCAL,"time/D");
  tree->Branch("usec", &usec,"usec/D");
  tree->Branch("Etot", &Etot_MCAL,"Etot/F");
  tree->Branch("contact", &contact,"contact/I");


  int c = 0;
  for ( x = 0; x < entries_tdata; x++) {
    tdata->GetEntry(x);
    tree->Fill();


    c += 1;
    if (c > entries_tdata/100) {
      cout << int((x/(entries_tdata + 0.0))*100) << endl;
      c = 0;
    }
  }
  delete file_MCAL;
  tree->BuildIndex("time", "usec");
  tree->Write();


  f->Write();
  f->Close();


  return 0;
}


int add_usec_MCAL(TString path_1, TString filename){
  double obt_MCAL, usec;
  float Etot_MCAL;
  int entries_tdata, contact{0}, y, x;

  // create a new ROOT file
  TFile *f = new TFile(path_1 + filename +"_usec_sorted.root","RECREATE");
  TTree *tree = new TTree("data","");
  tree->Branch("time", &obt_MCAL,"time/D");
  tree->Branch("usec", &usec,"usec/D");
  tree->Branch("Etot", &Etot_MCAL,"Etot/F");
  tree->Branch("contact", &contact,"contact/I");
  //#################################################################################

  //################################################################################
  TFile *file_MCAL = new TFile(path_1 + filename + "_sorted.root", "r");
  TTree *tdata = (TTree*) file_MCAL->Get("data");
  tdata->SetBranchAddress("time", &obt_MCAL);
  tdata->SetBranchAddress("Etot", &Etot_MCAL);

  entries_tdata = tdata->GetEntries();
  int counter = 0;
  for (size_t x = 0; x < entries_tdata; x++) {
    tdata->GetEntry(x);
    if (counter > entries_tdata/100) {
      cout << int((x/(float(entries_tdata + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;
    usec = round((obt_MCAL - int(obt_MCAL))*1000000);
    tree->Fill();

  }
  tree->Print();
  tree->Show(0);
  //tree->Scan("", "", "col=20.6f");
  f->Write();

  file_MCAL->Close();
  f->Close();




  return 0;
}


int C3_reduce_MCAL_time_Etot(){
  //################ < 2015 #############################

  TString path_less_2015 = "/scratch/Master/fresh_WWLLN_AGILE/AC_enabled_analyse/MCAL_folder/";
  //add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_2981-6050");
  /*
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_2981-6050");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_6051-7999");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_8000-9000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_9001-10000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_10001-10700");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_10701-11300");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_11301-12000_1");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_11301-12000_2");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_12001-12500");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_12501-13300_1");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_12501-13300_2");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_12501-13300_3");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_13301-14000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_14001-15000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_15001-16000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_16001-17000");

  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_17001-18000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_18001-19000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_19001-20000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_20001-20700");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_20701-21430");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_21431-22200");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_22201-22700");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_22701-23000");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_23001-23400");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_23401-23800");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_23801-24200");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_24201-25400");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_25401-30800");
  add_usec_MCAL(path_less_2015, "MCAL_time_Etot_2_30801-40900");
  */
  //################ 2015 #############################
  TString path_AC_disabled = "/scratch/Master/fresh_WWLLN_AGILE/AC_disabled_analyse/";
  //reduce_MCAL_time_Etot(path_AC_disabled, "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2015/", 40900, 42298);
  //sort_MCAL(path_AC_disabled, "MCAL_time_Etot_40900_42298", 354153660 - 60, 362188799 + 60 );



  //################ 2016 #############################
  //reduce_MCAL_time_Etot("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/", "/Volumes/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", 45604, 47597);
  //sort_MCAL("/Volumes/TOSHIBA/fresh_WWLLN_AGILE/Time_drift/", "MCAL_time_Etot_45604_47597", 385776000 - 60, 393811199 + 60);


  //################ Timedrift 2017 #############################
  //TString path_timedrift_2017 = "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017/";
  //reduce_MCAL_time_Etot(path_timedrift_2017, "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", 45000, 49700);
  //reduce_MCAL_time_Etot(path_timedrift_2017, "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", 49701, 51000);
  //reduce_MCAL_time_Etot(path_timedrift_2017, "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", 51001, 52600);
  //reduce_MCAL_time_Etot(path_timedrift_2017, "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", 52601, 53800);
  //reduce_MCAL_time_Etot(path_timedrift_2017, "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/MCAL_2016/", 53801, 54583);


  //sort_MCAL(path_timedrift_2017, "MCAL_time_Etot_45000_49700", 378691200 - 60, 441849599 + 60); //01.01.2016 clock 00.00.00 and 31.12.2017 clock 23.59.59
  //sort_MCAL(path_timedrift_2017, "MCAL_time_Etot_49701_51000", 378691200 - 60, 441849599 + 60); //01.01.2016 clock 00.00.00 and 31.12.2017 clock 23.59.59
  //sort_MCAL(path_timedrift_2017, "MCAL_time_Etot_51001_52600", 378691200 - 60, 441849599 + 60); //01.01.2016 clock 00.00.00 and 31.12.2017 clock 23.59.59
  //sort_MCAL(path_timedrift_2017, "MCAL_time_Etot_52601_53800", 378691200 - 60, 441849599 + 60); //01.01.2016 clock 00.00.00 and 31.12.2017 clock 23.59.59 //THINKPAD
  //sort_MCAL(path_timedrift_2017, "MCAL_time_Etot_53800_54583", 378691200 - 60, 441849599 + 60); //01.01.2016 clock 00.00.00 and 31.12.2017 clock 23.59.59

  //################ Timedrift 2017 extra  #############################

  TString path_timedrift_2017_extra = "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/Time_drift_2017_extra/";
  //reduce_MCAL_time_Etot(path_timedrift_2017_extra, "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/MCAL_55655_55935/", 55655, 55935);
  //
  sort_MCAL(path_timedrift_2017_extra, "MCAL_time_Etot_55655_55935", 443332000 - 500, 445029000 + 500);

  return 0;

}


int copy_WWLLN_test(){

  int year_W, month_W, day_W, hour_W, minute_W, sec_W, usec_W ;
  float lat_W_deg, lon_W_deg;
  double time_WWLLN;

  year_W = 2015; month_W = 6; day_W = 23; hour_W = 23; minute_W = 59; sec_W = 59; usec_W = 0;
  TTimeStamp *epoch = new TTimeStamp(2004, 1, 1, 0, 0, 0, 0, 1, 0);
  // Convert WWLLN time to second since epoch. This is to compare with AGILE onboard time
  TTimeStamp *t = new TTimeStamp(year_W, month_W, day_W, hour_W, minute_W, sec_W, 0,1,0);
  long time_WWLLN_long =  (t->GetSec()  - epoch->GetSec());
  time_WWLLN = time_WWLLN_long + usec_W*pow(10,-6);
  cout.precision(15);
  cout << time_WWLLN << endl;


  return 0;
}
