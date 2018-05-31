// Written by Anders Lindanger. Master thesis in spacephysics 2018

void add_usec_GRID(TString path_1, TString filename){
  gROOT->Reset();
  double obt_MCAL;
  int usec;
  float Etot_MCAL;
  int entries_tdata, contact, y, x;

  // create a new ROOT file
  TFile *f = new TFile(path_1 + filename +"_usec.root","RECREATE");
  TTree *tree = new TTree("data_GRID","");
  tree->Branch("obt", &obt_MCAL,"time/D");
  tree->Branch("usec", &usec,"usec/I");
  tree->Branch("contact", &contact,"contact/I");

  //################################################################################
  TFile *file_MCAL = new TFile(path_1 + filename + ".root", "r");
  TTree *tdata = (TTree*) file_MCAL->Get("data_GRID");

  //TChain *tdata = new TChain("data_GRID");
  //tdata->Add(path_1 + "GRID_2008_2015/*.root");
  tdata->SetBranchAddress("obt", &obt_MCAL);
  tdata->SetBranchAddress("contact", &contact);
  entries_tdata = tdata->GetEntries();
  cout << entries_tdata << endl;

  set <double> obt_list;
  int counter = 0;
  for (int x = 0; x < entries_tdata; x++) {
    tdata->GetEntry(x);


    if (counter > entries_tdata/100) {
      cout << int((x/(float(entries_tdata + 0.0)))*100) << endl;
      counter = 0;
    }
    counter += 1;

    usec = int(round((obt_MCAL - int(obt_MCAL))*1000000));

    //if (contact != 21680) {
      //obt_list.insert(obt_MCAL);
      tree->Fill();
    //}


  }
  tree->Print();
  tree->Show(0);
  tree->Scan("obt:usec:contact", "", "col=20.6f");
  cout << "Entries TChain: "<< entries_tdata << endl;
  f->Write();
}

void sort_GRID(TString path_1, TString filename){
  gROOT->Reset();
  double obt_MCAL;
  float Etot_MCAL;
  int entries_tdata, contact, usec,y, x;

  //#################################################################################
  // Create new tree to store data in
  // the structure to hold the variables for the branch
  struct MCAL {
  };

  MCAL data;

  // create a new ROOT file
  TFile *f = new TFile(path_1 + filename +"_sorted.root","RECREATE");
  TTree *tree = new TTree("data_GRID","");
  tree->Branch("obt", &obt_MCAL,"time/D");
  tree->Branch("usec", &usec,"usec/I");
  tree->Branch("contact", &contact,"contact/I");
  //#################################################################################

  //################################################################################
  TFile *file_MCAL = new TFile(path_1 + filename + ".root", "r");
  TTree *tdata = (TTree*) file_MCAL->Get("data_GRID");
  tdata->SetBranchAddress("obt", &obt_MCAL);
  tdata->SetBranchAddress("usec", &usec);
  tdata->SetBranchAddress("contact", &contact);

  entries_tdata = tdata->GetEntriesFast();

  tdata->BuildIndex("obt", "usec");
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

    //if (minimum <= obt_MCAL && obt_MCAL <= maximum) {

      tree->Fill();
    //}
  }

  int counter = 0; double a = 0; double b = 0; tree->SetBranchAddress("obt", &obt_MCAL);
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
  tree->Scan("obt:usec:contact", "", "col=20.6f");
  f->Write();
  f->Close();
}


int C_sort_GRID(){

  TString path = "/media/fer003/TOSHIBA/fresh_WWLLN_AGILE/GRID/2008_2015/";
  //add_usec_GRID(path, "GRID_contact_21680_40575");

  //sort_GRID(path,"GRID_contact_5720_21679_usec" );
  //
  sort_GRID(path,"GRID_contact_21680_40575_usec" );
  return 0;
}
