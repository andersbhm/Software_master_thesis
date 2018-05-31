// Written by Anders Lindanger. Master thesis in spacephysics 2018

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <fstream>

using namespace std;

// This software convert WWLLN files from .loc to .root for faster computations.

// input: WWLLN.loc files
// Output: WWLLN.root files

void create_2(TString inputfile){
    gROOT->Reset();
    // the structure to hold the variables for the branch
    struct WWLLN {
        int year;
        int month;
        int day;
        int hour;
        int minute;
        int second;
        int usec;
        float latitude;
        float longitude;
    };

    WWLLN data;
    // open the ASCII file

    FILE *fp = fopen(inputfile + ".loc","r");

    char line[81];

    // create a new ROOT file
    TFile *f = new TFile(inputfile + ".root","RECREATE");

    // create tree
    TTree *tree = new TTree("WWLLN data",inputfile );

    // create one branch with all information from the stucture
    tree->Branch("year", &data.year,"year/I");
    tree->Branch("month", &data.month,"month/I");
    tree->Branch("day", &data.day,"day/I");
    tree->Branch("hour", &data.hour,"hour/I");
    tree->Branch("minute", &data.minute,"minute/I");
    tree->Branch("second", &data.second,"second/I");
    tree->Branch("usec", &data.usec,"usec/I");
    tree->Branch("latitude", &data.latitude,"latitude/F");
    tree->Branch("longitude", &data.longitude,"longitude/F");

    // fill the tree from the values in ASCII file
    while (fgets(line,80,fp)) {
        sscanf(&line[0]," %d / %d / %d , %d : %d : %d . %d , %f , %f",&data.year, &data.month, &data.day, &data.hour, &data.minute, &data.second, &data.usec, &data.latitude,&data.longitude);
        tree->Fill();
    }

    // check what the tree looks like
    //tree->Print();
    //tree->Show(0);
    // tree->Scan("year:month:day:hour:minute:second:usec:latitude:longitude");

    fclose(fp);
    f->Write();
}


// For loop for running the actual create-tree program
void create_1(TString path, TString name, int days_in_month, int start_day){
    TString days_in_month_str = days_in_month;

    for (int i = start_day; i < days_in_month + 1 ; i++){
        TString day_str = Form("%02d", i);
        TString inputfile = path + name + day_str;
        create_2(inputfile);
    }
}

// main function which decide which year and month to open, and start and end of month.
int A3_create_tree_WWLLN(){

    TString path = "/media/fer003/TOSHIBA/WWLLN_AGILE/WWLLN_data/2018/";


        //create_1(path, "A201701",  31,  1);
        //create_1(path, "A201702",  28,  1);
        //create_1(path, "A201703",  31,  1);
        //create_1(path, "A201704",  30,  1);
        //create_1(path, "A201705",  31,  1);
        //create_1(path, "A201706",  30,  1);
        //create_1(path, "A201707",  31,  1);
        //create_1(path, "A201708",  31,  1);
        //create_1(path, "A201709",  30,  1);
        //create_1(path, "A201710",  31,  1);
        //create_1(path, "A201711",  30,  5);
        //create_1(path, "A201712",  1,  1);

        //
        create_1(path, "A201801",  31,  1);
        //
        create_1(path, "A201802",  3,  1);

    return 0;

}
