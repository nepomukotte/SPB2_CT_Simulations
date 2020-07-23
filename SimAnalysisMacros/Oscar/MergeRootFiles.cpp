#include <TChain.h>
#include <TVectorD.h>

void readDir(string dirname,vector<string> *names) {

    string file;
  void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname.c_str()));
    const char* strptr;
    if (dir) {
        while ((strptr = gSystem->GetDirEntry(dir))) {
          file = strptr;
          //if (file=="." || file==".." || file=="log") continue;
          //names->push_back(file);
          size_t found;
          found=file.find(".root");
          if (found==string::npos) continue;
          names->push_back(file);
        }
        gSystem->FreeDirectory(dir);
    }
     // for(size_t i(0); i<names->size(); i++)
     //   Printf("%s",(*names)[i].c_str());
}

int MergeRootFiles (){
	string DataDir = "/home/oscar/Documents/Research/NSB_Traces";
	vector<string> * filenames = new vector<string>; 
	readDir(DataDir,filenames);
	

	TFile *fl_merge = new TFile("NSB_Traces_Merged.root","RECREATE");
	TChain *ch_Glob = 0;
	TChain *ch_T0_Glob = 0;
	TDirectory *dir_events = fl_merge->mkdir("Events","Events");
	fl_merge->cd("Events");

	TList *lst_Files = new TList();

	for (int i = 0; i<filenames->size(); i++){
		stringstream fnamestream;
		fnamestream<<DataDir<<"/"<<(*filenames)[i].c_str();
		string fstring(fnamestream.str());
		cout<<fstring<<endl;
		lst_Files->Add(TFile::Open(fstring.c_str()));
	}

	TString path( (char*)strstr( dir_events->GetPath(), ":" ) );
	path.Remove( 0, 2 );
	

	cout<<"Working Path: "<<path<<endl;

	TFile *first_source = (TFile*)lst_Files->First();
	first_source->cd(path);
	gDirectory->pwd();
	TDirectory *current_sourcedir = gDirectory;
	//TString patha( (char*)strstr( current_sourcedir->GetPath(), ":" ) );
	//patha.Remove( 0, 2 );
	//cout<<"Working Path: "<<patha<<endl;

	ch_Glob = new TChain("Events/tSimulatedEvents");
	ch_Glob->Add(first_source->GetName());

	TFile *nextsource = (TFile*)lst_Files->After(first_source);
	
	while(nextsource){
		nextsource->cd(path);
		ch_Glob->Add(nextsource->GetName());
		nextsource = (TFile*)lst_Files->After(nextsource);
		cout<<"Hello"<<endl;
	}

	ch_T0_Glob = new TChain("Events/T0");
	ch_T0_Glob->Add(first_source->GetName());

	nextsource = (TFile*)lst_Files->After(first_source);
	
	while(nextsource){
		nextsource->cd(path);
		ch_Glob->Add(nextsource->GetName());
		nextsource = (TFile*)lst_Files->After(nextsource);
		cout<<"Hello"<<endl;
	}

	dir_events->cd();
	dir_events->GetFile();
	ch_Glob->Merge(dir_events->GetFile(),0,"keep");
	ch_T0_Glob->Merge(dir_events->GetFile(),0,"keep");
	dir_events->SaveSelf(true);
	return 0;


}