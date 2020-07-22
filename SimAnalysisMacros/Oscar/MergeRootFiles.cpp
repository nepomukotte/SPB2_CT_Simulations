#include <TChain>
#include <TVectorD>

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
	string DataDir = "";
	vector<string> * filenames = new vector<string>; 
	readDir(DataDir,filenames);

	TChain *ch_Tree = new TChain("Events");

	for(int i = 0; i<filenames->size(); i++){

		stringstream fnamestream;
		fnamestream<<DataDir<<"/"<<(*filenames)[i].c_str();
		string fstring(fnamestream.str());
		cout<<"opening file: "<<fstring<<endl;
		f = new TFile(fstring.c_str(),"READ");
		if( f->IsZombie() ){
		   cout << "error opening root input file: " <<fnamestream.str().c_str()  << endl;
		   cout << "...skipping" << endl;
		    //exit( -1 );
		}
		else{
			ch_Tree->Add(fstring.c_str());
		}
	}

	ch_Tree->Merge("NSB_Traces_Merged.root")

	return 0;


}