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
	string MainDir = "/storage/hive/project/phy-otte/shared/Merged_SPB2_CARE_NSB_Sims/NSB_Traces/F1Filter/";
	string DataDir = MainDir+"Event_Traces";
	vector<string> * filenames = new vector<string>; 
	readDir(DataDir,filenames);
	
	string outFile =  MainDir+"Merged_NSB_Traces.root";	

	TFile *fl_merge = new TFile(outFile.c_str(),"RECREATE");
	TChain *ch_Glob = 0;
	TChain *ch_T0_Glob = 0;
	TDirectory *dir_events = fl_merge->mkdir("Events","Events");
	//fl_merge->cd("Events");
	

	int iNumPixels = 512;

	TTree *t_merg_tSim = new TTree("tSimulatedEvents", "tSimulatedEvents");
	TTree *t_merg_T0 = new TTree("T0","T0");

	UInt_t uNumTelescopes;	
	Float_t energy;
	UInt_t eventNumber ; 
	Float_t xcore ;
	Float_t ycore ;
	Float_t azPrim ;
	Float_t znPrim ;
	Bool_t arrayTriggerBit ;
	std::vector< Bool_t > *vTelescopeTriggerBits =0 ;
	Float_t DeltaTL3 ;

	UInt_t uNumTelescopes_t;
	Float_t energy_t;
	UInt_t eventNumber_t ; 
	Float_t xcore_t;
	Float_t ycore_t;
	Float_t azPrim_t;
	Float_t znPrim_t;
	Bool_t arrayTriggerBit_t;
	std::vector< Bool_t > *vTelescopeTriggerBits_t =0 ;
	Float_t DeltaTL3_t;
	vector< vector<Int_t> *>   iFADCTraceInPixel;
	vector<Int_t>   *iPEInPixel;

	vector<int> *vTriggerCluster=0;
	vector<Float_t> *fTimeOverThreshold=0;
	vector<Float_t>   *fSumTimeInPixel=0; 
	vector<Int_t>   *iQDCInPixel=0;
	Int_t           iNumPhotonsInFocalPlane=0;
	Float_t         fAzTel=0 ;
	Float_t         fZnTel =0;
	vector<Bool_t>  *bInLoGain=0;
	iFADCTraceInPixel.assign(iNumPixels,0);
	cout<<&iFADCTraceInPixel[0]<<endl;

	vector<int> *vTriggerCluster_t=0;
	vector<Float_t> *fTimeOverThreshold_t=0;
	vector<Float_t>   *fSumTimeInPixel_t=0; 
	vector<Int_t>   *iQDCInPixel_t=0;
	Int_t           iNumPhotonsInFocalPlane_t=0;
	Float_t         fAzTel_t=0 ;
	Float_t         fZnTel_t=0;
	vector<Bool_t>  *bInLoGain_t=0;
	vector< vector<Int_t> *>   iFADCTraceInPixel_t;
	vector<Int_t>   *iPEInPixel_t=0;
	iFADCTraceInPixel_t.assign(iNumPixels,0);
	cout<<&iFADCTraceInPixel_t[0]<<endl;

	t_merg_tSim->Branch("energy",&energy);
	t_merg_tSim->Branch("ZnPrim",&znPrim);
	t_merg_tSim->Branch("AzPrim",&azPrim);
	t_merg_tSim->Branch("xcore",&xcore);
	t_merg_tSim->Branch("ycore",&ycore);
	t_merg_tSim->Branch("arrayTriggerBit",&arrayTriggerBit);
	t_merg_tSim->Branch("uNumTelescopes",&uNumTelescopes);
	t_merg_tSim->Branch("eventNumber",&eventNumber);
	t_merg_tSim->Branch("DeltaTL3",&DeltaTL3);
	t_merg_tSim->Branch("vTelescopeTriggerBits",&vTelescopeTriggerBits);

	t_merg_T0->Branch("vGroupsInTriggerCluster",&vTriggerCluster);
	t_merg_T0->Branch("vTimeOverThreshold",&fTimeOverThreshold);
	t_merg_T0->Branch("vSumTimeInPixel", &fSumTimeInPixel);
	t_merg_T0->Branch("vPEInPixel", &iPEInPixel);
	t_merg_T0->Branch("vQDCValue", &iQDCInPixel);
	t_merg_T0->Branch("iPhotonsInFocalPlane", &iNumPhotonsInFocalPlane);
	t_merg_T0->Branch("fAzTel", &fAzTel);
	t_merg_T0->Branch("fZnTel", &fZnTel);
	t_merg_T0->Branch("vHiLoGainBit", &bInLoGain);
	TString name;
	for(int g=0;g<iNumPixels;g++){
		name.Form("vFADCTraces%i",g);
		//cout<<name<<endl;
		t_merg_T0->Branch(name,&iFADCTraceInPixel[g]);
	}

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

	
	TFile *source = (TFile*)lst_Files->First();
	TTree *t_temp_tSim;
	TTree *t_temp_T0;
	while(source){
		t_temp_tSim = (TTree*)source->Get("Events/tSimulatedEvents");
		t_temp_T0 = (TTree*)source->Get("Events/T0");

		t_temp_tSim->SetBranchAddress("energy",&energy_t);
		t_temp_tSim->SetBranchAddress("ZnPrim",&znPrim_t);
		t_temp_tSim->SetBranchAddress("AzPrim",&azPrim_t);
		t_temp_tSim->SetBranchAddress("xcore",&xcore_t);
		t_temp_tSim->SetBranchAddress("ycore",&ycore_t);
		t_temp_tSim->SetBranchAddress("arrayTriggerBit",&arrayTriggerBit_t);
		t_temp_tSim->SetBranchAddress("uNumTelescopes",&uNumTelescopes_t);
		t_temp_tSim->SetBranchAddress("eventNumber",&eventNumber_t);
		t_temp_tSim->SetBranchAddress("DeltaTL3",&DeltaTL3_t);
		t_temp_tSim->SetBranchAddress("vTelescopeTriggerBits",&vTelescopeTriggerBits_t);
	
		t_temp_T0->SetBranchAddress("vGroupsInTriggerCluster",&vTriggerCluster_t);
		t_temp_T0->SetBranchAddress("vTimeOverThreshold",&fTimeOverThreshold_t);
		t_temp_T0->SetBranchAddress("vSumTimeInPixel", &fSumTimeInPixel_t);
		t_temp_T0->SetBranchAddress("vPEInPixel", &iPEInPixel_t);
		t_temp_T0->SetBranchAddress("vQDCValue", &iQDCInPixel_t);
		t_temp_T0->SetBranchAddress("iPhotonsInFocalPlane", &iNumPhotonsInFocalPlane_t);
		t_temp_T0->SetBranchAddress("fAzTel", &fAzTel_t);
		t_temp_T0->SetBranchAddress("fZnTel", &fZnTel_t);
		t_temp_T0->SetBranchAddress("vHiLoGainBit", &bInLoGain_t);
		TString name_t;
		for(int g=0;g<iNumPixels;g++){
			name_t.Form("vFADCTraces%i",g);
		    //cout<<name<<endl;
		    t_temp_T0->SetBranchAddress(name_t,&iFADCTraceInPixel_t[g]);
		}

		for(int k=0; k<t_temp_tSim->GetEntries(); k++){
			t_temp_tSim->GetEntry(k);
			t_temp_T0->GetEntry(k);

			energy = energy_t;
			znPrim = znPrim_t;
			azPrim = azPrim_t;
			xcore = xcore_t;
			ycore = ycore_t;
			arrayTriggerBit = arrayTriggerBit_t;
			uNumTelescopes = uNumTelescopes_t;
			eventNumber = eventNumber_t;
			DeltaTL3 = DeltaTL3_t;
			vTelescopeTriggerBits = vTelescopeTriggerBits_t;

			vTriggerCluster = vTriggerCluster_t;
			fTimeOverThreshold = fTimeOverThreshold_t;
			fSumTimeInPixel = fSumTimeInPixel_t;
			iPEInPixel = iPEInPixel_t;
			iQDCInPixel = iQDCInPixel_t;
			iNumPhotonsInFocalPlane = iNumPhotonsInFocalPlane_t;
			fAzTel = fAzTel_t;
			fZnTel = fZnTel_t;
			bInLoGain = bInLoGain_t;
			iFADCTraceInPixel = iFADCTraceInPixel_t;

			t_merg_tSim->Fill();
			t_merg_T0->Fill();
		}

		t_temp_T0->Delete();
		t_temp_tSim->Delete();
		source = (TFile*)lst_Files->After(source);
	}

	fl_merge->cd("Events");
	t_merg_T0->Write();
	t_merg_tSim->Write();

	//fl_merge->Write();
	fl_merge->Close();

		
	return 0;


}
