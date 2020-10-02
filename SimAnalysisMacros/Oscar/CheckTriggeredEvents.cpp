#include <TTree.h>

int CheckTriggeredEvents(){
	TFile *file = new TFile( "NSB_Traces_Merged_Fixed.root", "READ" );
	TTree *tSimulatedEvents = (TTree*)file->Get( "Events/tSimulatedEvents" );
	bool *isTriggered=0;
	tSimulatedEvents->SetBranchAddress("arrayTriggerBit",&isTriggered);
	int nTriggered=0;
	int nEvents;

	for (int i = 0; i<tSimulatedEvents->GetEntries(); i++){
		tSimulatedEvents->GetEntry(i);
		if(isTriggered){
			nTriggered++;
		}
	}
	cout<<nTriggered<<endl;
	return 0;
}