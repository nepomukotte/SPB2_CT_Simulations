#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TH2.h>
#include <TTimer.h>


void FindBin(int iPix,int *nx, int *ny)
{
	//MUSIC ID
	int MUSIC_ID = iPix/8;	
	//Channel in MUSIC
	int MUSIC_Channel = iPix%8; 	
	//row
	int MUSIC_row = MUSIC_ID/16;
	
	//column    
	int MUSIC_column = MUSIC_ID%16;	
	*ny = 4*MUSIC_row+MUSIC_Channel%4+1; 
	*nx = 2*MUSIC_column+MUSIC_Channel/4+1; 
	// cout<<"Pixel "<<iPix<<" x: "<<*nx<<" y: "<<*ny<<endl;

}

int SignalExtraction (){
	int nBins;
	vector <int> *iPePix=0;
	vector<vector<int>*> iTrace;

	int nx;
	int ny;

	int nPixels = 512;

	bool *isTriggered=0;
	int cEvent;

	double sum_base=0;
	double avg_base;
	double integ;
	int nSampleOffset=3;

	iTrace.assign(nPixels,0);



	TFile *file = new TFile( "/home/oscar/Documents/Research/Testing_Events1.root", "READ" );
	TH2 *h_camera = new TH2F ("h_camera", "Camera", 32,-0.5,32.5,16,-0.5,15.5);
	TH2 *h_camera_integ = new TH2F ("h_camera_integ", "Camera Charge", 32,-0.5,32.5,16,-0.5,15.5);
	TTree *tSimulatedEvents = (TTree*)file->Get( "Events/tSimulatedEvents" );
	TTree *tEvent = (TTree*)file->Get( "Events/T0" );

	for (int i = 0; i<nPixels; i++){
		tEvent->SetBranchAddress(TString::Format("vFADCTraces%i",i),&iTrace[i]);
	}
	tSimulatedEvents->SetBranchAddress("arrayTriggerBit",&isTriggered);
	tEvent->SetBranchAddress("vPEInPixel", &iPePix);

	h_camera->SetStats(0);
	

	cEvent = 0;
	tSimulatedEvents->GetEntry(0);
	tEvent->GetEntry(0);
	/*while (!isTriggered){
		cEvent++;
		tSimulatedEvents->GetEntry(cEvent);
		tEvent->GetEntry(cEvent);
	}*/
	cout << cEvent <<endl;

	tSimulatedEvents->GetEntry(9472);
	tEvent->GetEntry(9574);
	for (int i=0; i<nPixels; i++){
		sum_base=0;
		FindBin(i,&nx,&ny);
		for (int j = 0; j<nSampleOffset; j++){
			sum_base += iTrace[i]->at(j);
		}
		avg_base = sum_base*(0.333);
		//cout << avg_base << endl;
		h_camera->SetBinContent(nx,ny,iPePix->at(i));

		integ = 0;
		for (int j = 0; j<iTrace[i]->size(); j++){
			integ += iTrace[i]->at(j)-avg_base;
		}
		h_camera_integ->SetBinContent(nx,ny,integ);
	}

	TCanvas *c_camera = new TCanvas("c_camera", "Camera Canvas",800,500);
	h_camera->Draw("colz");

	TCanvas *c_camera_charge = new TCanvas("c_camera_charge", "Camera Charge Canvas",800,500);
	h_camera_integ->Draw("colz");
	h_camera_integ->SetStats(0);

	return 0;
}
