#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TH2.h>
#include <TTimer.h>
#include <cmath>


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

void PopulateADCPe(){

}

int FindFirstPix(int iMUSIC){
	return iMUSIC*8;
}

int SignalExtraction (){
	int nBins;
	vector <int> *iPePix=0;
	vector<vector<int>*> iTrace;
	vector<int> *iTrigMusic = 0;

	int nx;
	int ny;

	int nPixels = 512;

	bool *isTriggered=0;
	int cEvent;
	int nEvents=10000;
	int iPixSignal=0;
	int iPixSignalID;

	double dynamicRangeAGET = 750;//mV
	double mvPerPe = 1.92;
	double nBitsAGET = 12;


	double sum_base=0;
	double avg_base;
	double integ;
	int nSampleOffset=3;
	double px_rms;
	iTrace.assign(nPixels,0);



	TFile *file = new TFile( "/home/oscar/Documents/Research/TestAnalog.root", "READ" );
	TH2 *h_camera = new TH2F ("h_camera", "Camera", 32,-0.5,32.5,16,-0.5,15.5);
	TH2 *h_camera_integ = new TH2F ("h_camera_integ", "Camera Charge", 32,-0.5,32.5,16,-0.5,15.5);
	TH1 *h_base = new TH1F("h_base", "Pedestal", 4096,0,4096);
	TH1 *h_rms = new TH1F("h_rms", "RMS Of Pedestal", nPixels,0,nPixels);
	TH1 *h_rms_dist = new TH1F("h_rms_dist","Camera RMS", 1000,0,10);
	TH2 *h_pe_adc = new TH2F("h_pe_adc","ADC per Pe",500,0,500,4096,0,4096);
	TH1 *h_trace;
	TTree *tSimulatedEvents = (TTree*)file->Get( "Events/tSimulatedEvents" );
	TTree *tEvent = (TTree*)file->Get( "Events/T0" );
	TGraph *gr_mean = new TGraph();

	for (int i = 0; i<nPixels; i++){
		tEvent->SetBranchAddress(TString::Format("vFADCTraces%i",i),&iTrace[i]);
	}
	tSimulatedEvents->SetBranchAddress("arrayTriggerBit",&isTriggered);
	tEvent->SetBranchAddress("vPEInPixel", &iPePix);
	tEvent->SetBranchAddress("vGroupsInTriggerCluster", &iTrigMusic);

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

	//tSimulatedEvents->GetEntry(9574);
	//tEvent->GetEntry(9574);
	for (int k = 0; k<nEvents; k++){
		tSimulatedEvents->GetEntry(k);
		tEvent->GetEntry(k);
		if(isTriggered){
			h_base->Reset();
			h_camera_integ->Reset();

			//cout<<iTrigMusic->size()<<endl;
			for (int l = 0; l<iTrigMusic->size(); l++){
				for(int i=FindFirstPix(iTrigMusic->at(l)); i<FindFirstPix(iTrigMusic->at(l))+8; i++){
					sum_base=0;
					FindBin(i,&nx,&ny);
					cout<<"Hello"<<endl;
					for (int j = 0; j<nSampleOffset-1; j++){
						sum_base += iTrace[i]->at(j);
						//px_rms += iTrace[i]->at(j)*iTrace[i]->at(j);
						gr_mean->SetPoint(nSampleOffset*i+j,i,iTrace[i]->at(j));
						h_base->Fill(iTrace[i]->at(j));
					}
					avg_base = sum_base/nSampleOffset;

					h_camera->SetBinContent(nx,ny,iPePix->at(i));
					
					integ = 0;
					for (int j = 0; j<iTrace[i]->size(); j++){
						if(iTrace[i]->at(j)-avg_base>0){
							integ += iTrace[i]->at(j)-avg_base;
						}
						
					}

					h_camera_integ->SetBinContent(nx,ny,integ);
					h_pe_adc->Fill(iPePix->at(i),integ);

					if(iPixSignal<integ){
						iPixSignal=integ;
						iPixSignalID=i;
					}
				}

				cout<<"Signal is in Pixel: "<<iPixSignalID<<endl;
				iPixSignal=0;

			}
		}
			/*for (int i=0; i<nPixels; i++){
				sum_base=0;
				//px_rms=0;
				
				//px_rms = sqrt(px_rms/nSampleOffset);
				
				//cout << avg_base << endl;
				
		
				integ = 0;
				if(k==9999 && i==184){
					h_trace = new TH1F("h_trace","Trace", iTrace[i]->size(),0,iTrace[i]->size());
					for (int j = 0; j<iTrace[i]->size();j++){
						h_trace->Fill(j,iTrace[i]->at(j)-avg_base);
					}
				}
				for (int j = 0; j<iTrace[i]->size(); j++){
					integ += iTrace[i]->at(j)-avg_base;
					
				}
				h_camera_integ->SetBinContent(nx,ny,integ);
				h_pe_adc->Fill(iPePix->at(i),integ);
			}*/
			//cout<<k<<endl;
			//h_rms->Fill(,h_base->GetRMS());
			//h_rms->SetBinError(h_rms->FindBin(k),0);
		h_rms_dist->Fill(h_base->GetRMS());
			//cout<<h_base->GetRMS()<<endl;
			//cout<<k<<endl 257.5;
	}
	
	h_pe_adc->Fit("pol 1","Q","E");

	TCanvas *c_camera = new TCanvas("c_camera", "Camera Canvas",800,500);
	h_camera->Draw("colz");

	TCanvas *c_camera_charge = new TCanvas("c_camera_charge", "Camera Charge Canvas",800,500);
	h_camera_integ->Draw("colz");
	h_camera_integ->SetStats(0);

	TCanvas *c_mean = new TCanvas("c_mean", "Mean Dispersions", 800,500);
	gr_mean->GetHistogram()->Draw();
	gr_mean->GetHistogram()->GetXaxis()->SetTitle("Pixel Number");
	gr_mean->GetHistogram()->GetYaxis()->SetTitle("Pedestal");
	gr_mean->Draw("P");

	/*TCanvas *c_rms = new TCanvas("c_rms", "Canvas RMS", 800,500);
	h_rms->Draw("HIST");
	h_rms->SetStats(0);
	h_rms->GetXaxis()->SetTitle("Pixel");
	h_rms->GetYaxis()->SetTitle("RMS");
	h_rms->SetStats(0);
	h_rms->SetMarkerStyle(kFullCircle);*/

	TCanvas *c_base = new TCanvas("c_base", "Canvas Baseline", 800,500);
	h_base->Draw("HIST");
	h_base->SetStats(0);
	h_base->GetXaxis()->SetTitle("DC");
	h_base->GetYaxis()->SetTitle("Frequency");
	h_base->GetXaxis()->SetRangeUser(490,510);

	
	TCanvas *c_rms = new TCanvas("c_rms", "Canvas RMS", 800,500);
	h_rms_dist->Draw("HIST");
	h_rms_dist->SetStats(0);
	h_rms_dist->SetLineColor(kBlack);
	h_rms_dist->GetXaxis()->SetTitle("RMS");
	h_rms_dist->GetYaxis()->SetTitle("Frequency");
	h_rms_dist->GetXaxis()->SetRangeUser(0,10);	

	TCanvas *c_pe_chg = new TCanvas("c_pe_chg","Canvas Pe Charge",800,500);
	h_pe_adc->Draw("SCAT");

	TCanvas *c_trace = new TCanvas("c_trace","Trace",800,500);
	h_trace->Draw("HIST");
	h_trace->GetXaxis()->SetTitle("Sample [10 ns]");
	h_trace->GetYaxis()->SetTitle("DC");

	return 0;
}
