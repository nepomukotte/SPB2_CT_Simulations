#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TH2.h>
#include <TTimer.h>
#include <cmath>
#include <iostream>
#include <fstream>


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

bool isNSBTriggerOld(int *signalID){
	return TMath::Abs(signalID[0]-signalID[1])!=8;
}

int FindFirstPix(int iMUSIC){
	return iMUSIC*8;
}

bool isNSBTrigger(vector<int> *iTrigMusic, vector<vector<int>*> iTrace){
	int max_val=0;
	int max_val_Sample_Idx=0;
	int max_val_Pix_Idx;
	int max_temp = 0;
	int bifoc_Pix_Idx;
	int bifoc_max_val;
	int bifoc_Sample_Idx;
	bool isNoise;
	int tolerance = 5;
	//int tolerance = 2;

	if (iTrigMusic->size()!=0){
		for(int i=FindFirstPix(iTrigMusic->at(0)); i<FindFirstPix(iTrigMusic->at(0))+8; i++){
			max_temp = (int)*std::max_element(iTrace[i]->begin(),iTrace[i]->end());
			if(max_temp>max_val){
				max_val=max_temp;
				max_val_Sample_Idx = std::distance(iTrace[i]->begin(),std::max_element(iTrace[i]->begin(),iTrace[i]->end()));
				max_val_Pix_Idx = i;
			}
		}
	
		bifoc_Pix_Idx = max_val_Pix_Idx + 8;
		bifoc_Sample_Idx = std::distance(iTrace[bifoc_Pix_Idx]->begin(),std::max_element(iTrace[bifoc_Pix_Idx]->begin(),iTrace[bifoc_Pix_Idx]->end()));
		bifoc_max_val = *std::max_element(iTrace[bifoc_Pix_Idx]->begin(),iTrace[bifoc_Pix_Idx]->end());
		if(TMath::Abs(bifoc_Sample_Idx - max_val_Sample_Idx)<=1){
			isNoise = false;
			
		}
		else if(TMath::Abs(bifoc_max_val - max_val)<=tolerance){
			isNoise = false;
			
		}
		else{
			isNoise = true;
		}
		cout<<bifoc_Sample_Idx<<" "<<max_val_Sample_Idx<<endl;
		cout<<bifoc_max_val<<" "<<max_val<<endl;
		//cout<<iTrigMusic->at(0)<<" "<<iTrigMusic->at(1)<<endl;
	}
	else{
		isNoise = true;
	}

	
	return isNoise;
}

bool isAboveThreshold (double *max_val, int thresh){
	bool isLarger;
	
	if(max_val[0]>=thresh && max_val[1]>=thresh){
		isLarger=true;
	}
	else{
		isLarger=false;
	}

	return isLarger;
}


int FindPixelNumber(int nx, int ny){
	return nx*4+128*(ny/4)+ny%4;
}

int FindTraceWidth (int *maxSampleID, int *iPixSignalID, vector<vector<int>*>iTrace, int offset,double baseline){
	int full_width = 0;
	int i = maxSampleID[0];
	bool isFoundEnd = false;
	do{
		if (iTrace[iPixSignalID[0]]->at(i)<baseline){
			full_width = i;
			isFoundEnd =true;
		}
		else{
			i++;
		}
	}while(!isFoundEnd && i<iTrace[iPixSignalID[0]]->size());

	return full_width;
}

void ExtractSignalTraces (TH2 *h_charge, TH1 *h_trace_iso[15], int *signalID, vector<vector<int>*> iTrace, double baseline,int fw, int offset, TH2 *h_pe_adc, vector<int> *iPePix){
	int nx, ny;
	int nTrace=0;
	int iPixNum;
	double integ = 0;

	FindBin(signalID[0],&nx,&ny);
	//cout << " " << signalID[0] << " "<< nx << " "<<ny<<endl;
	for (int i = nx-1; i<=nx+3; i++){
		for(int j = ny-1; j<=ny+1; j++){
			if((i-1)>=0 && (i-1)<32 && (j-1)>=0 && (j-1)<16){
				iPixNum = FindPixelNumber(i-1, j-1);
				//cout<<i<<" "<<j<<" "<<iPixNum;
				//cout<<iPixNum<<endl;
				//cout << iTrace[iPixNum]->at(9) <<endl;
				h_trace_iso[nTrace]= new TH1F(TString::Format("h_trace_iso[%d]",nTrace),
				 	"Trace", iTrace[iPixNum]->size(),0,iTrace[iPixNum]->size());
				for (int k = 0; k<iTrace[iPixNum]->size();k++){
					if(iTrace[iPixNum]->at(k)-baseline>=0){
						h_trace_iso[nTrace]->Fill(k,iTrace[iPixNum]->at(k)-baseline);
					}
					else{
						h_trace_iso[nTrace]->Fill(k,0);
					}
				}
				integ = h_trace_iso[nTrace]->Integral(offset+1,fw+1,"width");
				if(integ == 0){
					integ = 0.1;
				}
				//cout<<" "<<"Integral: "<<integ<<endl;
				h_charge->SetBinContent(i,j,integ);
				if(iPePix->at(iPixNum)!= 0 && integ>0.1){
					h_pe_adc->Fill((double)iPePix->at(iPixNum),(double)integ);
				}
				
			}
			else{
				iPixNum = FindPixelNumber(i-1, j-1);
				h_trace_iso[nTrace]= new TH1F(TString::Format("h_trace_iso[%d]",nTrace),
				 	"Trace", iTrace[0]->size(),0,iTrace[0]->size());
				for (int k = 0; k<iTrace[0]->size();k++){
					h_trace_iso[nTrace]->Fill(k,0);
				}
				//cout<<endl;
			}
			
			
			nTrace++;
			integ = 0;
			//cout<<nTrace<<endl;
		}
	}

}

void ExtractSignalTracesFixed (TH2 *h_charge, TH1 *h_trace_iso[15], int *signalID, vector<vector<int>*> iTrace, double baseline,int fw, int offset, TH2 *h_pe_adc, vector<int> *iPePix, TProfile *p_dc){
	int nx, ny;
	int nTrace=0;
	int iPixNum;
	double integ = 0;

	FindBin(signalID[0],&nx,&ny);
	//cout << " " << signalID[0] << " "<< nx << " "<<ny<<endl;
	for (int i = nx-1; i<=nx+3; i++){
		for(int j = ny-1; j<=ny+1; j++){
			if((i-1)>=0 && (i-1)<32 && (j-1)>=0 && (j-1)<16){
				iPixNum = FindPixelNumber(i-1, j-1);
				//cout<<i<<" "<<j<<" "<<iPixNum;
				//cout<<iPixNum<<endl;
				//cout << iTrace[iPixNum]->at(9) <<endl;
				h_trace_iso[nTrace]= new TH1F(TString::Format("h_trace_iso[%d]",nTrace),
				 	"Trace", iTrace[iPixNum]->size(),0,iTrace[iPixNum]->size());
				for (int k = 0; k<iTrace[iPixNum]->size();k++){
					if(iTrace[iPixNum]->at(k)-baseline>=0){
						h_trace_iso[nTrace]->Fill(k,iTrace[iPixNum]->at(k)-baseline);
					}
					else{
						h_trace_iso[nTrace]->Fill(k,0);
					}
				}
				integ = h_trace_iso[nTrace]->Integral(offset+1,offset+1+6,"width");
				if(integ == 0){
					integ = 0.1;
				}
				//cout<<" "<<"Integral: "<<integ<<endl;
				h_charge->SetBinContent(i,j,integ);
				if(iPePix->at(iPixNum)!= 0 && integ>0.1){
					h_pe_adc->Fill((double)iPePix->at(iPixNum),(double)integ);
					p_dc->Fill((double)iPePix->at(iPixNum),(double)integ);
				}
				
			}
			else{
				iPixNum = FindPixelNumber(i-1, j-1);
				h_trace_iso[nTrace]= new TH1F(TString::Format("h_trace_iso[%d]",nTrace),
				 	"Trace", iTrace[0]->size(),0,iTrace[0]->size());
				for (int k = 0; k<iTrace[0]->size();k++){
					h_trace_iso[nTrace]->Fill(k,0);
				}
				//cout<<endl;
			}
			
			
			nTrace++;
			integ = 0;
			//cout<<nTrace<<endl;
		}
	}

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
	int iPixSignalID[2];
	int maxSampleID[2];
	int thresh = 14;
	double max_val[2];
	float totalPe;


	bool isNoise;
	bool isAbove;

	double dynamicRangeAGET = 750;//mV
	double mvPerPe = 1.92;
	double nBitsAGET = 12;
	int nPixIsol = 15;
	int full_width;
	double sum_base=0;
	double avg_base;
	double integ;
	int nSampleOffset=3;
	double px_rms;
	iTrace.assign(nPixels,0);
	int nNoise =0;
	int nSignal = 0;
	int nBelow = 0;


	int color_wheel [15] = {1,2,3,4,6,7,8,9,11,30,40,41,46,38,28};



	TFile *file = new TFile( "NSB_Traces_Merged_Fixed.root", "READ" );
	TH2 *h_camera = new TH2F ("h_camera", "Camera", 32,-0.5,31.5,16,-0.5,15.5);
	TH2 *h_camera_integ = new TH2F ("h_camera_integ", "Camera Charge", 32,-0.5,31.5,16,-0.5,15.5);
	TH1 *h_base = new TH1F("h_base", "Pedestal", 4096,0,4096);
	TH1 *h_rms = new TH1F("h_rms", "RMS Of Pedestal", nPixels,0,nPixels);
	TH1 *h_rms_dist = new TH1F("h_rms_dist","Camera RMS", 1000,0,10);
	TH2 *h_pe_adc = new TH2F("h_pe_adc","ADC per Pe",500,0,500,4096,0,4096);
	TH2 *h_pe_adc_fixed = new TH2F("h_pe_adc","ADC per Pe",500,0,500,4096,0,4096);
	TH1 *h_total_events = new TH1F("h_total_events","Total Events Bifocal",1002,9.9,1001);
	TH1 *h_total_events_acc = new TH1F("h_total_events_acc","Total Events Accepted",1002,9.9,1001);
	TProfile *p_dc = new TProfile("p_dc", "DC Profile", 300,0,300);
	TH1 *h_resolution = new TH1F("h_resolution","Resolution",300,0,300);
	TH1 *h_trace;
	TH1 *h_trace_iso[nPixIsol];
	TTree *tSimulatedEvents = (TTree*)file->Get( "Events/tSimulatedEvents" );
	TF1 *f_linear = new TF1("f_linear","[0]*x",0,4096);
	TTree *tEvent = (TTree*)file->Get( "Events/T0" );
	TGraph *gr_mean = new TGraph();

	for (int i = 0; i<nPixels; i++){
		tEvent->SetBranchAddress(TString::Format("vFADCTraces%i",i),&iTrace[i]);
	}
	tSimulatedEvents->SetBranchAddress("arrayTriggerBit",&isTriggered);
	tSimulatedEvents->SetBranchAddress("energy", &totalPe);
	tEvent->SetBranchAddress("vPEInPixel", &iPePix);
	tEvent->SetBranchAddress("vGroupsInTriggerCluster", &iTrigMusic);

	h_camera->SetStats(0);
	h_total_events_acc->SetStats(0);
	h_total_events->SetStats(0);
	h_resolution->SetStats(0);

	h_total_events_acc->Sumw2();
	h_total_events->Sumw2();
	h_resolution->Sumw2();

	p_dc->SetErrorOption("s");
	

	cEvent = 0;
	tSimulatedEvents->GetEntry(0);
	tEvent->GetEntry(0);

	ofstream file_NSB;
	file_NSB.open("NSBEvents.txt");

	/*while (!isTriggered){
		cEvent++;
		tSimulatedEvents->GetEntry(cEvent);
		tEvent->GetEntry(cEvent);
	}*/
	cout << cEvent <<endl;

	//tSimulatedEvents->GetEntry(9574);
	//tEvent->GetEntry(9574);

	bool isTrace = false;
	cout<<"N Events: "<<tSimulatedEvents->GetEntries()<<endl;
	for (int k = 0; k<tSimulatedEvents->GetEntries(); k++){
		tSimulatedEvents->GetEntry(k);
		tEvent->GetEntry(k);
		if(isTriggered){
			
			//cout<<k<<endl;

			h_base->Reset("ICES");
			h_camera_integ->Reset("ICES");
			h_camera->Reset("ICES");
			//cout<<iTrigMusic->size()<<endl;
			isNoise = isNSBTrigger(iTrigMusic,iTrace);
			if(!isNoise){
				h_total_events->Fill(totalPe);
				for (int l = 0; l<iTrigMusic->size(); l++){
					max_val[l]=0;
					maxSampleID[l]=0;
					for(int i=FindFirstPix(iTrigMusic->at(l)); i<FindFirstPix(iTrigMusic->at(l))+8; i++){
						
						sum_base=0;
						FindBin(i,&nx,&ny);
						//cout<<"Hello"<<endl;
						for (int j = 0; j<nSampleOffset-1; j++){
							sum_base += iTrace[i]->at(j);
							//px_rms += iTrace[i]->at(j)*iTrace[i]->at(j);
							gr_mean->SetPoint(nSampleOffset*i+j,i,iTrace[i]->at(j));
							h_base->Fill(iTrace[i]->at(j));
						}
						avg_base = sum_base/(nSampleOffset-1);
						//cout<<"Average: "<<avg_base<<endl;
	
						h_camera->SetBinContent(nx,ny,iPePix->at(i));
						
						integ = 0;
						//cout<<"Pixel ID: "<<i<<endl;
						if(!isTrace){
								h_trace = new TH1F("h_trace", "Trace", iTrace[i]->size(),0,iTrace[i]->size());
								isTrace=true;
						}
						
						for (int j = 0; j<iTrace[i]->size(); j++){
							/*if(iTrace[i]->at(j)-avg_base>0){
								integ += iTrace[i]->at(j)-avg_base;
							}*/
							if(max_val[l]<iTrace[i]->at(j)-avg_base){
								max_val[l]=iTrace[i]->at(j)-avg_base;
								maxSampleID[l]=j;
								iPixSignalID[l]=i;
							}
							//integ += iTrace[i]->at(j)-avg_base;
	
							//h_trace->SetBinContent(j+1,iTrace[i]->at(j)-avg_base);					
						}
	
						//cout<<"Integral: "<<integ<<endl;
	
						//h_camera_integ->SetBinContent(nx,ny,integ);
						//h_pe_adc->Fill(iPePix->at(i),integ);
	
						/*if(iPixSignal<integ){
							iPixSignal=integ;
							iPixSignalID[l]=i;
						}*/
					}
					//cout<<"Signal is in Pixel: "<<iPixSignalID<<endl;
					iPixSignal=0;
					

				}
				cout<<"Max Signal Pixels: "<< iPixSignalID[0]<< " "<<iPixSignalID[1]<<endl;
				isAbove = isAboveThreshold(max_val,thresh);
				if(isAbove){
					h_total_events_acc->Fill(totalPe);
					full_width = FindTraceWidth(maxSampleID,iPixSignalID,iTrace,nSampleOffset,avg_base);
					ExtractSignalTracesFixed(h_camera_integ,h_trace_iso, iPixSignalID, iTrace, avg_base,full_width,nSampleOffset, h_pe_adc_fixed, iPePix, p_dc);
					//ExtractSignalTraces(h_camera_integ,h_trace_iso, iPixSignalID, iTrace, avg_base,full_width,nSampleOffset, h_pe_adc, iPePix);
					nSignal++;
				}
				else{
					nBelow++;
				}
				
			}
			else{
				nNoise++;
				file_NSB<<k<<endl;
			
			}
			//IntegTrace(max_val,maxSampleID,iPixSignalID,iTrace,low_bd,up_bd);	
		}
	}

	/*for(int i=1; i<p_dc->GetNbinsX();i++){
		cout<<p_dc->GetBinError(i+1)/p_dc->GetBinContent(i+1)<<endl;
		h_resolution->SetBinContent(i+1,p_dc->GetBinError(i+1)/p_dc->GetBinContent(i+1));
	}
	
	//h_pe_adc->Fit("pol 1","Q","E");

	TCanvas *c_camera = new TCanvas("c_camera", "Camera Canvas",800,500);
	h_camera->Draw("colz");

	TCanvas *c_camera_charge = new TCanvas("c_camera_charge", "Camera Charge Canvas",800,500);
	h_camera_integ->Draw("colz");
	h_camera_integ->SetStats(0);

	TCanvas *c_traces_iso = new TCanvas("c_traces_iso","All Traces",800,500);
	for(int i=0;i<15;i++){
		h_trace_iso[i]->SetStats(0);
		//h_trace_iso[i]->SetError(0);
		if (i==0){
			h_trace_iso[i]->Draw("HIST CP");
			h_trace_iso[i]->GetXaxis()->SetTitle("Sample [10ns]");
			h_trace_iso[i]->GetYaxis()->SetTitle("DC");
			h_trace_iso[i]->GetYaxis()->SetRangeUser(0,4096);
			h_trace_iso[i]->SetLineColor(color_wheel[i]);
		}
		else{
			h_trace_iso[i]->Draw("SAME HIST CP");
			h_trace_iso[i]->SetLineColor(color_wheel[i]);
		}
	}

	h_total_events_acc->Divide(h_total_events);

	/*TCanvas *c_mean = new TCanvas("c_mean", "Mean Dispersions", 800,500);
	gr_mean->GetHistogram()->Draw();
	gr_mean->GetHistogram()->GetXaxis()->SetTitle("Pixel Number");
	gr_mean->GetHistogram()->GetYaxis()->SetTitle("Pedestal");
	gr_mean->Draw("P");

	TCanvas *c_rms = new TCanvas("c_rms", "Canvas RMS", 800,500);
	h_rms->Draw("HIST");
	h_rms->SetStats(0);
	h_rms->GetXaxis()->SetTitle("Pixel");
	h_rms->GetYaxis()->SetTitle("RMS");
	h_rms->SetStats(0);
	h_rms->SetMarkerStyle(kFullCircle);*/

	/*TCanvas *c_base = new TCanvas("c_base", "Canvas Baseline", 800,500);
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
	h_pe_adc->GetXaxis()->SetRangeUser(0,400);
	h_pe_adc->GetXaxis()->SetTitle("Pe");
	h_pe_adc->GetYaxis()->SetTitle("Charge");
	h_pe_adc->SetTitle("Baseline Determined Window");
	h_pe_adc->Fit(f_linear);
	h_pe_adc->Draw("SCAT");
	f_linear->Draw("SAME");

	TCanvas *c_pe_chg_fixed = new TCanvas("c_pe_chg_fixed","Canvas Pe Charge",800,500);
	h_pe_adc_fixed->GetXaxis()->SetRangeUser(0,400);
	h_pe_adc_fixed->GetXaxis()->SetTitle("Pe");
	h_pe_adc_fixed->GetYaxis()->SetTitle("Charge");
	h_pe_adc_fixed->SetTitle("Fixed Integration Window");
	h_pe_adc_fixed->Fit(f_linear);
	h_pe_adc_fixed->Draw("SCAT");
	f_linear->Draw("SAME");

	TCanvas *c_efficiency = new TCanvas("c_efficiency", "Canvas Efficienty",800,500);
	c_efficiency->SetLogx();
	h_total_events_acc->Draw("PE");
	h_total_events_acc->Sumw2();
	h_total_events_acc->GetXaxis()->SetTitle("Total Photoelectron Signal");
	h_total_events_acc->GetYaxis()->SetTitle("Efficiency");

	TCanvas *c_resolution = new TCanvas("c_resolution", "Resolution", 800,500);
	h_resolution->Draw("PE");
	h_resolution->GetXaxis()->SetTitle("Photoelectrons");
	h_resolution->GetYaxis()->SetTitle("Resolution");

	TCanvas *c_dcProfile = new TCanvas("c_dcProfile", "DC Profile",800,500);
	p_dc->Draw("PE");
	p_dc->GetXaxis()->SetTitle("Photoelectrons");
	p_dc->GetYaxis()->SetTitle("Charge");*/

	cout<<"NSB Events: "<<nNoise<<"\nSignal Events: "<<nSignal<<endl;
	cout<<"Threshold Events: "<<nBelow<<endl;
	cout<<"Dscarded Events \%: "<< (double) 100.0*(nNoise+nBelow)/(double)(nNoise+nBelow+nSignal)<<endl; 
	/*TCanvas *c_trace = new TCanvas("c_trace","Trace",800,500);
	h_trace->Draw("HIST");
	h_trace->GetXaxis()->SetTitle("Sample [10 ns]");
	h_trace->GetYaxis()->SetTitle("DC");*/
	file_NSB.close();
	return 0;
}
