#include <iostream>
#include <iterator>
#include <TH1F.h>
#include <TTimer.h>
#include <TGraph.h>
#include <TMath.h>
#include <TProfile.h>
#include <vector>

using namespace std;
using namespace std::chrono;

void FindNeighborPixels();
void LoadSPB2Events(string NameofFile);
void FindTrigEvents();
void PlotBaseline();
void PlotBaselineMeanRMS();
void PlotDCPE();
void PlotResolution();
void CalcEfficiency(int threshold);
void PlotEfficiencyPE();
void PlotSingleEvent();

float heffPE[8][20];
float SoftwarePE[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

// These are trees and parameters loaded from event file
TTree *tSimulatedEvents;
TTree *T0;
Float_t energy;
UInt_t eventNumber;
Float_t xcore;
Float_t ycore;
Bool_t arrayTriggerBit;
vector<int> *vTriggerCluster=0;
vector< vector<Int_t> *> iFADCTraceInPixel;
vector<Int_t> *iPEInPixel;

// Signal extraction paramters
const int iNumPixels = 512; 
const int SignalStart = 4;
const int SignalEnd = 10;
const int SignalWidth = SignalEnd - SignalStart;
const int CoincWindow = 2;
int PEThreshold = 0;
float DCPEFactor = 0;
double Baseline[iNumPixels];
double PixelCharge[iNumPixels];
vector<int> TriggeredEventsID;
vector<vector<int> > vFiredPixels;
vector<vector<double> > BaselineDist;
TProfile* DCProfile;

/*	// This is the calculated baseline for all pixels from NSB triggered events
double Baseline[iNumPixels] = { 499.502, 499.537, 499.435, 499.475, 499.525, 499.475, 499.506, 499.526, 499.564, 499.531, 499.512, 499.522, 499.534, 499.516, 499.538, 499.544, 499.515, 499.468, 499.483, 499.494,
								499.512, 499.566, 499.521, 499.494, 499.496, 499.481, 499.574, 499.488, 499.497, 499.528, 499.502, 499.5, 499.526, 499.464, 499.584, 499.527, 499.526, 499.524, 499.552, 499.529,
								499.497, 499.542, 499.546, 499.496, 499.566, 499.486, 499.534, 499.517, 499.525, 499.444, 499.524, 499.497, 499.511, 499.565, 499.555, 499.472, 499.531, 499.492, 499.494, 499.552,
								499.554, 499.506, 499.56, 499.481, 499.527, 499.488, 499.526, 499.529, 499.513, 499.469, 499.563, 499.522, 499.52, 499.548, 499.52, 499.477, 499.558, 499.502, 499.491, 499.513,
								499.481, 499.528, 499.556, 499.513, 499.465, 499.456, 499.551, 499.547, 499.524, 499.523, 499.541, 499.538, 499.472, 499.522, 499.513, 499.536, 499.471, 499.515, 499.575, 499.543,
								499.569, 499.498, 499.552, 499.545, 499.503, 499.557, 499.486, 499.48, 499.527, 499.543, 499.5, 499.499, 499.539, 499.525, 499.527, 499.511, 499.488, 499.527, 499.502, 499.501,
								499.518, 499.461, 499.533, 499.532, 499.484, 499.5, 499.565, 499.535, 499.482, 499.53, 499.492, 499.503, 499.549, 499.519, 499.467, 499.497, 499.469, 499.541, 499.545, 499.493,
								499.475, 499.516, 499.523, 499.573, 499.492, 499.519, 499.512, 499.517, 499.481, 499.534, 499.535, 499.539, 499.454, 499.495, 499.56, 499.459, 499.482, 499.492, 499.548, 499.484,
								499.595, 499.485, 499.537, 499.558, 499.509, 499.488, 499.546, 499.519, 499.543, 499.508, 499.575, 499.462, 499.526, 499.548, 499.587, 499.53, 499.462, 499.531, 499.473, 499.523,
								499.518, 499.5, 499.481, 499.508, 499.507, 499.509, 499.446, 499.509, 499.502, 499.469, 499.539, 499.523, 499.522, 499.483, 499.534, 499.537, 499.487, 499.508, 499.553, 499.54,
								499.52, 499.492, 499.552, 499.506, 499.514, 499.492, 499.552, 499.515, 499.525, 499.518, 499.507, 499.464, 499.534, 499.501, 499.496, 499.464, 499.551, 499.507, 499.532, 499.496,
								499.555, 499.537, 499.527, 499.51, 499.579, 499.506, 499.49, 499.523, 499.503, 499.557, 499.511, 499.559, 499.542, 499.464, 499.473, 499.491, 499.545, 499.557, 499.52, 499.517,
								499.526, 499.541, 499.5, 499.501, 499.49, 499.523, 499.564, 499.533, 499.495, 499.412, 499.522, 499.509, 499.505, 499.522, 499.476, 499.487, 499.468, 499.519, 499.491, 499.511,
								499.567, 499.51, 499.482, 499.536, 499.543, 499.531, 499.5, 499.568, 499.529, 499.521, 499.571, 499.494, 499.55, 499.537, 499.517, 499.516, 499.516, 499.502, 499.496, 499.509,
								499.489, 499.525, 499.526, 499.555, 499.516, 499.467, 499.487, 499.512, 499.469, 499.526, 499.495, 499.566, 499.489, 499.495, 499.511, 499.536, 499.473, 499.527, 499.454, 499.519,
								499.522, 499.517, 499.522, 499.553, 499.507, 499.549, 499.528, 499.542, 499.488, 499.557, 499.501, 499.504, 499.566, 499.539, 499.513, 499.541, 499.497, 499.525, 499.546, 499.532,
								499.462, 499.506, 499.511, 499.483, 499.497, 499.497, 499.504, 499.536, 499.484, 499.493, 499.519, 499.565, 499.531, 499.525, 499.491, 499.497, 499.534, 499.525, 499.463, 499.55,
								499.493, 499.491, 499.515, 499.524, 499.511, 499.498, 499.521, 499.567, 499.537, 499.572, 499.515, 499.519, 499.461, 499.528, 499.547, 499.464, 499.543, 499.559, 499.527, 499.509,
								499.572, 499.455, 499.505, 499.548, 499.559, 499.552, 499.502, 499.531, 499.537, 499.486, 499.506, 499.539, 499.524, 499.566, 499.563, 499.506, 499.49, 499.503, 499.459, 499.477,
								499.511, 499.533, 499.463, 499.524, 499.505, 499.523, 499.538, 499.55, 499.488, 499.537, 499.509, 499.5, 499.559, 499.533, 499.517, 499.514, 499.513, 499.523, 499.465, 499.6,
								499.509, 499.617, 499.542, 499.528, 499.576, 499.56, 499.532, 499.527, 499.512, 499.533, 499.471, 499.501, 499.513, 499.498, 499.52, 499.556, 499.478, 499.502, 499.506, 499.485,
								499.585, 499.569, 499.542, 499.563, 499.516, 499.529, 499.502, 499.538, 499.513, 499.543, 499.542, 499.547, 499.54, 499.507, 499.5, 499.481, 499.553, 499.467, 499.547, 499.559,
								499.514, 499.515, 499.524, 499.528, 499.506, 499.488, 499.535, 499.49, 499.548, 499.486, 499.509, 499.501, 499.493, 499.523, 499.524, 499.545, 499.529, 499.5, 499.446, 499.535,
								499.564, 499.536, 499.481, 499.544, 499.491, 499.544, 499.515, 499.426, 499.544, 499.496, 499.521, 499.52, 499.537, 499.481, 499.51, 499.557, 499.535, 499.547, 499.543, 499.528,
								499.495, 499.506, 499.532, 499.521, 499.508, 499.47, 499.533, 499.49, 499.551, 499.522, 499.495, 499.532, 499.497, 499.496, 499.524, 499.506, 499.529, 499.523, 499.504, 499.508,
								499.5, 499.575, 499.503, 499.504, 499.458, 499.469, 499.482, 499.5, 499.516, 499.512, 499.531, 499.556};
*/

// Histograms and Display parameters
int iLastPix = -1;
TLatex *text = 0;
TCanvas *cDisplay = 0;
TCanvas *eDisplay = 0;
TCanvas *gDisplay = 0;
TCanvas *iDisplay = 0;
TH1F *hPixelTrace =0;
TH1F *hBaseDist =0;


int main(){

	cDisplay = new TCanvas("cDisplay","Baseline",4000,2000);
	cDisplay->Divide(2,2);

	eDisplay = new TCanvas("eDisplay","Signal",4000,2000);
	eDisplay->Divide(2,2);

	gDisplay = new TCanvas("gDisplay","Resolution",4000,2000);
	gDisplay->Divide(2,2);

	iDisplay = new TCanvas("iDisplay","Efficiency",4000,2000);
	iDisplay->Divide(2,2);

	string FileName = "/home/mahdi/Programs/SPB2/SPB2_CT_Simulations/data/test_50.root";
	//string FileName = "/storage/hive/project/phy-otte/shared/Merged_SPB2_CARE_NSB_Sims/NSB_Traces_Merged_Manual.root";
	//string FileName = "/mnt/hgfs/Shared Folder/NSB_Traces_Merged_Manual.root";

	LoadSPB2Events(FileName);

	FindTrigEvents();
	PlotBaseline();
	PlotBaselineMeanRMS();

	PlotDCPE();
	PlotResolution();

	CalcEfficiency(10);
	FindNeighborPixels();
	PlotSingleEvent();
}

Bool_t HandleInput()
{
	TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
	while (1)
	{
		timer.TurnOn();
		TString input = Getline("Type 'q' to exit, <return> to go on:");
		timer.TurnOff();

		if (input=="q\n")
			return kFALSE;

		if (input=="\n")
			return kTRUE;
	}

	return kFALSE;
}

void LoadSPB2Events(string NameofFile)
{
	TFile *fO = new TFile(NameofFile.c_str(), "READ");

	tSimulatedEvents = (TTree*)fO->Get("Events/tSimulatedEvents");
	tSimulatedEvents->SetBranchAddress("energy",&energy);
	tSimulatedEvents->SetBranchAddress("xcore",&xcore);
	tSimulatedEvents->SetBranchAddress("ycore",&ycore);
	tSimulatedEvents->SetBranchAddress("arrayTriggerBit",&arrayTriggerBit);
	tSimulatedEvents->SetBranchAddress("eventNumber",&eventNumber);

	iFADCTraceInPixel.assign(iNumPixels,0);

	T0 = (TTree*)fO->Get("Events/T0");
	T0->SetBranchAddress("vGroupsInTriggerCluster",&vTriggerCluster);
	T0->SetBranchAddress("vPEInPixel", &iPEInPixel);
	TString name;
	for(int g=0;g<iNumPixels;g++)
	{
		name.Form("vFADCTraces%i",g);
		T0->SetBranchAddress(name,&(iFADCTraceInPixel[g]));
	}
}

void FindBin(int iPix,int *nx, int *ny)
{
	int MUSIC_ID = iPix/8;
	int MUSIC_Channel = iPix%8; 
	int MUSIC_row = MUSIC_ID/16;
	int MUSIC_column = MUSIC_ID%16;

	*ny = 4*MUSIC_row+MUSIC_Channel%4+1; 
	*nx = 2*MUSIC_column+MUSIC_Channel/4+1; 
}

int FindPixel(int x, int y)
{
	int MUSIC_column = x/2;
	int MUSIC_row = y/4;
	int MUSIC_ID = MUSIC_column+MUSIC_row*16;
	int MUSIC_Channel = y%4+4*(x%2);

	return MUSIC_row*8*16+MUSIC_column*8+MUSIC_Channel;
}

int FindPix(int nx, int ny)
{
	int PixelID;
	PixelID = ((nx-1)*4+((ny-1)%4)+((ny-1)/4)*128);
	return PixelID;
}

void FindNeighborPixels()
{
	for (int k=0; k<iNumPixels; k++)
	{
		int nx, ny;
		FindBin(k, &nx, &ny);
		vector<int> vNeighbor;
		for (int i=nx-1; i<nx+4 && i<33; i++)
		{
			for (int j=ny-1; j<ny+2 && j<17; j++)
			{
				if((i > 0) && (j > 0))
				{
					int tmp = FindPix(i, j);
					if((tmp < 512) && (tmp > -1))
						vNeighbor.push_back(tmp);
				}
			}
		}
		vFiredPixels.push_back(vNeighbor);			
	}
}

void DrawMUSICBoundaries()
{
	TBox *b = new TBox(-.5,-0.5,1.5,3.5);
	b->SetFillStyle(0);
	b->SetLineColor(kRed);
	b->Draw();
	TPoint p;
	for(int i=0;i<iNumPixels/8;i++)
	{
		TBox *bn = (TBox*)b->Clone();
		bn->SetX1((i%16)*2-0.5);
		bn->SetX2((i%16)*2+1.5);
		bn->SetY1((i/16)*4-0.5);
		bn->SetY2((i/16)*4+3.5);
		bn->Draw();
	}
}

void ShowInfoAtCursor(int x, int y)
{
	int MUSIC_column = x/2;
	int MUSIC_row = y/4;
	int MUSIC_ID = MUSIC_column+MUSIC_row*16;
	int MUSIC_Channel = y%4+4*(x%2);
	int PixID = MUSIC_row*8*16+MUSIC_column*8+MUSIC_Channel;

	TString statusline;
	statusline.Form("MUSIC_ID: %i MUSIC_Channel: %i Pixel ID: %i  PEs in Pixel: %i",MUSIC_ID,MUSIC_Channel,PixID,iPEInPixel->at(PixID));
	if(text!=0)
		text->Delete();
	TLatex T1;
	text = T1.DrawLatexNDC(0.1,0.95,statusline.Data());
	gPad->Modified();
	gPad->Update();
}

void FindTrigEvents()
{
	int TotalEvents = tSimulatedEvents->GetEntries();
	for (int i=0; i<TotalEvents; i++)
	{
		tSimulatedEvents->GetEntry(i);
		if(arrayTriggerBit)
		{
			TriggeredEventsID.push_back(i);
		}
	}
	cout<<"Total number of events: "<<TotalEvents<<endl;
	cout<<"Total triggered events: "<<TriggeredEventsID.size()<<endl;
}

void PlotBaseline()
{
	int BaselineWidth = 4;
	float SingleMFactor = (1.0/BaselineWidth);
	float TotalMFactor = 0;

	cout<<"Calculating the baseline distribution"<<endl;

	for (int i=0; i<TriggeredEventsID.size(); i++)
	{
		tSimulatedEvents->GetEntry(TriggeredEventsID[i]);
		T0->GetEntry(TriggeredEventsID[i]);
		vector<double> tmp;
		for(int j=0;j<iNumPixels;j++)
		{
			double tmpbase;
			for (int k=0; k<BaselineWidth; k++)
			{
				tmpbase += ((iFADCTraceInPixel[j])->at(k));
			}
			tmp.push_back(tmpbase*SingleMFactor);
			Baseline[j] += tmpbase;
			tmpbase = 0;
		}
		BaselineDist.push_back(tmp);
	}

	TotalMFactor = (1.0/(BaselineWidth*TriggeredEventsID.size()));

	for(int i=0;i<iNumPixels;i++)
	{
		Baseline[i] = Baseline[i]*TotalMFactor;
	}


	cDisplay->cd(1);
	gPad->AddExec("ex","BaselinePixelClicked()");
	TH2F *h1Display = new TH2F("h1Display","Average Baseline",32,-0.5,31.5,16,-0.5,15.5);
	h1Display->SetStats(0);
	h1Display->Draw("colz");
	DrawMUSICBoundaries();
	h1Display->Clear();

	for(int g=0;g<iNumPixels;g++)
	{
		int nx, ny;
		FindBin(g,&nx,&ny);
		h1Display->SetBinContent(nx,ny,Baseline[g]);
	}

	cDisplay->cd(1)->Modified();
	cDisplay->cd(1)->Update();
	cout<<"Finished the baseline calculation"<<endl;
}

void PlotBaselineMeanRMS()
{
	TH1F *hBaseTemp = new TH1F("hBaseTemp","Temporary Histogram for Baseline",100,495.0,505.0);

	TH1F *hBaseMean = new TH1F("hBaseMean","Pixel Baseline Mean",50,499.0,500.0);
	hBaseMean->SetStats(0);
	hBaseMean->GetXaxis()->SetTitle("Mean");
	hBaseMean->GetYaxis()->SetTitle("Mean Frequency");

	TH1F *hBaseRMS = new TH1F("hBaseRMS","Pixel Baseline RMS",50,1.0,1.5);
	hBaseRMS->SetStats(0);
	hBaseRMS->GetXaxis()->SetTitle("RMS");
	hBaseRMS->GetYaxis()->SetTitle("RMS Frequency");

	for (int i=0; i<iNumPixels; i++)
	{
		for(int s=0;s<BaselineDist.size();s++)
			{
				hBaseTemp->Fill(BaselineDist[s][i]);
			}
		hBaseRMS->Fill(hBaseTemp->GetRMS());
		hBaseMean->Fill(hBaseTemp->GetMean());
		hBaseTemp->Reset();
	}
	
	cDisplay->cd(2);
	hBaseMean->Draw();
	gPad->Modified();
	gPad->Update();

	cDisplay->cd(4);
	hBaseRMS->Draw();
	gPad->Modified();
	gPad->Update();
}

void ShowBaseDist(int iPix)
{
	if(hBaseDist==0)
		hBaseDist = new TH1F("hBaseDist","Pixel Baseline Dist.",80,490.0,510.0);
		hBaseDist->SetStats(0);
		hBaseDist->GetXaxis()->SetTitle("Baseline");
		hBaseDist->GetYaxis()->SetTitle("Baseline Frequency");

	hBaseDist->Reset();
	TString title;
	title.Form("Baseline Dist. of Pixel %i",iPix);
	hBaseDist->SetTitle(title);
	for(int s=0;s<BaselineDist.size();s++)
	{
		hBaseDist->Fill(BaselineDist[s][iPix]);
	}

	cDisplay->cd(3);
	hBaseDist->Draw();
	gPad->Modified();
	gPad->Update();
}

void BaselinePixelClicked()
{
	int event = gPad->GetEvent();
	TObject *o = gPad->GetSelected();
	if (!o) return;
	if (!(o->InheritsFrom("TH2")))
	   return;
	TH2F *h = (TH2F*)o;
	int px = gPad->GetEventX();
	int py = gPad->GetEventY();
	Float_t xx = gPad->AbsPixeltoX(px);
	Float_t yy = gPad->AbsPixeltoY(py);
	Float_t x = 0.5+gPad->PadtoX(xx);
	Float_t y = 0.5+gPad->PadtoY(yy);

	int pix = FindPixel((int)x,(int)y);
	if (event == 11)
		ShowBaseDist(pix);
}

void PlotDCPE()
{
	vector<int> DCValue;
	vector<int> PEValue;
	int BeginPixelID =0;
	int TotalCharge =0;

	// This section calculates the DC and PE for each event and its fired pixels
	for (int i=0; i<TriggeredEventsID.size(); i++)
	{
		tSimulatedEvents->GetEntry(TriggeredEventsID[i]);
		T0->GetEntry(TriggeredEventsID[i]);
		BeginPixelID = vTriggerCluster->at(0)*8;
		for(int j=BeginPixelID;j<BeginPixelID+16;j++)
		{
			if((iPEInPixel->at(j)) != 0)
			{
				PEValue.push_back(iPEInPixel->at(j));
				for(int k=SignalStart; k<SignalEnd; k++)
				{
					TotalCharge += ((iFADCTraceInPixel[j])->at(k));
				}
				TotalCharge = TotalCharge - (SignalWidth*Baseline[j]);
				DCValue.push_back(TotalCharge);
				TotalCharge=0;
			}
		}
	}

	// This section fills up a 2D histogram with DC and PE and finds the ratio by fitting it to f(x) = x
	gDisplay->cd(1);
	TH2D *hScatDCPE = new TH2D("hScatDCPE","DC vs. PE Histogram",100,-20,300,100,-20,600);
	hScatDCPE->GetXaxis()->SetTitle("PE");
	hScatDCPE->GetYaxis()->SetTitle("DC");
	hScatDCPE->SetStats(0);
	for (int i=0; i<PEValue.size(); i++)
	{
	  hScatDCPE->Fill(PEValue[i], DCValue[i]);
	}
	TF1 *fDCPE = new TF1("fDCPE","[0]*x",0,300);
	fDCPE->SetParNames("slope");
	hScatDCPE->Fit(fDCPE);
	hScatDCPE->Draw("SCAT");
	DCPEFactor = (1.0/fDCPE->GetParameter(0));
	gDisplay->cd(1)->Modified();
	gDisplay->cd(1)->Update();

	// This section creates a TProfile to calculate the standard deviation of the mean of DC values.
	int NumofBins = 50;
	gDisplay->cd(2);
	DCProfile = new TProfile("DCProfile","Profile of extracted digital counts",2*NumofBins,0.0,300.0);
	DCProfile->GetXaxis()->SetTitle("PE");
	DCProfile->GetYaxis()->SetTitle("DC");
	DCProfile->SetStats(0);
	for (int i=0; i<PEValue.size(); i++)
	{
	  DCProfile->Fill(PEValue[i], DCValue[i]);
	}
	DCProfile->SetErrorOption("s");
	DCProfile->Draw();
	gDisplay->cd(2)->Modified();
	gDisplay->cd(2)->Update();
}

void PlotResolution()
{
	int NumofBins = 100;
	gDisplay->cd(3);
	TH1F *hResolution = new TH1F("hResolution","Resolution vs. PE",NumofBins,0.0,300.0);
	hResolution->GetXaxis()->SetTitle("PE");
	hResolution->GetYaxis()->SetTitle("Resolution");
	hResolution->GetYaxis()->SetRange(0.0,2.0);
	hResolution->SetStats(0);
	for (int i=0; i<NumofBins; i++)
	{
		float ResValue = (DCProfile->GetBinError(i+1))/(DCProfile->GetBinContent(i+1));
		if (ResValue > 0)
			hResolution->SetBinContent(i+1, ResValue);
	}
	hResolution->Draw();

	TF1 *fa1 = new TF1("fa1","1/sqrt(x)",0.1,300);
	fa1->Draw("same");
	gDisplay->cd(3)->Modified();
	gDisplay->cd(3)->Update();
}

void ShowPixeltrace(int iPix)
{
	if(hPixelTrace==0)
	{
		hPixelTrace = new TH1F("hPixelTrace","Pixel Trace",25,-0.5,24.5);
		hPixelTrace->SetStats(0);
		hPixelTrace->GetXaxis()->SetTitle("ADC sample");
		hPixelTrace->GetYaxis()->SetTitle("ADC counts");
	}

	hPixelTrace->Clear();
	TString title;
	title.Form("Trace of Pixel %i",iPix);
	hPixelTrace->SetTitle(title);
	for(int s=0;s<iFADCTraceInPixel[iPix]->size();s++)
	{
		hPixelTrace->SetBinContent(s+1,(iFADCTraceInPixel[iPix])->at(s));
	}

	eDisplay->cd(3);
	hPixelTrace->Draw();
	gPad->Modified();
	gPad->Update();
}

void PixelClicked()
{
	//this action function is called whenever you move the mouse
	//it just prints the id of the picked triangle
	//you can add graphics actions instead
	int event = gPad->GetEvent();
	TObject *o = gPad->GetSelected();
	if (!o) return;
	if (!(o->InheritsFrom("TH2")))
	   return;
	TH2F *h = (TH2F*)o;
	int px = gPad->GetEventX();
	int py = gPad->GetEventY();
	Float_t xx = gPad->AbsPixeltoX(px);
	Float_t yy = gPad->AbsPixeltoY(py);
	Float_t x = 0.5+gPad->PadtoX(xx);
	Float_t y = 0.5+gPad->PadtoY(yy);
	int pix = FindPixel((int)x,(int)y);
	if(pix!=iLastPix)
	{
		iLastPix=pix;
		ShowInfoAtCursor((int) x, (int) y);
	}

	if (event == 11)
		ShowPixeltrace(pix);
}

void CalcEfficiency(int threshold)
{
	int NumofBins = 50;
	int FirstPixelID = 0;
	int FiredPixelID =0;
	int ReconstructedEvents = 0;
	float AlgorithmEfficiency =0;

	PEThreshold = threshold;

	TH1F *hEfficiency = new TH1F("hEfficiency","Reconstruction Efficiency",2*NumofBins,1,1001);
	hEfficiency->GetXaxis()->SetTitle("Total Shower Photoelectron Signal");
	hEfficiency->GetYaxis()->SetTitle("Reconstruction Efficiency");
	hEfficiency->SetStats(0);
	hEfficiency->SetLineWidth(2);
	hEfficiency->Sumw2();		

	TH1F *hTriggered = new TH1F("hTriggered","Trigger Efficiency",2*NumofBins,1,1001);
	hTriggered->GetXaxis()->SetTitle("Total Shower Photoelectron Signal");
	hTriggered->GetYaxis()->SetTitle("Trigger Efficiency");
	hTriggered->SetStats(0);
	hTriggered->SetLineWidth(2);
	hTriggered->Sumw2();		
	
	TH1F *hTotalEfficiency = new TH1F("hTotalEfficiency","Overall Efficiency",2*NumofBins,1,1001);
	hTotalEfficiency->GetXaxis()->SetTitle("Total Shower Photoelectron Signal");
	hTotalEfficiency->GetYaxis()->SetTitle("Overall Efficiency");
	hTotalEfficiency->SetStats(0);
	hTotalEfficiency->SetLineWidth(2);
	hTotalEfficiency->Sumw2();		

	for(int n=0;n<TriggeredEventsID.size();n++)
	{
		FiredPixelID = -1;
		tSimulatedEvents->GetEntry(TriggeredEventsID[n]);
		T0->GetEntry(TriggeredEventsID[n]);

		hTriggered->Fill(energy);

		FirstPixelID = vTriggerCluster->at(0)*8;
		for(int i=FirstPixelID; i<FirstPixelID+16; i++)
		{
			for(int j=SignalStart; j<SignalEnd; j++)
				{
					PixelCharge[i] += ((iFADCTraceInPixel[i])->at(j));
				}
		PixelCharge[i] = PixelCharge[i] - (SignalWidth*Baseline[i]);
		}

		vector<int> vPixCandidate;
		vector<int> vPeakTimeIndex;
		for(int i=FirstPixelID; i<FirstPixelID+8; i++)
		{
			int M1PeakTimeIndex =0;
			int M2PeakTimeIndex =0;
			float M1PeakValue = 0;
			float M2PeakValue = 0;
			for(int j=SignalStart; j<SignalEnd; j++)
			{
				if((iFADCTraceInPixel[i]->at(j) - Baseline[i]) > M1PeakValue)
				{
					M1PeakTimeIndex = j;
					M1PeakValue = iFADCTraceInPixel[i]->at(j) - Baseline[i];
				}

				if((iFADCTraceInPixel[i+8]->at(j) - Baseline[i+8]) > M2PeakValue)
				{
					M2PeakTimeIndex = j;
					M2PeakValue = iFADCTraceInPixel[i+8]->at(j) - Baseline[i+8];
				}
			}

			if((abs(M2PeakTimeIndex - M1PeakTimeIndex) < CoincWindow) && (M1PeakValue > PEThreshold) && (M2PeakValue > PEThreshold))
			{
				vPixCandidate.push_back(i);
				vPixCandidate.push_back(i+8);
				vPeakTimeIndex.push_back(M1PeakTimeIndex);
				vPeakTimeIndex.push_back(M2PeakTimeIndex);
			}
		}

		int PixMaxCharge = 0;
		int PeakTimeIndex = 0;
		if(vPixCandidate.size() > 0)
		{
			for (int i=0; i<vPixCandidate.size(); i++)
			{
				if(PixelCharge[vPixCandidate[i]] > PixMaxCharge)
				{
					PeakTimeIndex = vPeakTimeIndex[i];
					FiredPixelID = vPixCandidate[i];
					PixMaxCharge = PixelCharge[vPixCandidate[i]];
				}
			}

			if (FiredPixelID > -1)
			{
				hEfficiency->Fill(energy);
				ReconstructedEvents++;
			}
		}
		else
		{
			//cout<<"Did not meet bi-focal requirements: "<<TriggeredEventsID[n]<<endl;
		}

		for(int g=0; g<iNumPixels; g++)
		{
			PixelCharge[g] = 0;
		}
	}

	iDisplay->cd(1);
	hEfficiency->Divide(hTriggered);
	hEfficiency->Draw();

	iDisplay->cd(2);
	hTriggered->Scale(1e-3*(4*105*55)/(199.8*99.8));
	hTriggered->Draw();

	iDisplay->cd(3);
	hTotalEfficiency->Multiply(hEfficiency,hTriggered);
	hTotalEfficiency->Draw();

	cout<<"Total Reconstructed Events: "<<ReconstructedEvents<<endl;
	AlgorithmEfficiency = (ReconstructedEvents*100.0/TriggeredEventsID.size());
	cout<<"Your Algorithm Efficiency is: "<<AlgorithmEfficiency<<"%"<<endl;
	
	heffPE[0][threshold-1] = hEfficiency->GetBinContent(3)*100;
	heffPE[1][threshold-1] = hEfficiency->GetBinContent(4)*100;
	heffPE[2][threshold-1] = hEfficiency->GetBinContent(7)*100;
	heffPE[3][threshold-1] = hEfficiency->GetBinContent(10)*100;
	heffPE[4][threshold-1] = hEfficiency->GetBinContent(16)*100;
	heffPE[5][threshold-1] = hEfficiency->GetBinContent(26)*100;
	heffPE[6][threshold-1] = hEfficiency->GetBinContent(40)*100;
	heffPE[7][threshold-1] = hEfficiency->GetBinContent(63)*100;

	//hEfficiency->Reset();
	//hTriggered->Reset();
}

void PlotEfficiencyPE()
{
	for (int i=1; i<21; i++)
	{
		CalcEfficiency(i);
	}

	iDisplay->cd(4);
	TGraph *RecEfficiency1 = new TGraph(20, SoftwarePE, heffPE[0]);
	TGraph *RecEfficiency2 = new TGraph(20, SoftwarePE, heffPE[1]);
	TGraph *RecEfficiency3 = new TGraph(20, SoftwarePE, heffPE[2]);
	TGraph *RecEfficiency4 = new TGraph(20, SoftwarePE, heffPE[3]);
	TGraph *RecEfficiency5 = new TGraph(20, SoftwarePE, heffPE[4]);
	TGraph *RecEfficiency6 = new TGraph(20, SoftwarePE, heffPE[5]);
	TGraph *RecEfficiency7 = new TGraph(20, SoftwarePE, heffPE[6]);
	TGraph *RecEfficiency8 = new TGraph(20, SoftwarePE, heffPE[7]);

	RecEfficiency1->SetMarkerStyle(20);
	RecEfficiency2->SetMarkerStyle(21);
	RecEfficiency3->SetMarkerStyle(22);
	RecEfficiency4->SetMarkerStyle(23);
	RecEfficiency5->SetMarkerStyle(29);
	RecEfficiency6->SetMarkerStyle(33);
	RecEfficiency7->SetMarkerStyle(43);
	RecEfficiency8->SetMarkerStyle(47);
	
	RecEfficiency1->SetLineColor(1);
	RecEfficiency2->SetLineColor(2);
	RecEfficiency3->SetLineColor(3);
	RecEfficiency4->SetLineColor(4);
	RecEfficiency5->SetLineColor(5);
	RecEfficiency6->SetLineColor(6);
	RecEfficiency7->SetLineColor(7);
	RecEfficiency8->SetLineColor(9);

	TMultiGraph  *RecEfficiency  = new TMultiGraph();
	RecEfficiency->Add(RecEfficiency1);
	RecEfficiency->Add(RecEfficiency2);
	RecEfficiency->Add(RecEfficiency3);
	RecEfficiency->Add(RecEfficiency4);
	RecEfficiency->Add(RecEfficiency5);
	RecEfficiency->Add(RecEfficiency6);
	RecEfficiency->Add(RecEfficiency7);
	RecEfficiency->Add(RecEfficiency8);
	RecEfficiency->Draw("ALP");
	RecEfficiency->GetXaxis()->SetTitle("Software PE threshold");
	RecEfficiency->GetYaxis()->SetTitle("Reconstruction Efficiency (%)");
	gPad->Modified();
	gPad->Update();
}

void PlotSingleEvent()
{
	eDisplay->cd(1);
	gPad->AddExec("ex","PixelClicked()");
	TH2F *h3Display = new TH2F("h3Display","Charge",32,-0.5,31.5,16,-0.5,15.5);
	h3Display->SetStats(0);
	h3Display->Draw("colz");
	DrawMUSICBoundaries();

	eDisplay->cd(2);
	TH2F *h4Display = new TH2F("h4Display","Signal after extraction",32,-0.5,31.5,16,-0.5,15.5);
	h4Display->SetStats(0);
	h4Display->Draw("colz");
	DrawMUSICBoundaries();

	TRandom3 rand;
	int FirstPixelID, FiredPixelID;
	double ChargeThreshold = 0.0;

	while(1)
	{
		int n = rand.Integer(TriggeredEventsID.size());
		tSimulatedEvents->GetEntry(TriggeredEventsID[n]);
		T0->GetEntry(TriggeredEventsID[n]);
		cout<<"Event "<<TriggeredEventsID[n]<<" is triggered"<<endl;

		// This section fills up and plots a 2D histogram
		// using the simulated number of PEs that hit each pixel.
  		eDisplay->cd(1);
  		h3Display->Clear();
  		for(int g=0;g<iNumPixels;g++)
  		{
  			int nx, ny;
  			FindBin(g,&nx,&ny);
  			h3Display->SetBinContent(nx,ny,iPEInPixel->at(g));
  		}
		eDisplay->cd(1)->Modified();
		eDisplay->cd(1)->Update();

		cout<<"First triggered Music is: "<<vTriggerCluster->at(0)<<endl;
		cout<<"Second triggered Music is: "<<vTriggerCluster->at(1)<<endl;

		eDisplay->cd(2);
		h4Display->Reset();

		// This section calcualtes the total charge of pixels on the
		// two triggered music chips from bi-focal coincidence.
		FiredPixelID = -1;
		FirstPixelID = vTriggerCluster->at(0)*8;
		for(int i=FirstPixelID; i<FirstPixelID+16; i++)
		{
			for(int j=SignalStart; j<SignalEnd; j++)
				{
					PixelCharge[i] += ((iFADCTraceInPixel[i])->at(j));
				}
			PixelCharge[i] = PixelCharge[i] - (SignalWidth*Baseline[i]);
		}

		// This section scans through 8 pairs of pixels on the two triggered music
		// chips and finds the possible candidates for the bi-focal pairs by applying
		// a few selection rules. The rules are: 1) Peaks from two pixels happen in 50 ns
		// of each other, and 2) Both peaks are above a certain threshold. At the end, it stores
		// the pixel ID of both pixels from possible candidates in a vector.
		vector<int> vPixCandidate;
		vector<int> vPeakTimeIndex;
		for(int i=FirstPixelID; i<FirstPixelID+8; i++)
		{
			int M1PeakTimeIndex =0;
			int M2PeakTimeIndex =0;
			float M1PeakValue = 0;
			float M2PeakValue = 0;
			for(int j=SignalStart; j<SignalEnd; j++)
			{
				if((iFADCTraceInPixel[i]->at(j) - Baseline[i]) > M1PeakValue)
				{
					M1PeakTimeIndex = j;
					M1PeakValue = iFADCTraceInPixel[i]->at(j) - Baseline[i];
				}

				if((iFADCTraceInPixel[i+8]->at(j) - Baseline[i+8]) > M2PeakValue)
				{
					M2PeakTimeIndex = j;
					M2PeakValue = iFADCTraceInPixel[i+8]->at(j) - Baseline[i+8];
				}
			}

			if((abs(M2PeakTimeIndex - M1PeakTimeIndex) < CoincWindow) && (M1PeakValue > PEThreshold) && (M2PeakValue > PEThreshold))
			{
				vPixCandidate.push_back(i);
				vPixCandidate.push_back(i+8);
				vPeakTimeIndex.push_back(M1PeakTimeIndex);
				vPeakTimeIndex.push_back(M2PeakTimeIndex);
			}
		}

		cout<<"Number of candidates: "<<vPixCandidate.size()<<endl;

		// At this point, we have a list of pair candidates. We use the extarcted charge 
		// of those pair pixels to find out which one is the main fired pixel. If fired pixel
		// is in the 2nd music, we still record its pair in the first music as main fired pixel,
		// because later, vFiredPixels expects the ID of the pixel in the 1st music to give us the neighbors.
		int PixMaxCharge = 0;
		int PeakTimeIndex = 0;

		if(vPixCandidate.size() > 0)
		{
			for (int i=0; i<vPixCandidate.size(); i++)
			{
				cout<<"candidate number "<<i+1<<" is pixel#:"<<vPixCandidate[i]<<endl;
				if(PixelCharge[vPixCandidate[i]] > PixMaxCharge)
				{
					PeakTimeIndex = vPeakTimeIndex[i];
					FiredPixelID = vPixCandidate[i];
					PixMaxCharge = PixelCharge[vPixCandidate[i]];
				}
			}
			if (FiredPixelID > FirstPixelID + 7)
			{
				cout<<"Fired pixel with maximum charge is in the 2nd music chip."<<endl;
				FiredPixelID = FiredPixelID - 8;
			}
		}
		
		cout<<"Chosing Pixel ID: "<<FiredPixelID<<" as Fired Pixel"<<endl;

		cout<<"Peak time index: "<<PeakTimeIndex<<endl;

		cout<<"Total charge of fired pixel: "<<PixMaxCharge<<endl;

		// At this point, we have determined the ID of the main fired pixel.
		// So, we use previously calculated neighbor IDs (stored in vFiredPixels)
		// to fill up a histogram and plot the extracted signal.
		if (FiredPixelID > -1)
		{
			for (int m=0; m<vFiredPixels[FiredPixelID].size(); m++)
			{
  				for(int j=SignalStart; j<SignalEnd; j++)
  					{
  						PixelCharge[m] += ((iFADCTraceInPixel[vFiredPixels[FiredPixelID][m]])->at(j));
  					}
				PixelCharge[m] = PixelCharge[m] - (SignalWidth*Baseline[vFiredPixels[FiredPixelID][m]]);

  				int nx, ny;
  				FindBin(vFiredPixels[FiredPixelID][m], &nx, &ny);
				if(PixelCharge[m] > ChargeThreshold)
				{
	  				h4Display->SetBinContent(nx,ny,PixelCharge[m]*DCPEFactor);
				}
				else
				{
					h4Display->SetBinContent(nx,ny,-1);
				}
			}				
		}
		else
		{
			cout<<"This event is skipped because it did not meet bi-focal requirements."<<endl;
		}

  		// Resetting pixel charge array for next event
  		for(int g=0; g<iNumPixels; g++)
  		{
			PixelCharge[g] = 0;
		}

		h4Display->SetMinimum(0);
		eDisplay->cd(2)->Modified();
		eDisplay->cd(2)->Update();

		if(!HandleInput())
			break;
  	}
}
