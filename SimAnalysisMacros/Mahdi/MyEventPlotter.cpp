#include <iostream>
#include <iterator>
#include <TH1F.h>
#include <TTimer.h>
#include <TGraph.h>
#include <TMath.h>
#include <TProfile.h>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

void FindNeighborPixels();
void PlotSPB2Events(string fInputFileName);

int iLastPix = -1;
const int iNumPixels = 512; 
const int SignalStart = 4;
const int SignalEnd = 10;
const int SignalWidth = SignalEnd - SignalStart;
const int CoincWindow = 5;
double Baseline[iNumPixels];
double PixelCharge[iNumPixels];
vector< vector<Int_t> *>  iFADCTraceInPixel;
vector<Int_t>  *iPEInPixel;
vector<vector<int> > vFiredPixels;
vector<vector<double> > BaselineDist;
TLatex *text = 0;
TCanvas *cDisplay = 0;
TH1F *hPixelTrace =0;
TH1F *hBaseDist =0;
TH1F *hBaseTemp =0;
TH1F *hBaseMean =0;
TH1F *hBaseRMS =0;

int main(){

	FindNeighborPixels();

	string FileName = "/home/mahdi/Programs/SPB2/SPB2_CT_Simulations/data/test_50.root";

	//auto start = high_resolution_clock::now();
	PlotSPB2Events(FileName);
	//auto stop = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(stop - start);
	//cout << duration.count() << " microseconds" << endl;

}

Bool_t HandleInput()
{
  TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
  while (1)
    {
      //
      // While reading the input process gui events asynchronously
      //
      timer.TurnOn();
      TString input = Getline("Type 'q' to exit, <return> to go on:");
      timer.TurnOff();

      if (input=="q\n")
        return kFALSE;

      if (input=="\n")
        return kTRUE;
    };

  return kFALSE;
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

	cDisplay->cd(3);
	hPixelTrace->Draw();
	gPad->Modified();
	gPad->Update();
}

void ShowBaseDist(int iPix)
{
  if(hBaseDist==0)
  {
    hBaseDist = new TH1F("hBaseDist","Pixel Baseline Dist.",50,495.0,505.0);
    hBaseDist->SetStats(0);
    hBaseDist->GetXaxis()->SetTitle("Baseline");
    hBaseDist->GetYaxis()->SetTitle("Baseline Frequency");
  }

  hBaseDist->Reset();
  TString title2;
  title2.Form("Baseline Dist. of Pixel %i",iPix);
  hBaseDist->SetTitle(title2);
  for(int s=0;s<BaselineDist.size();s++)
  {
    hBaseDist->Fill(BaselineDist[s][iPix]);
  }
  
  cDisplay->cd(7);
  hBaseDist->Draw();
  gPad->Modified();
  gPad->Update();
}

void GetPixMeanRMS()
{
	if(hBaseTemp==0)
	{
		hBaseTemp = new TH1F("hBaseTemp","Temporary Histogram for Baseline",100,495.0,505.0);
	}

	if(hBaseMean==0)
	{
		hBaseMean = new TH1F("hBaseMean","Pixel Baseline Mean",50,499.0,500.0);
		hBaseMean->SetStats(0);
		hBaseMean->GetXaxis()->SetTitle("Mean");
		hBaseMean->GetYaxis()->SetTitle("Mean Frequency");
	}

	if(hBaseRMS==0)
	{
		hBaseRMS = new TH1F("hBaseRMS","Pixel Baseline RMS",50,1.0,1.5);
		hBaseRMS->SetStats(0);
		hBaseRMS->GetXaxis()->SetTitle("RMS");
		hBaseRMS->GetYaxis()->SetTitle("RMS Frequency");
	}

	hBaseTemp->Reset();
	hBaseMean->Reset();
	hBaseRMS->Reset();

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
	
	cDisplay->cd(5);
	hBaseMean->Draw();
	gPad->Modified();
	gPad->Update();

	cDisplay->cd(6);
	hBaseRMS->Draw();
	gPad->Modified();
	gPad->Update();
}

int FindPixel(int x, int y)
{
	//cout<<x<<"  "<<y<<endl;
	//column
	int MUSIC_column = x/2;

	//cout<<"MUSIC column "<<MUSIC_column<<endl;
	//row
	int MUSIC_row = y/4;
	//cout<<"MUSIC row "<<MUSIC_row<<endl;

	//MUSIC ID
	int MUSIC_ID = MUSIC_column+MUSIC_row*16;
	//cout<<"MUSIC_ID "<<MUSIC_ID<<endl;

	//Channel in MUSIC
	int MUSIC_Channel = y%4+4*(x%2);
	//cout<<"MUSIC_Channel "<<MUSIC_Channel<<endl;

	return MUSIC_row*8*16+MUSIC_column*8+MUSIC_Channel;
}

void ShowInfoAtCursor(int x, int y)
{
	//cout<<x<<"  "<<y<<endl;
	//column
	int MUSIC_column = x/2;

	//cout<<"MUSIC column "<<MUSIC_column<<endl;
	//row
	int MUSIC_row = y/4;
	//cout<<"MUSIC row "<<MUSIC_row<<endl;

	//MUSIC ID
	int MUSIC_ID = MUSIC_column+MUSIC_row*16;
	//cout<<"MUSIC_ID "<<MUSIC_ID<<endl;

	//Channel in MUSIC
	int MUSIC_Channel = y%4+4*(x%2);
	//cout<<"MUSIC_Channel "<<MUSIC_Channel<<endl;

	int PixID = MUSIC_row*8*16+MUSIC_column*8+MUSIC_Channel;

	//cout<<"MUSIC_ID: "<<MUSIC_ID<<" MUSIC_Channel: "<<MUSIC_Channel<<" Pixel ID: "<<PixID<<endl;
	TString statusline;
	statusline.Form("MUSIC_ID: %i MUSIC_Channel: %i Pixel ID: %i  PEs in Pixel: %i",MUSIC_ID,MUSIC_Channel,PixID,iPEInPixel->at(PixID));
	if(text!=0)
		text->Delete();
	TLatex T1;
	text = T1.DrawLatexNDC(0.1,0.95,statusline.Data());
	gPad->Modified();
	gPad->Update();
}

int FindPix(int nx, int ny)
{
	int PixelID;
	PixelID = ((nx-1)*4+((ny-1)%4)+((ny-1)>>2)*128);
	return PixelID;
}

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
	//cout<<"Pixel "<<iPix<<" x: "<<*nx<<" y: "<<*ny<<endl;
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
					{
						vNeighbor.push_back(tmp);
					}
				}
			}
		}
		vFiredPixels.push_back(vNeighbor);			
	}
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
	//cout<<px<<"  "<<py<<endl;
	Float_t xx = gPad->AbsPixeltoX(px);
	Float_t yy = gPad->AbsPixeltoY(py);
	Float_t x = 0.5+gPad->PadtoX(xx);
	Float_t y = 0.5+gPad->PadtoY(yy);
	//cout<<x<<"  "<<y<<endl;
	int pix = FindPixel((int)x,(int)y);
	if(pix!=iLastPix)
	{
	  iLastPix=pix;
	  ShowInfoAtCursor((int) x, (int) y);
	}

	if (event == 11){
	  //cout<<pix<<endl;  
	  ShowPixeltrace(pix);
	}
}

void PixelClicked2()
{
	int event = gPad->GetEvent();
	TObject *o = gPad->GetSelected();
	if (!o) return;
	if (!(o->InheritsFrom("TH2")))
	   return;
	TH2F *h = (TH2F*)o;
	int px = gPad->GetEventX();
	int py = gPad->GetEventY();
	//cout<<px<<"  "<<py<<endl;
	Float_t xx = gPad->AbsPixeltoX(px);
	Float_t yy = gPad->AbsPixeltoY(py);
	Float_t x = 0.5+gPad->PadtoX(xx);
	Float_t y = 0.5+gPad->PadtoY(yy);
	//cout<<x<<"  "<<y<<endl;
	int pix = FindPixel((int)x,(int)y);
	if(pix!=iLastPix)
	{
	  iLastPix=pix;
	  //ShowInfoAtCursor((int) x, (int) y);
	}

	if (event == 11){
	  //cout<<pix<<endl;
	  ShowBaseDist(pix);
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

void PlotSPB2Events(string fInputFileName)
{
	TFile *fO = new TFile( fInputFileName.c_str(), "READ" );
	if( fO->IsZombie() )
	{
		cout << "error opening root input file: " << fInputFileName << endl;
		cout << "...exiting" << endl;
		exit( -1 );
	}

	cout<<"Have opened the file with the simulated events: "<<fInputFileName.c_str()<<endl;

    UInt_t uNumTelescopes;
    Float_t energy;
    UInt_t eventNumber;
    Float_t xcore ;
    Float_t ycore ;
    Float_t azPrim ;
    Float_t znPrim ;
    Bool_t arrayTriggerBit ;
    std::vector< Bool_t > *vTelescopeTriggerBits =0 ;
    Float_t DeltaTL3 ;

    // read tree with all the event by event information
    cout<<"Looking for Tree tSimulatedEvents"<<endl;
    TTree *tSimulatedEvents = (TTree*)fO->Get( "Events/tSimulatedEvents" );
	if( !tSimulatedEvents )
	{
		cout << "error: tree tSimulatedEvents not found in " << fInputFileName << endl;
	    cout << "...exiting" << endl;
	    exit( -1 );
	}

    tSimulatedEvents->SetBranchAddress("energy",&energy);
    tSimulatedEvents->SetBranchAddress("ZnPrim",&znPrim);
    tSimulatedEvents->SetBranchAddress("AzPrim",&azPrim);
    tSimulatedEvents->SetBranchAddress("xcore",&xcore);
    tSimulatedEvents->SetBranchAddress("ycore",&ycore);
    tSimulatedEvents->SetBranchAddress("arrayTriggerBit",&arrayTriggerBit);
    tSimulatedEvents->SetBranchAddress("uNumTelescopes",&uNumTelescopes);
    tSimulatedEvents->SetBranchAddress("eventNumber",&eventNumber);
    tSimulatedEvents->SetBranchAddress("DeltaTL3",&DeltaTL3);
    tSimulatedEvents->SetBranchAddress("vTelescopeTriggerBits",&vTelescopeTriggerBits);

    cout<<"Looking for tree with camera output"<<endl;
    TTree *T0 = (TTree*)fO->Get( "Events/T0" );
    if( !T0 )
    {
    	cout << "error: tree T0 not found in " << fInputFileName << endl;
    	cout << "...exiting" << endl;
    	exit( -1 );
    }

    vector<int> *vTriggerCluster = 0;
    vector<Float_t> *fTimeOverThreshold;
    vector<Float_t>   *fSumTimeInPixel;
    vector<Int_t>   *iQDCInPixel;
    Int_t           iNumPhotonsInFocalPlane;
    Float_t         fAzTel ;
    Float_t         fZnTel ;
    vector<Bool_t>  *bInLoGain;
    iFADCTraceInPixel.assign(iNumPixels,0);

    T0->SetBranchAddress("vGroupsInTriggerCluster",&vTriggerCluster);
    //T0->SetBranchAddress("vTimeOverThreshold",&fTimeOverThreshold);
    //T0->SetBranchAddress("vSumTimeInPixel", &fSumTimeInPixel);
    T0->SetBranchAddress("vPEInPixel", &iPEInPixel);
    //T0->SetBranchAddress("vQDCValue", &iQDCInPixel);
    T0->SetBranchAddress("iPhotonsInFocalPlane", &iNumPhotonsInFocalPlane);
    T0->SetBranchAddress("fAzTel", &fAzTel);
    T0->SetBranchAddress("fZnTel", &fZnTel);
    //T0->SetBranchAddress("vHiLoGainBit", &bInLoGain);
    TString name;
    for(int g=0;g<iNumPixels;g++)
    {
    	name.Form("vFADCTraces%i",g);
    	T0->SetBranchAddress(name,&(iFADCTraceInPixel[g]));
    }

    int BaselineWidth = 4;
	float MultiplyFactor1 = 0;
	float MultiplyFactor2 = (1.0/BaselineWidth);
	int TotalEvents = tSimulatedEvents->GetEntries();
	vector<int> TriggeredEventsID;

	cout<<"Total number of events: "<<TotalEvents<<endl;

	for (int i=0; i<TotalEvents; i++)
	{
		tSimulatedEvents->GetEntry(i);
		if(arrayTriggerBit)
		{
			T0->GetEntry(i);
			vector<double> tmp;
			for(int j=0;j<iNumPixels;j++)
			{
				double tmpbase;
				for (int k=0; k<BaselineWidth; k++)
				{
					tmpbase += ((iFADCTraceInPixel[j])->at(k));
				}
				tmp.push_back(tmpbase*MultiplyFactor2);
				Baseline[j] += tmpbase;
				tmpbase = 0;
			}
			BaselineDist.push_back(tmp);
			TriggeredEventsID.push_back(i);
		}
	}

	cout<<"Total triggered events: "<<TriggeredEventsID.size()<<endl;
	MultiplyFactor1 = (1.0/(BaselineWidth*TriggeredEventsID.size()));

	for(int i=0;i<iNumPixels;i++)
	{
		Baseline[i] = Baseline[i]*MultiplyFactor1;
	}

	cDisplay = new TCanvas("cDisplay","Display",2000,1000);
	cDisplay->Divide(3,3);

	GetPixMeanRMS();

	cDisplay->cd(4);
	gPad->AddExec("ex","PixelClicked2()");
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

	cDisplay->cd(4)->Modified();
	cDisplay->cd(4)->Update();

	cDisplay->cd(2);
	TH2F *h4Display = new TH2F("h4Display","Signal after extraction",32,-0.5,31.5,16,-0.5,15.5);
	h4Display->SetStats(0);
	h4Display->Draw("colz");
	DrawMUSICBoundaries();

	cDisplay->cd(1);
	gPad->AddExec("ex","PixelClicked()");
	TH2F *h3Display = new TH2F("h3Display","Charge",32,-0.5,31.5,16,-0.5,15.5);
	h3Display->SetStats(0);
	h3Display->Draw("colz");
	DrawMUSICBoundaries();

	// This section calculates the DC and PE for each event and its fired pixels
	vector<int> DCValue;
	vector<int> PEValue;
	int BeginPixelID =0;
	int TotalCharge =0;

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
	cDisplay->cd(8);
	float MultiplyFactor3;
	TH2D *hScatDCPE = new TH2D("hScatDCPE","DC vs. PE Histogram",100,-20,300,100,-20,600);
	hScatDCPE->Clear();
	hScatDCPE->SetStats(0);
	hScatDCPE->GetXaxis()->SetTitle("PE");
	hScatDCPE->GetYaxis()->SetTitle("DC");
	for (int i=0; i<PEValue.size(); i++)
	{
	  hScatDCPE->Fill(PEValue[i], DCValue[i]);
	}
	TF1 *fDCPE = new TF1("fDCPE","[0]*x",0,300);
	fDCPE->SetParNames("slope");
	hScatDCPE->Fit(fDCPE);
	hScatDCPE->Draw("SCAT");
	MultiplyFactor3 = (1.0/fDCPE->GetParameter(0));
	cDisplay->cd(8)->Modified();
	cDisplay->cd(8)->Update();

	// This section calculates the standard deviation of mean for the DC/PE histogram.
	cDisplay->cd(9);
	TProfile* DCProfile;
	DCProfile = new TProfile("DCProfile","Profile of extracted digital counts",50,0.0,300.0);
	DCProfile->GetXaxis()->SetTitle("PE");
	DCProfile->GetYaxis()->SetTitle("DC");
	DCProfile->SetStats(0);
	for (int i=0; i<PEValue.size(); i++)
	{
	  DCProfile->Fill(PEValue[i], DCValue[i]);
	}
	DCProfile->Draw();
	cDisplay->cd(9)->Modified();
	cDisplay->cd(9)->Update();	

	
	TRandom3 rand;
	int FirstPixelID, FiredPixelID;
	int DCThreshold = 510;
	double ChargeThreshold = 0.0;


	while(1)
	{
		int n = rand.Integer(tSimulatedEvents->GetEntries());
  		tSimulatedEvents->GetEntry(n);
  		if(arrayTriggerBit)
  		{
  			// This section fills up and plots a 2D histogram
  			// using the simulated number of PEs that hit each pixel.
	  		h3Display->Clear();
	  		cDisplay->cd(1);
	  		cout<<"Event "<<n<<" is triggered"<<endl;
	  		T0->GetEntry(n);
	  		for(int g=0;g<iNumPixels;g++)
		  		{
  					int nx, ny;
  					FindBin(g,&nx,&ny);
  					h3Display->SetBinContent(nx,ny,iPEInPixel->at(g));
  				}
  			cDisplay->cd(1)->Modified();
  			cDisplay->cd(1)->Update();

  			h4Display->Reset();
  			cDisplay->cd(2);

  			// This section calcualtes the total charge of pixels on the
  			// two triggered music chips from bi-focal coincidence.
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
  			for(int i=FirstPixelID; i<FirstPixelID+8; i++)
  			{
  				int M1PeakTimeIndex =0;
  				int M2PeakTimeIndex =0;
  				int M1PeakValue = 0;
  				int M2PeakValue = 0;
  				for(int j=SignalStart; j<SignalEnd; j++)
  				{
  					if(iFADCTraceInPixel[i]->at(j) > M1PeakValue)
  					{
  						M1PeakTimeIndex = j;
  						M1PeakValue = iFADCTraceInPixel[i]->at(j);
  					}
  					if(iFADCTraceInPixel[i+8]->at(j) > M2PeakValue)
  					{
  						M2PeakTimeIndex = j;
  						M2PeakValue = iFADCTraceInPixel[i+8]->at(j);
  					}
  				}

				if((abs(M2PeakTimeIndex - M1PeakTimeIndex) < CoincWindow) && (M1PeakValue > DCThreshold) && (M2PeakValue > DCThreshold))
				{
					vPixCandidate.push_back(i);
					vPixCandidate.push_back(i+8);
				}
			}

			// At this point, we have a list of pair candidates. We use the extarcted charge 
			// of those pair pixels to find out which one is the main fired pixel. If fired pixel
			// is in the 2nd music, we still record its pair in the first music as main fired pixel,
			// because later, vFiredPixels expects the ID of the pixel in the 1st music to give us the neighbors.
			int PixMaxCharge = 0;
			for (int i=0; i<vPixCandidate.size(); i++)
			{
				if(PixelCharge[vPixCandidate[i]] > PixMaxCharge)
				{
					FiredPixelID = vPixCandidate[i];
					PixMaxCharge = PixelCharge[vPixCandidate[i]];
				}
			}
			if (FiredPixelID > FirstPixelID + 7)
			{
				cout<<"Fired pixel with maximum charge is in the 2nd music chip."<<endl;
				FiredPixelID = FiredPixelID - 8;
			}
			cout<<"Chosing Pixel ID: "<<FiredPixelID<<" as Fired Pixel"<<endl;

			// At this point, we have determined the ID of the main fired pixel.
			// So, we use previously calculated neighbor IDs (stored in vFiredPixels)
			// to fill up a histogram and plot the extracted signal.
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
	  				h4Display->SetBinContent(nx,ny,PixelCharge[m]*MultiplyFactor3);
				}
				else
				{
					h4Display->SetBinContent(nx,ny,-1);
				}
			}

	  		// Resetting pixel charge array for next event
	  		for(int g=0; g<iNumPixels; g++)
	  		{
				PixelCharge[g] = 0;
			}

  			h4Display->SetMinimum(0);
  			cDisplay->cd(2)->Modified();
  			cDisplay->cd(2)->Update();

  			if(!HandleInput())
  				break;
  		}
  	}
}