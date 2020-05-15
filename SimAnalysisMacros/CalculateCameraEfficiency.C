void CalculateCameraEfficiency(string fInputFileName = "../data/TriggerEfficiency.root")
{
       //Open the CARE file
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
       UInt_t eventNumber ; 
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

      //Stuff for visualizing
      TCanvas *cDisplay = new TCanvas("cDisplay","Display",2000,500);
      cDisplay->Divide(2,1);
      cDisplay->cd(1);
      TH1F *hEfficiency = new TH1F("hEfficiency","Trigger Efficiency",1002,9.9,1001);
      hEfficiency->GetXaxis()->SetTitle("Total Shower Photoelectron Signal");
      hEfficiency->GetYaxis()->SetTitle("Trigger Efficiency");
      hEfficiency->SetStats(0);
      hEfficiency->SetLineWidth(2);
      hEfficiency->Sumw2();
      
      TH2F *hTriggerInEfficiencyAcrossCamera = new TH2F("hTriggerEfficiencyAcrossCamera","Trigger Inefficiency across the camera",214,-107,107,112,-56,56);
      //TH2F *hTriggerInEfficiencyAcrossCamera = new TH2F("hTriggerEfficiencyAcrossCamera","Trigger Inefficiency across the camera",234,-122,112,130,-60,70);
      hTriggerInEfficiencyAcrossCamera->SetStats(0);


       //looping over events doing something
       for(int n=0;n<tSimulatedEvents->GetEntries();n++)
         {
           tSimulatedEvents->GetEntry( n );
           if(arrayTriggerBit)
             {
               hEfficiency->Fill(energy);
             }
           if(arrayTriggerBit==0 && energy>100)
             hTriggerInEfficiencyAcrossCamera->Fill(xcore,ycore);
         }
       hEfficiency->Scale(1e-3*(4*105*55)/(199.8*99.8)); //divide by the number of simulated events per PE bin, multiply with the camera area of which the events have bin simulated and divide by the actual area of the camera
       hEfficiency->Draw();
       gPad->SetLogx();
       cDisplay->cd(2);
       hTriggerInEfficiencyAcrossCamera->Draw("colz");
       TBox *b = new TBox(-96.8-3.1,-46.8-3.1,96.8+3.1,46.8+3.1);
       b->SetFillStyle(0);
       b->SetLineColor(kRed);
       b->Draw(); 
}
