#include <iostream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TGraph.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TTimer.h"
#include "TMath.h"

#include "TMatrixD.h"
#include "Math/Vector3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/Rotation3D.h"
#include "TVector3.h"

#include "GPilot.h"
#include "GRootWriter.h"

TRandom3 TR3;

using namespace std;

/*! structure to hold command line entries.
 */
struct Cline {
  string pilotfile;    //!< name of pilot file
  string outFileName;  //!< name of output file
};

/*! structure to hold pilot file entries.
 */
struct Pilot {
  string pilotfile;     //!< name of pilot file
  string outFileName;  //!< name of output file
  string outFileHeaderTree; //!< name of header tree
  string outFileTelTreeName; //!< base name of telescope tree

  int nShower;  //!< number of showers to process (all if < 0)
  int nPhoton;  //!< number of photons to process (all if < 0)

  UInt_t seed;  //!< random number generator seed
  
  unsigned iNInitEvents; //!< vector capacity, see pilot file in Config directory
  bool debugBranchesFlag; //!< if true, create debug branches in output root file
  bool outFileDCos;     //!< add DCos branches to outFileTelTree 
};

// ostream for all logging info. set to cerr in main
ostream *oLog; //!< stream for log file, same in all files.

int commandLineHelp() {

  *oLog << "SimulatePointsourceInFP:  commandline options with default values" << endl;
  *oLog << "    -of <outputFileName> " << endl;
  *oLog << "    -p  <pilotFileName = ./Config/opticsSimulation.pilot>  ";
  return 1;
};
/********************** end of commandLineHelp **************************************/


int readCommandLine(int argc, char *argv[], Cline *cline) {

  // move arguments into strings, easier to work with strings
  vector<string> clArg;
  for (int i = 0;i<argc;i++) {
    clArg.push_back(argv[i]);
  }

  for (int i = 1;i<argc;i++) {
    *oLog << "COMMANDLINE " << clArg.at(i) << endl;
    if (clArg.at(i)=="-p" ) {
      cline->pilotfile = clArg.at(i+1);
      i++;
    }
    else if (clArg.at(i)=="-of" ) {
      cline->outFileName = clArg.at(i+1);
      i++;
    }
    else if (clArg.at(i)=="-h" ) {
      commandLineHelp();
      exit(0);
    }
    else {
      cout << "   INCORRECT COMMMANDLINE SWITCHES" << endl;
      commandLineHelp();
      exit(0);
    }
  }
  return 1;
};
/********************** end of readCommandLine ************************/

int readPilot(Pilot *pilot) {

  pilot->outFileName = "";
  pilot->outFileHeaderTree = "";
  pilot->outFileTelTreeName = "";
  pilot->nShower = -1;
  pilot->nPhoton = -1;
  pilot->iNInitEvents = 100;
  pilot->debugBranchesFlag = false;
  pilot->outFileDCos = false;

  pilot->seed = 0;

  vector<string> tokens;
  string spilotfile = pilot->pilotfile;

  GPilot *pi = new GPilot(spilotfile);

  string flag = "FILEOUT";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->outFileName = tokens.at(0);
    pilot->outFileHeaderTree = tokens.at(1);
    pilot->outFileTelTreeName = tokens.at(2);
    int fileDCos = atoi(tokens.at(3).c_str() );
    if (fileDCos > 0) pilot->outFileDCos = true;
  }
  flag = "NSHOWER";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->nShower = atoi(tokens.at(0).c_str());
    if (tokens.size() == 2) {
      pilot->nPhoton = atoi(tokens.at(1).c_str());
    }
  }  
  flag = "DEBUGBRANCHES";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    int tmp = (Int_t)atoi(tokens.at(0).c_str());
    if (tmp > 0) pilot->debugBranchesFlag = true;
  }  
  flag = "VECCAPACITY";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->iNInitEvents = (UInt_t)atoi(tokens.at(0).c_str());
  }
  flag = "SEED";
  pi->set_flag(flag);
  while (pi->get_line_vector(tokens) >=0) {
    pilot->seed = (UInt_t)atoi(tokens.at(0).c_str());
  }
    
  delete pi;
  return 1;
};
/******************************* end of readPilot *****************************/

int updatePilot(const Cline &cline,Pilot *pilot) {

  if (cline.outFileName != "") {
    pilot->outFileName = cline.outFileName;
  }
  return 1;
};
/****************************** end of updatePilot ******************************/


void WriteAllTelTree(GRootWriter *RootWriter, string treeH)
{
  // cd to writer root file
  RootWriter->cdToWriteRootFile();

  // use a string instead of a character array for the
  // file header.
  string sFileHeader = "Empty";
  string *strP = &sFileHeader;

  vector<int> telTypeVector;
  vector<double> telPlateSFactorVector;
  
  vector<int> telIDVector;
  vector<float> transitTimeVector;
  vector<float> telLocXGCVector;
  vector<float> telLocYGCVector;
  vector<float> telLocZGCVector;

  double telXLocGC = 0.0;
  double telYLocGC = 0.0;
  double telZLocGC = 0.0;

  string sVersion = "ToySim for SPB2 V0";
  double fObsHgt = 0;      //!< observatory height from input record
  double fGlobalEffic = 1;   //!< global efficiency from input record


  int telid = 1;
  double avgTransitTime = 0.0;
  telPlateSFactorVector.push_back(1.0);
  telTypeVector.push_back(0);
  telIDVector.push_back(telid);
  transitTimeVector.push_back( (float)avgTransitTime );

  telLocXGCVector.push_back( (float)telXLocGC);
  telLocYGCVector.push_back( (float)telYLocGC);
  telLocZGCVector.push_back( (float)telZLocGC);

  // create allTel tree, make branches, fill, and write the tree.
  TTree *allTel = new TTree(treeH.c_str(),treeH.c_str());
 
  allTel->Branch("fileHeader",&strP);
  allTel->Branch("GrOpticsVersion",&sVersion);
  allTel->Branch("globalEffic",&fGlobalEffic,"globalEffic/D");
  allTel->Branch("obsHgt",&fObsHgt,"obsHgt/D");
  allTel->Branch("telIDVector",&telIDVector);
  allTel->Branch("telTypeVector", &telTypeVector);
  allTel->Branch("telPlateSFacVector",&telPlateSFactorVector);
  allTel->Branch("telLocXVector",&telLocXGCVector);
  allTel->Branch("telLocYVector",&telLocYGCVector);
  allTel->Branch("telLocZVector",&telLocZGCVector);
  allTel->Branch("transitTimeVector",&transitTimeVector);

  // fill and write tree:  what about transit time map?
  // I added ProcessLine("#include map") in grOptics.cpp
  allTel->Fill();
  allTel->Write();

  //delete[] cstr;
}
/************** end of fillAllTelTree ******************/

int main( int argc, char **argv )
{
 cout<<"This program simulates a point source in the focal plane of the SPB2 camera"<<endl;

 ////////////////////////////////////////////////////////////
 /////// read command line and pilot file 
 /////// set random number generator seed
 /////// open and write to log file
 ////////////////////////////////////////////////////////////
  
  Pilot pilot;
  Cline cline;

  //set command line defaults
  cline.pilotfile = "Simulation.pilot"; // default value
  cline.outFileName   = "SPB2PointSourceSims.root";
  
  // command line printed to cerr in this function
  readCommandLine(argc,argv,&cline);

  // set pilotfile name, read pilot and update from command line
  // structure
  pilot.pilotfile = cline.pilotfile;
  readPilot(&pilot);
  updatePilot(cline,&pilot);

  // set seed so can print out seed if from machine clock
  TR3.SetSeed(pilot.seed); 


 string outFileName = pilot.outFileName;  //!< name of output file
 TFile *fO;
 fO = new TFile(outFileName.c_str(),"RECREATE");
 if (fO->IsZombie() ) {
      *oLog << "error opening root output file: " << outFileName << endl;
      *oLog << "...exiting" << endl;
      exit(-1);
 }

 ///////////////// set up the root writer //////////
 GRootWriter *RootWriter = new GRootWriter(fO,1,pilot.outFileTelTreeName,
					   pilot.outFileDCos, pilot.iNInitEvents,
					   pilot.debugBranchesFlag);


//Filling the tree with the information about all the telescopes in the output
//file. Only one in this case

 WriteAllTelTree(RootWriter,pilot.outFileHeaderTree);

 double fPhotonToCameraTime = 10;
 double fPhotWaveLgt = 300;

 // get ray tracing results
 ROOT::Math::XYZVector vPhotonCameraLoc;
 ROOT::Math::XYZVector vPhotonCameraDcos;

 //configuration parameters of simulation run
 int nEventsPerPELevel=1000;
 float fBifocalOffset = 12.5; //mm
 float xCamPosMax = 105.0; //measured from camera center, will be used +/-
 float yCamPosMax = 55.0; //measured from camera center, will be used +/-

 unsigned int fEventNumber = 0;

for(float logpe=1;logpe<=3;logpe+=0.2)
 {
   float fMeanPhotoElectrons = pow(10,logpe)/2.; //the signal from one shower is split in two

   for(int n=0;n<nEventsPerPELevel;n++)
    {
      //float xCamPos = TR3.Uniform(-0.5*xCamPosMax-fBifocalOffset,xCamPosMax);
      float xCamPos = TR3.Uniform(-xCamPosMax,xCamPosMax);
      float yCamPos = TR3.Uniform(-yCamPosMax,yCamPosMax);
      for(int p=0;p<2;p++) //two points
        {
          int iNumPhotons = TR3.Poisson(fMeanPhotoElectrons/0.127214); //division is by the PDE at 300nm, which is assumed all photons have here. This converts back to photons. Converting back to PEs is done in CARE. So make sure the PDE is the same.
          if(p==1)//yeah this is how to simulate the bifocaloptics
            xCamPos+= fBifocalOffset;
          for(int i=0;i<iNumPhotons;i++)
          { 
            vPhotonCameraLoc.SetX(xCamPos);   //units are mm
            vPhotonCameraLoc.SetY(yCamPos);  //[mm] note that the Y coordinates are multiplied with -1 before they are written to the file.

            // mRootWriter is passed in as a pointer so have to dereference
            // add photon to the appropriate writer if photon strikes the camera 
            RootWriter->addPhoton(vPhotonCameraLoc,vPhotonCameraDcos,
	   				fPhotonToCameraTime,
  					fPhotWaveLgt);
          }
       	} 
	   
     // add event
     fEventNumber++;
     unsigned int fPrimaryType = 1; //gamma ray. 
     double       fPrimaryEnergy = fMeanPhotoElectrons*2;


      // primary details
      ROOT::Math::XYZVector vSCore; //!< core loc.vec. ground coors.(meters)
      vSCore.SetX(xCamPos-fBifocalOffset);
      vSCore.SetY(yCamPos);
      ROOT::Math::XYZVector vSDcosGd; //!< core dir.cosines; grd coors.
      double fZnPrim = 0;                   //!< primary zenith angle
      double fAzPrim = 0;                   //!< primary azimuthal angle

      double fWobbleTE = 0; //!< wobble offset E, determined from fWobbleR
      double fWobbleTN = 0; //!< wobble offset N, determined from fWobbleR

      double fFirstIntHgt = 0;
      double fFirstIntDpt = 0;
      unsigned int iShowerID = 0;


      //unsigned int tel = RootWriter->getTelID();
      double fDelay = 0.0;
      double azTel = 0.0;
      double znTel = 0.0;
      double srcX = 0.0;
      double srcY = 0.0;
      // get telescope zn/az and src relative to telescope for writer
      //(*mArrayTel)[tel]->getAzZnTelescope(&azTel,&znTel);
      //(*mArrayTel)[tel]->getSrcRelativeToCamera(&srcX,&srcY);
      // get core locations for debug branches
      ROOT::Math::XYZVector vSCoreTC;
      ROOT::Math::XYZVector vSDcosTC;
      ROOT::Math::XYZVector vSCoreSC;  //!< primary coreHit primary(shower) coor.
      ROOT::Math::XYZVector vSDcosSC;  //!< primary dirCos. primary(shower) coor.
      ROOT::Math::XYZVector vTelLocTC;
      //(*mArrayTel)[tel]->getCoreLocDCosTC(&vSCoreTC,&vSDcosTC);
      //(*mArrayTel)[tel]->getCoreLocDCosSC(&vSCoreSC,&vSDcosSC);
      //(*mArrayTel)[tel]->getTelLocTC(&vTelLocTC);
      RootWriter->addEvent(fEventNumber,
					 fPrimaryType,
					 fPrimaryEnergy,
					 vSCore,
					 vSDcosGd,
					 fWobbleTE,
					 fWobbleTN,
					 fDelay,
					 fPhotonToCameraTime,
                                         azTel,znTel,
                                         fAzPrim,fZnPrim,
                                         srcX,srcY,
                                         fFirstIntHgt,fFirstIntDpt,iShowerID,
                                         vSCoreTC,vSDcosTC,
                                         vSCoreSC,vSDcosSC,
                                         vTelLocTC
                                         );
     }
 }
   RootWriter->write(); 
   fO->Close();
}
/********************** end of main ***************************************/
