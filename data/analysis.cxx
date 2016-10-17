#define Mkaon 0.493677
#define Mproton 0.93827231
#define Melectron 0.00051099907
#define PI 3.1415927

//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "sys/types.h"
#include "dirent.h"

#include "math.h"
#include "string.h"
//Add the data structure
#include <iomanip>
//#include <stdio.h> 

#ifndef __CINT__  
#include "TROOT.h"
#include "TFile.h" 
#include "TChain.h"
#include "TMath.h"
#include "TH1.h" 
#include "TH2.h"   
#include "TH3.h" 
#include "TF1.h" 
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "pevent.h"
#include "anaCuts.h"
#include "runnumber.h"
using std::cout;
using std::endl;
using std::setw;
#endif 

enum pairMode{
	OneEID=1,
	TwoEID=2
};

const int dsmHotTowers[106] = {
	31,51,114,275,293,479,509,533,555,561,639,681,740,743,749,772,779,799,840,860,880,893,897,982,986,993,1063,1142,1160,1200,1224,1232,1237,1241,1256,1263,1280,1284,1306,1313,1318,1337,1438,1486,1537,1592,1709,1713,1732,1823,1850,1856,1879,1945,1976,1984,2043,2145,2162,2164,2190,2202,2290,2299,2300,2313,2339,2414,2439,2459,2529,2580,2633,2652,2834,2863,2865,2874,3005,3020,3061,3137,3154,3420,3515,3532,3692,3720,3737,3838,3861,3925,3948,4013,4130,4169,4262,4316,4430,4458,4534,4560,4595,4684,4766,4781
};

void bookHistograms(char* outFile);
bool passEvent(pevent* event);
bool passTrack(pevent* event, int i);
void makeRealPairs();
void fillHistograms(std::string  unlikeOrlike,TLorentzVector JPSI);
void fill3DHistograms(std::string unlikeOrlike,TLorentzVector JPSI,int i,int j,int pairs);
void writeHistograms(char* outFile);
void printCuts();
void setPairingMode(pairMode mode);
double getHadronWt(double pt, double eta);
void Polarization(int icharge,int jcharge,TLorentzVector ivector,TLorentzVector jvector,TH2F *hist0, TH2F *hist1, TH2F *hist2, TH2F *hist3, bool isSignal);

const int nHPtCuts = 3;

TH1F *checkcut;
TH2F *hNe;
TH1F *hnTrack;
TH1F *hnEvent;

TH2F *hJpsiCosThetaPt;
TH2F *hJpsiCosThetaPtBG;

TH2F *hJpsiPhiPt;
TH2F *hJpsiPhiPtBG;

TH1F *hJpsiMass;
TH1F *hJpsiMassBG;
TH1F *hJpsiEta;
TH1F *hJpsiEtaBG;
TH1F *hJpsiRapidity;
TH1F *hJpsiRapidityBG;
TH1F *hJpsiPt;
TH1F *hJpsiPtBG;

TH2F *hJpsiMassPt;
TH2F *hJpsiMassPtBG;

TH2F *hJpsiCosThetaPtCS;
TH2F *hJpsiCosThetaPtCSBG;

TH2F *hJpsiPhiPtCS;
TH2F *hJpsiPhiPtCSBG;

TH2F *hJpsiPhiPtHX;
TH2F *hJpsiPhiPtHXBG;

TH3F *hJpsiCosThetaInvMPt;
TH3F *hJpsiCosThetaInvMPtBG;
TH3F *hJpsiCosThetaInvMPtCS;
TH3F *hJpsiCosThetaInvMPtCSBG;

TH3F *hJpsiPhiInvMPt;
TH3F *hJpsiPhiInvMPtBG;
TH3F *hJpsiPhiInvMPtCS;
TH3F *hJpsiPhiInvMPtCSBG;

//TH3F *dsmAdcInvMPt;
//TH3F *dsmAdcInvMPtBG;
//TH3F *AdcInvMPt;
//TH3F *AdcInvMPtBG;

TH3F *hJpsiCosThetaPhiPt;
TH3F *hJpsiCosThetaPhiPtBG;
TH3F *hJpsiCosThetaPhiPtCS;
TH3F *hJpsiCosThetaPhiPtCSBG;

TH3F *hNumJpsiPtMTofMult;
TH3F *hDenJpsiPtMTofMult;

TH1F *hnsigmaE;
TH1F *hemcPtMin;
TH1F *hDCA;
TH1F *hnhitsDEDX;

TTree *tree;
Float_t positron_theta_hx=-99.,positron_theta_cs=-99.,positron_phi_hx=-99.,positron_phi_cs=-99.;
Float_t electron_theta_hx=-99.,electron_theta_cs=-99.,electron_phi_hx=-99.,electron_phi_cs=-99.;
Float_t jpsi_pt,jpsi_eta,jpsi_phi,jpsi_InvM;
Float_t lepton1_pt,lepton1_eta,lepton1_phi,lepton1_InvM;
Float_t lepton2_pt,lepton2_eta,lepton2_phi,lepton2_InvM;
TLorentzVector lepton1,lepton2;
bool lepton1_isTPCe,lepton2_isTPCe,lepton1_isTOFe,lepton2_isTOFe,lepton1_isEMCe,lepton2_isEMCe;

pairMode mPairingMode = OneEID;
const int MaxEmcElectron = 1000, MaxTpcElectron = 1000, MaxJpsi = 20, MaxHadron = 1000;
int iran = 0;
int current_nEmce, current_nTpce, current_nJpsi, current_nJpsiLS, current_nHadron, current_nTrigE, current_nJpsiSide;
int current_centrality = 0;
TLorentzVector current_emce[MaxEmcElectron];
TLorentzVector current_tpce[MaxTpcElectron];

int current_jpsiDau1Id[MaxJpsi], current_jpsiDau2Id[MaxJpsi];
TLorentzVector current_jpsi[MaxJpsi];
int current_jpsiLSDau1Id[MaxJpsi], current_jpsiLSDau2Id[MaxJpsi];
TLorentzVector current_jpsiLS[MaxJpsi];
int current_jpsiSideDau1Id[MaxJpsi], current_jpsiSideDau2Id[MaxJpsi];
TLorentzVector current_jpsiSide[MaxJpsi];

int current_hadronId[MaxHadron];
TVector3 current_hadron[MaxHadron];

char current_tpceTofTag[MaxTpcElectron];
int current_emceId[MaxEmcElectron], current_tpceId[MaxTpcElectron];
char current_emceq[MaxEmcElectron], current_tpceq[MaxTpcElectron];
int current_emcdsmadc[MaxEmcElectron], current_tpcdsmadc[MaxTpcElectron];
int current_emcadc0[MaxEmcElectron], current_tpcadc0[MaxTpcElectron];
char current_emcIsTrg[MaxEmcElectron];
bool  current_EmcHits[MaxEmcElectron];
bool current_TpcHits[MaxTpcElectron];

int vzBufferPointer, cenBufferPointer;
const int nVzBin = 100;
const int nCenBin = 1;
const int nMaxEventsInBuffer = 10;
int nEventsInBuffer[nCenBin][nVzBin];
bool buffer_fullFlag[nCenBin][nVzBin];
int buffer_nEmce[nCenBin][nVzBin][nMaxEventsInBuffer], buffer_nTpce[nCenBin][nVzBin][nMaxEventsInBuffer];
TLorentzVector buffer_emce[nCenBin][nVzBin][nMaxEventsInBuffer][MaxEmcElectron];
TLorentzVector buffer_tpce[nCenBin][nVzBin][nMaxEventsInBuffer][MaxTpcElectron];
char buffer_emceq[nCenBin][nVzBin][nMaxEventsInBuffer][MaxEmcElectron];
char buffer_emcIsTrg[nCenBin][nVzBin][nMaxEventsInBuffer][MaxEmcElectron];
char buffer_tpceq[nCenBin][nVzBin][nMaxEventsInBuffer][MaxTpcElectron];
char buffer_tpceTofTag[nCenBin][nVzBin][nMaxEventsInBuffer][MaxTpcElectron];

int buffer_nJpsi[nCenBin][nVzBin][nMaxEventsInBuffer], buffer_nHadron[nCenBin][nVzBin][nMaxEventsInBuffer];
TLorentzVector buffer_jpsi[nCenBin][nVzBin][nMaxEventsInBuffer][MaxJpsi];
TVector3 buffer_hadron[nCenBin][nVzBin][nMaxEventsInBuffer][MaxHadron];
int nHadrons[4] = {0,0,0,0};

bool  isBHT0;
bool  isBHT1;
bool  isBHT2;
int mHTselect = 5;
int mSelectTOFHadrons = 0;
int runPointer;
TRandom3 *myRandom;
TFile *mFile;
double tpcEffPars[20][3];
double tofEffPars[20][3];
double zdcRate = 0.;
int tofMult = 0;
TF1 *fTpcEff;
TF1 *fTofEff;

int mNHitsFit;
int mDsmAdc = 1;
int mPOE = 0;

int main(int argc, char** argv) 
{
	if(argc!=1&&argc!=7) return -1;
	char *inFile="testana.list";
	char outFile[1024];
	sprintf(outFile,"test");
	if(argc==7){
		inFile = argv[1];
		//sprintf(outFile,"%s_%s",mHeader,argv[2]);
		sprintf(outFile,"%s_%s",argv[2],argv[6]);
		pairMode mode = static_cast<pairMode>(atoi(argv[3]));
		cout<<"input arguements: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<endl;
		if(mode>TwoEID||mode<OneEID){ 
			cout<<"wrong input mode : "<<argv[3]<<endl;
			return -1;
		}
		setPairingMode(mode);
		mHTselect = atoi(argv[4]);
		mSelectTOFHadrons = atoi(argv[5]);
		mNHitsFit = atoi(argv[6]);
		if(atoi(argv[6])==29) mDsmAdc = 1.05;
		if(atoi(argv[6])==30) mDsmAdc = 0.95;
		if(atoi(argv[6])==35) mPOE = 1;

		if(mSelectTOFHadrons!=0 && mSelectTOFHadrons!=1 ){ 
			cout<<"wrong input mode mSelectTOFHadrons = "<<mSelectTOFHadrons<<endl;
			return -1;
		}
		cout<<"input HT trigger = "<<mHTselect<<endl;
	}

	myRandom = new TRandom3();

	bookHistograms(outFile);
	hnEvent->Fill(10);

	//+---------------------------------+
	//| open files and add to the chain |
	//+---------------------------------+
	TChain *chain = new TChain("pevent");

	Int_t ifile=0;
	char filename[512];
	ifstream *inputStream = new ifstream;
	inputStream->open(inFile);
	if (!(inputStream)) {
		printf("can not open list file\n");
		return 0;
	}
	for(;inputStream->good();){
		inputStream->getline(filename,512);
		if(inputStream->good()) {
			TFile *ftmp = new TFile(filename);
			if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
				cout<<"something wrong"<<endl;
			} else {
				cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
				chain->Add(filename);
				ifile++;
			}
			delete ftmp;
		}
	}
	delete inputStream;
	pevent *eventtree = new pevent(chain);

	ifstream inf("/global/homes/h/huangbc/pwg/eliza14_bk/jpsi/emb/pikp_run12/txt/comb_eff_fit.txt");
	if(!inf.good()){ 
		cout<<"No input TPC efficiency tables!"<<endl;
		return -1;
	}
	for(int i=0;i<20;i++){
		for(int j=0;j<3;j++) inf>>tpcEffPars[i][j];
	}
	fTpcEff = new TF1("fTpcEff","[0]*exp(-pow([1]/x,[2]))",0.2,15);
	inf.close();

	inf.open("/global/homes/h/huangbc/pwg/eliza14_bk/jpsi/pp200run12/tofeff/txt/comb_tofeff_fit.txt");
	if(!inf.good()){ 
		cout<<"No input TOF efficiency tables!"<<endl;
		return -1;
	}
	for(int i=0;i<20;i++){
		for(int j=0;j<3;j++) inf>>tofEffPars[i][j];
		//for(int j=0;j<3;j++) cout<<tofEffPars[i][j]<<" ";
		//cout<<endl;
	}
	fTofEff = new TF1("fTofEff","[0]*exp(-pow([1]/x,[2]))",0.2,15);
	inf.close();

	hnEvent->Fill(11);

	//+-------------+
	//| loop events |
	//+-------------+
	Int_t n = chain->GetEntries();
	cout<<n<<" events"<<endl;
	Double_t process,process1,process2;
	for(int i=0;i<n;i++){ 
		process1=i;
		process2=n;
		process = process1/process2*100;
		//if(i%10000==0) cout << "begin " << process << "\% entry...." << endl;
//		cout << "begin " << i << " th entry...." << endl;
		eventtree->GetEntry(i);

		current_nEmce = 0;
		current_nTpce = 0;
		current_nJpsi = 0;
		current_nJpsiSide = 0;
		current_nJpsiLS = 0;
		current_nHadron = 0;
		current_nTrigE = 0;

		tofMult = 0;

		for(int j=0;j<4;j++) nHadrons[j] = 0;

		isBHT0      = kFALSE;
		isBHT1      = kFALSE;
		isBHT2      = kFALSE;
		current_centrality = 0;

		for(Int_t m=0;m<MaxEmcElectron;m++)
		{
			current_EmcHits[m] = kTRUE;
			current_TpcHits[m] = kFALSE;
		}
		hnEvent->Fill(12);

		if(!passEvent(eventtree)) continue;
		zdcRate = eventtree->zdcRate;
		if(eventtree->vz<mVzCut[0]||eventtree->vz>mVzCut[1]) continue;
		for(int j=0;j<eventtree->npTracks;j++) passTrack(eventtree, j);

		vzBufferPointer=int((eventtree->vz-mVzCut[0])/(mVzCut[1]-mVzCut[0])*nVzBin);
		if(vzBufferPointer<0||vzBufferPointer>=nVzBin) {
			cout<<"vz = "<<eventtree->vz<<"  out of range!!!!"<<endl;
			continue;
		}

		cenBufferPointer= current_centrality;
		hNe->Fill(current_nEmce,current_nTpce);
		if(current_nEmce+current_nTpce<2) continue;
		
		makeRealPairs();
	}//end of event loop	
	cout<<"Writing histograms into "<<outFile<<".ana.root"<<endl;
	writeHistograms(outFile);
	delete chain;
	printCuts();
	cout<<"end of program"<<endl; 
	return 0;
} 
//===============================================================================================
void makeRealPairs()
{

	TLorentzVector jpsi(0,0,0,0);
	//+-------------+
	//| TPCe + TPCe |
	//+-------------+
	if(!isBHT0&&!isBHT1&&!isBHT2){
		for(int i=0;i<current_nTpce;i++) {
			for(int j=i+1;j<current_nTpce;j++) {
				if(current_tpceId[i]==current_tpceId[j]) continue;
				if(mPairingMode==TwoEID){
					if(!current_tpceTofTag[i]||!current_tpceTofTag[j]) continue;
				}
				if(mPairingMode==OneEID){
					if(!current_tpceTofTag[i]&&!current_tpceTofTag[j]) continue;
				}
				jpsi = current_tpce[i] + current_tpce[j];
				if(jpsi.Rapidity()>mPairYCut[0]&&jpsi.Rapidity()<=mPairYCut[1]) {
					if(current_tpceq[i]+current_tpceq[j]==0) {//unlike-sign
						Polarization(current_tpceq[i],current_tpceq[j],current_tpce[i],current_tpce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
						fill3DHistograms("unlike",jpsi,i,j,0);							
						//tree->Fill();	
						if(current_centrality==0){ 
							hJpsiMassPt->Fill(jpsi.M(),jpsi.Pt());
              				hNumJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
							if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
								current_jpsi[current_nJpsi]=jpsi;
								current_jpsiDau1Id[current_nJpsi] = current_tpceId[i];
								current_jpsiDau2Id[current_nJpsi] = current_tpceId[j];
								current_nJpsi++;
								fillHistograms("unlike",jpsi);
							}
						}
					} else if(current_tpceq[i]+current_tpceq[j]==2) {//like-sign positive
						Polarization(current_tpceq[i],current_tpceq[j],current_tpce[i],current_tpce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
						fill3DHistograms("like",jpsi,i,j,0);	
						if(current_centrality==0){ 
							hJpsiMassPtBG->Fill(jpsi.M(),jpsi.Pt());
              				hDenJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
							if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
								current_jpsiLS[current_nJpsiLS]=jpsi;
								current_jpsiLSDau1Id[current_nJpsiLS] = current_tpceId[i];
								current_jpsiLSDau2Id[current_nJpsiLS] = current_tpceId[j];
								current_nJpsiLS++;
								fillHistograms("like",jpsi);
							}
						}
					} else if(current_tpceq[i]+current_tpceq[j]==-2) {//like-sign negative
						Polarization(current_tpceq[i],current_tpceq[j],current_tpce[i],current_tpce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
						fill3DHistograms("like",jpsi,i,j,0);	
						if(current_centrality==0){ 
							hJpsiMassPtBG->Fill(jpsi.M(),jpsi.Pt());
              				hDenJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
							if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
								current_jpsiLS[current_nJpsiLS]=jpsi;
								current_jpsiLSDau1Id[current_nJpsiLS] = current_tpceId[i];
								current_jpsiLSDau2Id[current_nJpsiLS] = current_tpceId[j];
								current_nJpsiLS++;
								fillHistograms("like",jpsi);
							}
						}
					}
				}
			}//end of TPC e loop
		}//end of TPC e loop
	}

	//+-------------+
	//| EMCe + TPCe |
	//+-------------+
	for(int i=0;i<current_nEmce;i++) {
		for(int j=0;j<current_nTpce;j++) {
			if(current_emceId[i]==current_tpceId[j]) continue;
			if(current_emcIsTrg[i]<0) continue;
			if(mPairingMode==TwoEID){
				if(!current_tpceTofTag[j]) continue;
			}
			jpsi = current_emce[i] + current_tpce[j];
			if(jpsi.Rapidity()>mPairYCut[0]&&jpsi.Rapidity()<=mPairYCut[1]) {
				if(current_emceq[i]+current_tpceq[j]==0) {//unlike-sign
					Polarization(current_emceq[i],current_tpceq[j],current_emce[i],current_tpce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
					fill3DHistograms("unlike",jpsi,i,j,1);						
					//tree->Fill();	
					if(current_centrality==0){ 
						hJpsiMass->Fill(jpsi.M());
						hJpsiMassPt->Fill(jpsi.M(),jpsi.Pt());
              			hNumJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
						if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
							current_jpsi[current_nJpsi]=jpsi;
							current_jpsiDau1Id[current_nJpsi] = current_emceId[i];
							current_jpsiDau2Id[current_nJpsi] = current_tpceId[j];
							current_nJpsi++;
							fillHistograms("unlike",jpsi);
						}
					}
				} else if(current_emceq[i]+current_tpceq[j]==2) {//like-sign positive
					Polarization(current_emceq[i],current_tpceq[j],current_emce[i],current_tpce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
					fill3DHistograms("like",jpsi,i,j,1);	
					if(current_centrality==0){ 
						hJpsiMassBG->Fill(jpsi.M());
						hJpsiMassPtBG->Fill(jpsi.M(),jpsi.Pt());
              			hDenJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
						if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
							current_jpsiLS[current_nJpsiLS]=jpsi;
							current_jpsiLSDau1Id[current_nJpsiLS] = current_emceId[i];
							current_jpsiLSDau2Id[current_nJpsiLS] = current_tpceId[j];
							current_nJpsiLS++;
							fillHistograms("like",jpsi);
						}
					}
				} else if(current_emceq[i]+current_tpceq[j]==-2) {//like-sign negative
					Polarization(current_emceq[i],current_tpceq[j],current_emce[i],current_tpce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
					fill3DHistograms("like",jpsi,i,j,1);	
					if(current_centrality==0){ 
						hJpsiMassBG->Fill(jpsi.M());
						hJpsiMassPtBG->Fill(jpsi.M(),jpsi.Pt());
              			hDenJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
						if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
							current_jpsiLS[current_nJpsiLS]=jpsi;
							current_jpsiLSDau1Id[current_nJpsiLS] = current_emceId[i];
							current_jpsiLSDau2Id[current_nJpsiLS] = current_tpceId[j];
							current_nJpsiLS++;
							fillHistograms("like",jpsi);
						}
					}
				}
			}
		}//end of TPC e loop
	}//end of EMC e loop

	//+-----------+
	//| EMC + EMC |
	//+-----------+
	for(int i=0;i<current_nEmce;i++) {
		for(int j=i+1;j<current_nEmce;j++) {
			if(current_emceId[i]==current_emceId[j]) continue;
			//	if(current_emcIsTrg[i]>=0&&current_emcIsTrg[j]>=0) continue;
			if(current_emcIsTrg[i]<0&&current_emcIsTrg[j]<0) continue;
			jpsi = current_emce[i] + current_emce[j];
			if(jpsi.Rapidity()>mPairYCut[0]&&jpsi.Rapidity()<=mPairYCut[1]) {
				if(current_emceq[i]+current_emceq[j]==0) {//unlike-sign
					Polarization(current_emceq[i],current_emceq[j],current_emce[i],current_emce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
					fill3DHistograms("unlike",jpsi,i,j,2);	
					//tree->Fill();	
					if(current_centrality==0){ 
						hJpsiMass->Fill(jpsi.M());
						hJpsiMassPt->Fill(jpsi.M(),jpsi.Pt());
              			hNumJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
						if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
							current_jpsi[current_nJpsi]=jpsi;
							current_jpsiDau1Id[current_nJpsi] = current_emceId[i];
							current_jpsiDau2Id[current_nJpsi] = current_emceId[j];
							current_nJpsi++;
							fillHistograms("unlike",jpsi);
						}
					}

				} else if(current_emceq[i]+current_emceq[j]==2) {//like-sign positive
					Polarization(current_emceq[i],current_emceq[j],current_emce[i],current_emce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
					fill3DHistograms("like",jpsi,i,j,2);	
					if(current_centrality==0){ 
						hJpsiMassBG->Fill(jpsi.M());
						hJpsiMassPtBG->Fill(jpsi.M(),jpsi.Pt());
              			hDenJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
						if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
							current_jpsiLS[current_nJpsiLS]=jpsi;
							current_jpsiLSDau1Id[current_nJpsiLS] = current_emceId[i];
							current_jpsiLSDau2Id[current_nJpsiLS] = current_emceId[j];
							current_nJpsiLS++;
							fillHistograms("like",jpsi);
						}
					}
				} else if(current_emceq[i]+current_emceq[j]==-2) {//like-sign negative
					Polarization(current_emceq[i],current_emceq[j],current_emce[i],current_emce[j],hJpsiCosThetaPt,hJpsiPhiPt,hJpsiCosThetaPtCS,hJpsiPhiPtCS,1);
					fill3DHistograms("like",jpsi,i,j,2);	
					if(current_centrality==0){ 
						hJpsiMassBG->Fill(jpsi.M());
						hJpsiMassPtBG->Fill(jpsi.M(),jpsi.Pt());
              			hDenJpsiPtMTofMult->Fill(tofMult,jpsi.Pt(),jpsi.M());
						if(jpsi.M()>mJpsiMassCut[0]&&jpsi.M()<mJpsiMassCut[1]){
							current_jpsiLS[current_nJpsiLS]=jpsi;
							current_jpsiLSDau1Id[current_nJpsiLS] = current_emceId[i];
							current_jpsiLSDau2Id[current_nJpsiLS] = current_emceId[j];
							current_nJpsiLS++;
							fillHistograms("like",jpsi);
						}
					}
				}
			}
		}//end of TPC e loop
	}//end of EMC e loop
}

void printCuts() 
{
	cout<<endl;
	cout<<"------- print all configures ----------"<<endl;
	cout<<"Event Cuts:"<<endl;
	cout<<mVzCut[0]<<" < vzTpc < "<<mVzCut[1]<<endl;
	//cout<<"|vzVpd - vzTpc| < "<<mVzDiffCut<<endl;
	//cout<<"ranking > "<<mRankingCut<<endl;
	cout<<"trigger Ids: ";
	for(int i=0;i<nTrigIds;i++){
		printf("0x%X \t",mTriggerId[i]);
	}
	cout<<endl;

	cout<<"Track Cuts:"<<endl;
	cout<<"+++TPC+++"<<endl;
	cout<<mPtCut[0]<<" < pT < "<<mPtCut[1]<<endl;
	cout<<mEtaCut[0]<<" < eta < "<<mEtaCut[1]<<endl;
	cout<<"0 < Dca < "<<mDcaCut<<endl;
	cout<<"nHitsFit >= "<<mHitsFitCut<<endl;
	cout<<"nHitsdEdx >= "<<mHitsDedxCut<<endl;
	cout<<"nHitsRatio > "<<mHitsRatioCut<<endl;
	cout<<mnSigECut[0]<<" < nSigmaE < "<<mnSigECut[1]<<endl;
	cout<<mnSigELooseCut[0]<<" < loose nSigmaE < "<<mnSigELooseCut[1]<<endl;

	cout<<"+++TOF+++"<<endl;
	cout<<"|yLocal| < "<<myLocalCut<<endl;
	cout<<mBetaCut[0]<<" < 1/beta < "<<mBetaCut[1]<<endl;

	cout<<"+++EMC+++"<<endl;
	cout<<mPveCut[0]<<" < p/E < "<<mPveCut[1]<<endl;

	cout<<"+++Pair Level+++"<<endl;
	cout<<mPairYCut[0]<<" < yee < "<<mPairYCut[1]<<endl;
	cout<<"--------- end of configures -----------"<<endl;
	cout<<endl;
}
//================================================================================================
bool passEvent(pevent* event)
{
	UShort_t triggerId = event->triggerId;
	int htId[nTrigIds];
	for(int i=0;i<nTrigIds;i++){
		htId[i] = -1;
		if(triggerId & mTriggerId[i]) htId[i] = i;
	}
	//0, bht0-vpdmb; 1, bht1-vpdmb; 2, bht2; 3, bht0-bbcmb-tof0; 4, bht1-bbcmb-tof0; 5, bht2-bbcmb; 
	float ranking = event->ranking;
	//if(ranking<=0) return false;

	int run = event->runnumber;
	hnEvent->Fill(0);
	for(int i=0;i<47;i++){
		if(run==mBadRuns[i]) return false;
	}
	hnEvent->Fill(1);
	runPointer = TMath::BinarySearch(nRuns, runIds, run);
	if(runPointer<0) {cout<<"run "<<run<<"    not in the list ! Ignore this run"<<endl;}

	//cout<<"triggerId = "<<triggerId<<" htId = "<<htId<<endl;
	if(mHTselect==-1){
		for(int i=0;i<nMBTrigIds;i++){
			if(triggerId & mMBTriggerId[i]) return true;
		}	
	}
	if(mHTselect==0){
		//		if(htId[0]==0&&htId[3]!=3){ isBHT0=kTRUE; return true;}
		if(htId[0]==0&&htId[1]!=1&&htId[2]!=2&&htId[3]!=3&&htId[4]!=4&&htId[5]!=5){ isBHT0=kTRUE; return true;}
	}
	if(mHTselect==1){
		//		if(htId[1]==1&&htId[4]!=4){ isBHT1=kTRUE; return true;}
		if(htId[1]==1&&htId[2]!=2&&htId[3]!=3&&htId[4]!=4&&htId[5]!=5){ isBHT1=kTRUE; return true;}
	}
	if(mHTselect==2){
		//	if(htId[2]==2&&htId[5]!=5){ isBHT2=kTRUE; return true;}
		if(htId[2]==2&&htId[3]!=3&&htId[4]!=4&&htId[5]!=5){ isBHT2=kTRUE; return true;}
	}
	if(mHTselect==3){
		if(htId[3]==3){ isBHT0=kTRUE; return true;}
//		if(htId[3]==3&&htId[4]!=4&&htId[5]!=5){ isBHT0=kTRUE; return true;}
	}
	if(mHTselect==4){
		if(htId[4]==4){ isBHT1=kTRUE; return true;}
//		if(htId[4]==4&&htId[5]!=5){ isBHT1=kTRUE; return true;}
	}
	if(mHTselect==5){
		if(htId[5]==5){ isBHT2=kTRUE; return true;}
	}
	hnEvent->Fill(4);
	return false;
}

//========================================================================
bool passTrack(pevent* event, Int_t i)
{
	Int_t q = (event->track_NHitsFit[i]>0) ? 1 : -1;
	checkcut->Fill(0);
	if(q!=1&&q!=-1) return kFALSE;
	hnTrack->Fill(0);
	Int_t nHitsFit = abs(event->track_NHitsFit[i]);
	Int_t nHitsDedx = event->track_NHitsDedx[i];
	Float_t beta = event->track_mBTofBeta[i]/20000.;
	Float_t ylocal = event->track_mBTofylocal[i]/1000.;
	Float_t nSigmaE = event->track_nSigmaE[i]/100.;
	Float_t nSigmaPi = event->track_nSigmaPi[i]/100.;
	Float_t nSigmaK = event->track_nSigmaK[i]/100.;
	Float_t nSigmaP = event->track_nSigmaP[i]/100.;
	Float_t px = event->track_px[i];
	Float_t py = event->track_py[i];
	Float_t pz = event->track_pz[i];
	if(sqrt(px*px+py*py)<0.2) return kFALSE;
//	Float_t pve = event->track_pve[i]/1000.;
	Float_t dca = event->track_dca[i]/1000.;
	int tofmatchflag = 0;
	if(beta>0) tofmatchflag = 1;
	Int_t trackId = i;
	TVector3 mom(0,0,0);
	mom.SetXYZ(px,py,pz);
	if(fabs(ylocal)<10) tofMult++;

	//cout<<"mom = "<<px<<","<<py<<","<<pz<<endl;
	Float_t pt = mom.Perp();
	Float_t p  = mom.Mag();
	Float_t phi = mom.Phi();
	Float_t eta = mom.Eta();
//	Char_t tpcFlag = event->track_tpcFlag[i];
//	Char_t tofFlag = event->track_tofFlag[i];
//	Char_t emcFlag = event->track_emcFlag[i];
//	Char_t smdFlag = event->track_smdFlag[i];

	Float_t energy = 0;
	energy = event->track_e0[i]/1000.;
	Float_t pve = energy>0?p/energy:0;
	//	if(pve>0) energy = p/pve;
	Float_t adc = event->track_adc0[i];
	Float_t dsmadc = event->track_dsmAdc[i];
	Short_t towerId = event->track_towerId[i];

	int hotTowFlag = 0;
	for(int k=0;k<106;k++){
		if(towerId==dsmHotTowers[k]) {hotTowFlag = 1; break;}
	}
	checkcut->Fill(1);
	if(eta<=mEtaCut[0] || eta>=mEtaCut[1]) return kFALSE;
	checkcut->Fill(2);
	if((nHitsFit<mEmcHitsFitCut[0] && mNHitsFit!=20 && mNHitsFit!=21 && mNHitsFit!=22 && mNHitsFit!=23) 
			||(nHitsFit<mEmcHitsFitCut[1] && mNHitsFit==20)
			||(nHitsFit<mEmcHitsFitCut[2] && mNHitsFit==21)
			||(nHitsFit<mEmcHitsFitCut[3] && mNHitsFit==22)
			||(nHitsFit<mEmcHitsFitCut[4] && mNHitsFit==23)) return kFALSE;
	checkcut->Fill(3);
	hnTrack->Fill(4);
	if(nHitsDedx<mEmcHitsDedxCut) return kFALSE;
	checkcut->Fill(4);
	hnhitsDEDX->Fill(nHitsDedx);	
	if(dca>=mDcaCut) return kFALSE;
	hnTrack->Fill(5);
	checkcut->Fill(5);
	hDCA->Fill(dca);

	if(nSigmaE>mnSigELooseCut[0]&&nSigmaE<mnSigELooseCut[1]) {//select loose electron
		int isEmcTrk = 0;
		checkcut->Fill(6);
		hnTrack->Fill(8);
//    cout<<"nSigmaE = "<<nSigmaE<<" pve = "<<pve<<" pt= "<<pt<<" mnSigECut[0]"<<mnSigECut[0]<<"   "<<mnSigECut[1]<<"mEmcMnPt"<<mEmcMinPt<<endl;
		if(nSigmaE>mnSigECut[0]&&nSigmaE<mnSigECut[1]&&((pve>mPveCut[0]&&pve<mPveCut[1] && mPOE==0)||(pve>mPveCut[2]&&pve<mPveCut[3] && mPOE==1))&&pt>mEmcMinPt && hotTowFlag==0){
//			hnsigmaE->Fill(nSigmaE);
					hnTrack->Fill(11);	
			hemcPtMin->Fill(pt);
			current_emce[current_nEmce].SetPtEtaPhiM(pt,eta,phi,Melectron);
			current_emceId[current_nEmce] = trackId;
			current_emceq[current_nEmce] = q;
			current_emcIsTrg[current_nEmce] = -1;
			current_emcdsmadc[current_nEmce] = dsmadc;
			current_emcadc0[current_nEmce] = adc;
//			if(isBHT0 && pt>2.5 && adc>mAdc0Cut[0] && dsmadc>11*dsmadcfactor) current_emcIsTrg[current_nEmce]  = 0;
//			if(isBHT1 && pt>3.6 && adc>mAdc0Cut[1] && dsmadc>15*dsmadcfactor) current_emcIsTrg[current_nEmce]  = 1;
//			if(isBHT2 && pt>4.3 && adc>mAdc0Cut[2] && dsmadc>18*dsmadcfactor) current_emcIsTrg[current_nEmce]  = 2;
	
			if(isBHT0 && pt>2.5 && adc>mAdc0Cut[0]*adcfactor) current_emcIsTrg[current_nEmce]  = 0;
			if(isBHT1 && pt>3.6 && adc>mAdc0Cut[1]*adcfactor) current_emcIsTrg[current_nEmce]  = 1;
		checkcut->Fill(7);
			if(isBHT2 && pt>4.3 && adc>mAdc0Cut[2]*adcfactor) current_emcIsTrg[current_nEmce]  = 2;
		checkcut->Fill(8);
/*	
			if(isBHT0 && pt>2.5 && dsmadc>11) current_emcIsTrg[current_nEmce]  = 0;
			if(isBHT1 && pt>3.6 && dsmadc>15) current_emcIsTrg[current_nEmce]  = 1;
		checkcut->Fill(7);
			if(isBHT2 && pt>4.3 && dsmadc>18) current_emcIsTrg[current_nEmce]  = 2;
		checkcut->Fill(8);
*/
			if(current_emcIsTrg[current_nEmce]>=0){ 
				current_nTrigE++;
			}
			current_nEmce++;
			isEmcTrk = 1;
		}

		//exclude Emc track
		if(!isEmcTrk){
			current_tpce[current_nTpce].SetPtEtaPhiM(pt,eta,phi,Melectron);
			current_tpceId[current_nTpce] = trackId;
			current_tpceq[current_nTpce] = q;
			current_tpceTofTag[current_nTpce] = 0;
			current_tpcdsmadc[current_nTpce] = dsmadc;
			current_tpcadc0[current_nTpce] = adc;
			if(nSigmaE>mnSigECut[0]&&nSigmaE<mnSigECut[1]&&1./beta>mBetaCut[0]&&1./beta<mBetaCut[1]&&fabs(ylocal)<myLocalCut){
				current_tpceTofTag[current_nTpce] = 1;
				hnTrack->Fill(15);
			}
			current_nTpce++;
		checkcut->Fill(9);
		}
	}
	return kTRUE;
}
//========================================================================
void bookHistograms(char* outFile){
	char buf[1024];
	sprintf(buf,"%s.ana.root",outFile);
	mFile = new TFile(buf,"recreate");
	TH1::SetDefaultSumw2();

	hnTrack = new TH1F("hnTrack","hnTrack",30,0,30);
	hnEvent = new TH1F("hnEvent","hnEvent",30,0,30);
  	int nZdcXBins = 100;
  float zdcXMin = 0;
  float zdcXMax = 20; // k

  int nBbcXBins = 150;
  float bbcXMin = 0;
  float bbcXMax = 1500; // K 

  int nTrkBins = 100;
  float nTrkMin = 0;
  float nTrkMax = 100;

  int nJpsiBins = 100;
  float nJpsiMin = -10;
  float nJpsiMax = 10;

  int nPtBins = 60;
  float ptMin = 0;
  float ptMax = 15;

  int nMassBins = 40; 
  float massMin = 2.; 
  float massMax = 4.;

  float nTpcT = 3.5;

  int nPvEBins = 100;
  float pveMin = 0;
 
	checkcut = new TH1F("checkcut","checkcut",10,-0.5,9.5);
  	hNe = new TH2F("hNe","#EMC e vs. #TPC e;#EMC e candidate;#TPC e candidate;Counts",50,0,50,50,0,50);
	char label[10][100] = {"Wrong","BHT0","BHT1*VPDMB","BHT2","BHT0*BBCMB*TOF0","BHT1*BBCMB*TOF0","BHT2*BBCMB","","",""};

	hJpsiPhiPt = new TH2F("hJpsiPhiPt","Jpsi Phi vs #Phi;#Phi;P_{T} GeV/c",360,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiPhiPtBG = new TH2F("hJpsiPhiPtBG","Jpsi Phi vs #Phi;#Phi;P_{T} GeV/c",360,-TMath::Pi(),TMath::Pi(),120,0,30);

	hJpsiCosThetaPt = new TH2F("hJpsiCosThetaPt","Jpsi Pt vs Cos(#theta); Cos(#theta); Jpsi Pt;",10,-1,1,120,0,30);
	hJpsiCosThetaPt->Sumw2();
	hJpsiCosThetaPtBG = new TH2F("hJpsiCosThetaPtBG","Jpsi Pt vs Cos(#theta); Cos(#theta); Jpsi Pt;",10,-1,1,120,0,30);
	hJpsiCosThetaPtBG->Sumw2();

	hJpsiPt = new TH1F("hJpsiPt","Jpsi Pt distribution; Jpsi Pt",120,0.,30.);
	hJpsiPtBG = new TH1F("hJpsiPtBG","Jpsi Pt distribution; Jpsi Pt",120,0.,30.);

	hJpsiMass = new TH1F("hJpsiMass","Invariant mass distribution ; M_{ee}",50,2,4);
	hJpsiMassBG = new TH1F("hJpsiMassBG","Invariant mass distribution ; M_{ee}",50,2,4);

	hJpsiEta = new TH1F("hJpsiEta","Jpsi Eta distribution; Jpsi #eta",30,-1.5,1.5);
	hJpsiEtaBG = new TH1F("hJpsiEtaBG","Jpsi Eta distribution; Jpsi #eta",30,-1.5,1.5);

	hJpsiRapidity = new TH1F("hJpsiRapidity","Jpsi Rapidity distribution; Jpsi Rapidity",30,-1.5,1.5);
	hJpsiRapidityBG = new TH1F("hJpsiRapidityBG","Jpsi Rapidity distribution; Jpsi Rapidity",30,-1.5,1.5);

	hnsigmaE = new TH1F("hnsigmaE","hnsigmaE",100,-5,5);
	hemcPtMin = new TH1F("hemcPtMin","hemcPtMin",100,0,10);
	hDCA = new TH1F("hDCA","hDCA",100,0,5);
	hnhitsDEDX = new TH1F("hnhitsDEDX","hnhitsDEDX",100,0,40);

	hJpsiMassPt = new TH2F("hJpsiMassPt","Jpsi Pt vs mass distribution; Jpsi mass; Jpsi Pt",40,2,4,120,0,30);
	hJpsiMassPt->Sumw2();
	hJpsiMassPtBG = new TH2F("hJpsiMassPtBG","Jpsi Pt vs mass distribution background; Jpsi mass; Jpsi Pt",40,2,4,120,0,30);
	hJpsiMassPtBG->Sumw2();

	hJpsiPhiPtHX = new TH2F("hJpsiPhiPtHX","Jpsi Pt vs #phi;#phi;Jpsi Pt",10,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiPhiPtHX->Sumw2();
	hJpsiPhiPtHXBG = new TH2F("hJpsiPhiPtHXBG","Jpsi Pt vs #phi;#phi;Jpsi Pt",10,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiPhiPtHXBG->Sumw2();

	hJpsiCosThetaPtCS = new TH2F("hJpsiCosThetaPtCS","Jpsi Pt vs Cos(#theta);Cos(#theta);Jpsi Pt",10,-1,1,120,0,30);
	hJpsiCosThetaPtCS->Sumw2();
	hJpsiCosThetaPtCSBG = new TH2F("hJpsiCosThetaPtCSBG","Jpsi Pt vs Cos(#theta);Cos(#theta);Jpsi Pt",10,-1,1,120,0,30);
	hJpsiCosThetaPtCSBG->Sumw2();

	hJpsiPhiPtCS = new TH2F("hJpsiPhiPtCS","Jpsi Pt vs #phi;#phi;Jpsi Pt",10,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiPhiPtCS->Sumw2();
	hJpsiPhiPtCSBG = new TH2F("hJpsiPhiPtCSBG","Jpsi Pt vs #phi;#phi;Jpsi Pt",10,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiPhiPtCSBG->Sumw2();

	hJpsiCosThetaInvMPt = new TH3F("hJpsiCosThetaInvMPt","hJpsiCosThetaInvMPt",40,-1,1,40,2,4,120,0,30);
	hJpsiCosThetaInvMPtBG= new TH3F("hJpsiCosThetaInvMPtBG","hJpsiCosThetaInvMPtBG",40,-1,1,40,2,4,120,0,30);
	hJpsiCosThetaInvMPtCS= new TH3F("hJpsiCosThetaInvMPtCS","hJpsiCosThetaInvMPtCS",40,-1,1,40,2,4,120,0,30);
	hJpsiCosThetaInvMPtCSBG= new TH3F("hJpsiCosThetaInvMPtCSBG","hJpsiCosThetaInvMPtCSBG",40,-1,1,40,2,4,120,0,30);

	hJpsiPhiInvMPt = new TH3F("hJpsiPhiInvMPt","hJpsiPhiInvMPt",40,-TMath::Pi(),TMath::Pi(),40,2,4,120,0,30);
	hJpsiPhiInvMPtBG = new TH3F("hJpsiPhiInvMPtBG","hJpsiPhiInvMPtBG",40,-TMath::Pi(),TMath::Pi(),40,2,4,120,0,30);
	hJpsiPhiInvMPtCS = new TH3F("hJpsiPhiInvMPtCS","hJpsiPhiInvMPtCS",40,-TMath::Pi(),TMath::Pi(),40,2,4,120,0,30);
	hJpsiPhiInvMPtCSBG = new TH3F("hJpsiPhiInvMPtCSBG","hJpsiPhiInvMPtCSBG",40,-TMath::Pi(),TMath::Pi(),40,2,4,120,0,30);

	hJpsiCosThetaInvMPt->Sumw2();
	hJpsiCosThetaInvMPtBG->Sumw2();
	hJpsiCosThetaInvMPtCS->Sumw2();
	hJpsiCosThetaInvMPtCSBG->Sumw2();

	hJpsiPhiInvMPt->Sumw2();
	hJpsiPhiInvMPtBG->Sumw2();
	hJpsiPhiInvMPtCS->Sumw2();
	hJpsiPhiInvMPtCSBG->Sumw2();

  hNumJpsiPtMTofMult = new TH3F("hNumJpsiPtMTofMult","J/#psi pT M vs tofMult; tofMult; J/#psi pT (GeV/c); M_{ee} (GeV/c^{2})",nTrkBins,nTrkMin,nTrkMax,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);

  hDenJpsiPtMTofMult = new TH3F("hDenJpsiPtMTofMult","J/#psi pT M vs tofMult; tofMult; J/#psi pT (GeV/c); M_{ee} (GeV/c^{2})",nTrkBins,nTrkMin,nTrkMax,nPtBins,ptMin,ptMax,nMassBins,massMin,massMax);
//	dsmAdcInvMPt = new TH3F("dsmAdcInvMPt","dsmAdcInvMPt;dsmAdc;M_{ee};p_{T}",65,0,65,40,2,4,120,0,30);
//	dsmAdcInvMPtBG = new TH3F("dsmAdcInvMPtBG","dsmAdcInvMPtBG;dsmAdc;M_{ee};p_{T}",65,0,65,40,2,4,120,0,30);
//	AdcInvMPt = new TH3F("AdcInvMPt","AdcInvMPt;adc;M_{ee};p_{T}",800,0,800,40,2,4,120,0,30);
//	AdcInvMPtBG = new TH3F("AdcInvMPtBG","AdcInvMPtBG;adc;M_{ee};p_{T}",800,0,800,40,2,4,120,0,30);

	//dsmAdcInvMPt->Sumw2();
	//dsmAdcInvMPtBG->Sumw2();
	//AdcInvMPt->Sumw2();
//	AdcInvMPtBG->Sumw2();

	hJpsiCosThetaPhiPt = new TH3F("hJpsiCosThetaPhiPt","hJpsiCosThetaPhiPt;cos#theta;#phi;p_{T}",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiCosThetaPhiPtBG = new TH3F("hJpsiCosThetaPhiPtBG","hJpsiCosThetaPhiPtBG;cos#theta;#phi;p_{T}",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiCosThetaPhiPt->Sumw2();
	hJpsiCosThetaPhiPtBG->Sumw2();

	hJpsiCosThetaPhiPtCS = new TH3F("hJpsiCosThetaPhiPtCS","hJpsiCosThetaPhiPtCS;cos#theta;#phi;p_{T}",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiCosThetaPhiPtCSBG = new TH3F("hJpsiCosThetaPhiPtCSBG","hJpsiCosThetaPhiPtCSBG;cos#theta;#phi;p_{T}",40,-1,1,40,-TMath::Pi(),TMath::Pi(),120,0,30);
	hJpsiCosThetaPhiPtCS->Sumw2();
	hJpsiCosThetaPhiPtCSBG->Sumw2();


	tree = new TTree("tree","J/psi polarization");
	tree->SetAutoSave(100000);
	tree->Branch("jpsi_pt",&jpsi_pt,"jpsi_pt/F");
	tree->Branch("jpsi_eta",&jpsi_eta,"jpsi_eta/F");
	tree->Branch("jpsi_phi",&jpsi_phi,"jpsi_phi/F");
	tree->Branch("jpsi_InvM",&jpsi_InvM,"jpsi_InvM/F");
	tree->Branch("lepton1_pt",&lepton1_pt,"lepton1_pt/F");
	tree->Branch("lepton1_eta",&lepton1_eta,"lepton1_eta/F");
	tree->Branch("lepton1_phi",&lepton1_phi,"lepton1_phi/F");
	tree->Branch("lepton1_InvM",&lepton1_InvM,"lepton1_InvM/F");
	tree->Branch("lepton2_pt",&lepton2_pt,"lepton2_pt/F");
	tree->Branch("lepton2_eta",&lepton2_eta,"lepton2_eta/F");
	tree->Branch("lepton2_phi",&lepton2_phi,"lepton2_phi/F");
	tree->Branch("lepton2_InvM",&lepton2_InvM,"lepton2_InvM/F");
	tree->Branch("positron_theta_hx",&positron_theta_hx,"positron_theta_hx/F");
	tree->Branch("positron_theta_cs",&positron_theta_cs,"positron_theta_cs/F");
	tree->Branch("positron_phi_hx",&positron_phi_hx,"positron_phi_hx/F");
	tree->Branch("positron_phi_cs",&positron_phi_cs,"positron_phi_cs/F");
	tree->Branch("electron_theta_hx",&electron_theta_hx,"electron_theta_hx/F");
	tree->Branch("electron_theta_cs",&electron_theta_cs,"electron_theta_cs/F");
	tree->Branch("electron_phi_hx",&electron_phi_hx,"electron_phi_hx/F");
	tree->Branch("electron_phi_cs",&electron_phi_cs,"electron_phi_cs/F");
	tree->Branch("lepton1_isTPCe",&lepton1_isTPCe,"lepton1_isTPCe/B");
	tree->Branch("lepton2_isTPCe",&lepton2_isTPCe,"lepton2_isTPCe/B");

}
//===============================================================================================

void fillHistograms(std::string unlikeOrlike, TLorentzVector JPSI)
{
	if(unlikeOrlike.compare("unlike")==0){
		hJpsiCosThetaPt->Fill(TMath::Cos(positron_theta_hx),JPSI.Pt());
		hJpsiPhiPtHX->Fill(positron_phi_hx,JPSI.Pt());
		hJpsiCosThetaPtCS->Fill(TMath::Cos(positron_theta_cs),JPSI.Pt());
		hJpsiPhiPtCS->Fill(positron_phi_cs,JPSI.Pt());
		hJpsiCosThetaPhiPt->Fill(TMath::Cos(positron_theta_hx),positron_phi_hx,JPSI.Pt());	
		hJpsiCosThetaPhiPtCS->Fill(TMath::Cos(positron_theta_cs),positron_phi_cs,JPSI.Pt());	
	}
	else{
		hJpsiCosThetaPtBG->Fill(TMath::Cos(positron_theta_hx),JPSI.Pt(),0.5);
		hJpsiPhiPtHXBG->Fill(positron_phi_hx,JPSI.Pt(),0.5);
		hJpsiCosThetaPtCSBG->Fill(TMath::Cos(positron_theta_cs),JPSI.Pt(),0.5);
		hJpsiPhiPtCSBG->Fill(positron_phi_cs,JPSI.Pt(),0.5);
		hJpsiCosThetaPtBG->Fill(TMath::Cos(electron_theta_hx),JPSI.Pt(),0.5);
		hJpsiPhiPtHXBG->Fill(electron_phi_hx,JPSI.Pt(),0.5);
		hJpsiCosThetaPtCSBG->Fill(TMath::Cos(electron_theta_cs),JPSI.Pt(),0.5);
		hJpsiPhiPtCSBG->Fill(electron_phi_cs,JPSI.Pt(),0.5);
		hJpsiCosThetaPhiPtBG->Fill(TMath::Cos(positron_theta_hx),positron_phi_hx,JPSI.Pt(),0.5);
		hJpsiCosThetaPhiPtBG->Fill(TMath::Cos(electron_theta_hx),electron_phi_hx,JPSI.Pt(),0.5);
		hJpsiCosThetaPhiPtCSBG->Fill(TMath::Cos(positron_theta_cs),positron_phi_cs,JPSI.Pt(),0.5);
		hJpsiCosThetaPhiPtCSBG->Fill(TMath::Cos(electron_theta_cs),electron_phi_cs,JPSI.Pt(),0.5);
	}
}

void fill3DHistograms(std::string unlikeOrlike, TLorentzVector JPSI,int i,int j,int pairs){
	if(unlikeOrlike.compare("unlike")==0){
		hJpsiCosThetaInvMPt->Fill(TMath::Cos(positron_theta_hx),JPSI.M(),JPSI.Pt());
		hJpsiCosThetaInvMPtCS->Fill(TMath::Cos(positron_theta_cs),JPSI.M(),JPSI.Pt());
		hJpsiPhiInvMPt->Fill(positron_phi_hx,JPSI.M(),JPSI.Pt());
		hJpsiPhiInvMPtCS->Fill(positron_phi_cs,JPSI.M(),JPSI.Pt());
		if(pairs==0){
			//dsmAdcInvMPt->Fill(current_tpcdsmadc[i],JPSI.M(),JPSI.Pt());
			//AdcInvMPt->Fill(current_tpcadc0[i],JPSI.M(),JPSI.Pt());	
			//dsmAdcInvMPt->Fill(current_tpcdsmadc[j],JPSI.M(),JPSI.Pt());
			//AdcInvMPt->Fill(current_tpcadc0[j],JPSI.M(),JPSI.Pt());	
		}
		if(pairs==1){
			if(current_emcIsTrg[i]>=0){	
				//dsmAdcInvMPt->Fill(current_emcdsmadc[i],JPSI.M(),JPSI.Pt());
				//AdcInvMPt->Fill(current_emcadc0[i],JPSI.M(),JPSI.Pt());	
			}
		}
		if(pairs==2){
			if(current_emcIsTrg[i]>=0){	
				//dsmAdcInvMPt->Fill(current_emcdsmadc[i],JPSI.M(),JPSI.Pt());
				//AdcInvMPt->Fill(current_emcadc0[i],JPSI.M(),JPSI.Pt());	
			}
			if(current_emcIsTrg[j]>=0){	
				//dsmAdcInvMPt->Fill(current_emcdsmadc[j],JPSI.M(),JPSI.Pt());
				//AdcInvMPt->Fill(current_emcadc0[j],JPSI.M(),JPSI.Pt());	
			}
		}
	}
	else{	
		hJpsiCosThetaInvMPtBG->Fill(TMath::Cos(positron_theta_hx),JPSI.M(),JPSI.Pt(),0.5);
		hJpsiCosThetaInvMPtBG->Fill(TMath::Cos(electron_theta_hx),JPSI.M(),JPSI.Pt(),0.5);
		hJpsiCosThetaInvMPtCSBG->Fill(TMath::Cos(positron_theta_cs),JPSI.M(),JPSI.Pt(),0.5);
		hJpsiCosThetaInvMPtCSBG->Fill(TMath::Cos(electron_theta_cs),JPSI.M(),JPSI.Pt(),0.5);
		hJpsiPhiInvMPtBG->Fill(positron_phi_hx,JPSI.M(),JPSI.Pt(),0.5);
		hJpsiPhiInvMPtBG->Fill(electron_phi_hx,JPSI.M(),JPSI.Pt(),0.5);
		hJpsiPhiInvMPtCSBG->Fill(positron_phi_cs,JPSI.M(),JPSI.Pt(),0.5);
		hJpsiPhiInvMPtCSBG->Fill(electron_phi_cs,JPSI.M(),JPSI.Pt(),0.5);
		if(pairs==0){
			//dsmAdcInvMPtBG->Fill(current_tpcdsmadc[i],JPSI.M(),JPSI.Pt());
//			dsmAdcInvMPtBG->Fill(current_tpcdsmadc[j],JPSI.M(),JPSI.Pt());
//			AdcInvMPtBG->Fill(current_tpcadc0[i],JPSI.M(),JPSI.Pt());	
//			AdcInvMPtBG->Fill(current_tpcadc0[j],JPSI.M(),JPSI.Pt());	
		}
		if(pairs==1){
			if(current_emcIsTrg[i]>=0){	
				//dsmAdcInvMPtBG->Fill(current_emcdsmadc[i],JPSI.M(),JPSI.Pt());
//				AdcInvMPtBG->Fill(current_emcadc0[i],JPSI.M(),JPSI.Pt());	
			}
		}
		if(pairs==2){
			if(current_emcIsTrg[i]>=0){	
				//dsmAdcInvMPtBG->Fill(current_emcdsmadc[i],JPSI.M(),JPSI.Pt());
				//AdcInvMPtBG->Fill(current_emcadc0[i],JPSI.M(),JPSI.Pt());	
			}
			if(current_emcIsTrg[j]>=0){	
				//dsmAdcInvMPtBG->Fill(current_emcdsmadc[j],JPSI.M(),JPSI.Pt());
				//AdcInvMPtBG->Fill(current_emcadc0[j],JPSI.M(),JPSI.Pt());	
			}
		}
	}
}


void writeHistograms(char* outFile)
{
	mFile->cd();
	mFile->Write();
	mFile->Close();
} 

void setPairingMode(pairMode mode){
	mPairingMode = mode;
}

void Polarization(int icharge,int jcharge,TLorentzVector ivector,TLorentzVector jvector,TH2F *hist0, TH2F *hist1, TH2F *hist2, TH2F *hist3, bool isSignal){
	TLorentzVector mpositron,melectron,JPSI;
	TLorentzVector Proton1(0.,0.,100.,100.),Proton2(0.,0.,-100.,100.);
	TVector3 XX,YY,ZZ;
	JPSI = ivector+jvector;
	mpositron = icharge>=jcharge? ivector:jvector;
	melectron = jcharge<=icharge? jvector:ivector;

	Int_t NFRAME =4;
	Double_t theta[NFRAME];
	Double_t phi[NFRAME];
	TVector3 XXHX,YYHX,ZZHX;
	ZZHX = JPSI.Vect();
	YYHX = JPSI.Vect().Cross(Proton1.Vect());
	XXHX = YYHX.Cross(ZZHX);

	mpositron.Boost(-JPSI.Px()/JPSI.E(),-JPSI.Py()/JPSI.E(),-JPSI.Pz()/JPSI.E());
	melectron.Boost(-JPSI.Px()/JPSI.E(),-JPSI.Py()/JPSI.E(),-JPSI.Pz()/JPSI.E());
	Proton1.Boost(-JPSI.Px()/JPSI.E(),-JPSI.Py()/JPSI.E(),-JPSI.Pz()/JPSI.E());
	Proton2.Boost(-JPSI.Px()/JPSI.E(),-JPSI.Py()/JPSI.E(),-JPSI.Pz()/JPSI.E());

	theta[0]= mpositron.Angle(JPSI.Vect());
	phi[0] = TMath::ATan2(mpositron.Vect().Dot(YYHX.Unit()),mpositron.Vect().Dot(XXHX.Unit()));
	electron_theta_hx = melectron.Angle(JPSI.Vect());
	electron_phi_hx = TMath::ATan2(melectron.Vect().Dot(YYHX.Unit()),(melectron.Vect().Dot(XXHX.Unit())));

	ZZ = Proton1.Vect()*(1/(Proton1.Vect()).Mag())-Proton2.Vect()*(1/(Proton2.Vect()).Mag());
	YY = Proton1.Vect().Cross(Proton2.Vect());
	XX = Proton1.Vect()*(1/(Proton1.Vect()).Mag())+Proton2.Vect()*(1/(Proton2.Vect()).Mag());

	theta[1] = mpositron.Angle(ZZ);
	phi[1] = TMath::ATan2(mpositron.Vect().Dot(YY.Unit()),mpositron.Vect().Dot(XX.Unit()));

	positron_theta_hx = theta[0];
	positron_theta_cs = theta[1];
	positron_phi_hx = phi[0];
	positron_phi_cs = phi[1];

	electron_theta_cs = melectron.Angle(ZZ);
	electron_phi_cs = TMath::ATan2(melectron.Vect().Dot(YY.Unit()),melectron.Vect().Dot(XX.Unit()));

	jpsi_pt = JPSI.Pt();
	jpsi_eta = JPSI.Eta();
	jpsi_phi = JPSI.Phi();
	jpsi_InvM = JPSI.M();

	lepton1_pt = mpositron.Pt();
	lepton1_eta = mpositron.Eta();
	lepton1_phi = mpositron.Phi();
	lepton1_InvM = mpositron.M();

	lepton2_pt = melectron.Pt();
	lepton2_eta = melectron.Eta();
	lepton2_phi = melectron.Phi();
	lepton2_InvM = melectron.M();
}





