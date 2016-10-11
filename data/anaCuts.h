const Double_t mVzCut[2] = {-50., 50.};
const Double_t mVzDiffCut = 6.;
const Double_t mRankingCut = 0;

const int nMBTrigIds = 1;
int mMBTriggerId[nMBTrigIds] = {0x1};
//ht 370501,370511,370531,370542,370546,370522
const int nTrigIds = 6;
int mTriggerId[nTrigIds] = {0x8,0x10,0x20,0x40,0x80,0x100};

const Double_t mHadronDcaCut = 1.;
const Double_t mHadronnHitsFitCut = 15;

const Double_t mEtaCut[2] = {-1, 1};
const Double_t mHadEtaCut[2] = {-0.5, 0.5};
const Double_t mPtCut[2] = {0.2,50};
const Double_t mDcaCut = 1.;
const Double_t mHitsRatioCut = 0.52;
const Double_t mHitsFitCut = 20;
const Double_t mHitsDedxCut = 11;
const Double_t mnSigECut[2] = {-1.9, 3.};
const Double_t mnSigELooseCut[2] = {-2.5, 3.};
const Double_t myLocalCut = 2;
const Double_t mBetaCut[2] = {0.97,1.03};

const Double_t mEmcMinPt = 1.;
const Double_t mPveCut[4] = {0.3,1.5,0.2,2.2};
const Double_t mEmcHitsFitCut[5] = {20,19,21,18,22};
const Double_t mEmcHitsDedxCut = 11;
const Double_t mnEtaCut = 2;
const Double_t mnPhiCut = 2;
const Double_t mPhiDistCut = 0.01;
const Double_t mzDistCut = 2;
const Double_t mAdc0Cut[4] = {180,250,300,400};
const Double_t mEmcPtCut[4] = {2.5,3.6,4.3,5.7};
const Double_t adcfactor = 1.;
const Double_t mPairYCut[2] = {-1,1}; // 20160106
const Double_t mJpsiMassCut[2] = {3,3.2};
const Double_t mPairDcaCut = 3;

const int nSysCuts = 16;

const int mBadRuns[47]={
	//pp200 run12 VPDMB + BHT0,1,2 from Xiaozhi
	// http://portal.nersc.gov/project/star/xiaozhi/Run_QA/RunQA_V4.pdf
13041002,	
13041010,	
13042003,	
13044118,	
13045145,	
13047004,	
13047014,	
13047044,	
13047046,	
13047050,	
13048031,	
13049052,	
13052063,	
13053021,	
13054020,	
13054057,	
13057053,	
13058047,	
13058048,	
13059079,	
13060001,	
13061024,	
13061025,	
13061026,	
13064030,	
13066101,	
13066102,	
13066104,	
13066109,	
13067001,	
13068084,	
13069068,	
13070019,	
13070030,	
13070056,	
13070057,	
13070058,	
13070059,	
13070060,	
13070061,	
13071003,	
13071004,	
13071006,
13071006,	
13071008,	
13071034,
13071037
};
