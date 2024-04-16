//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 27 11:30:39 2024 by ROOT version 6.30/02
// from TChain Events/
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
//							//
//   File created by blehtela on 27.03.2024,		//
//   based on TChain of three 2023D files.		//
//   							//
//   NOTE: Removed some of the branches not used in 	//
//	   photon+jet analysis.				//
//	   e.g. electron, (boosted) tau, fatjet, mu...	//
//	   also remove L1 branches etc.			//
//	   only keep HLT_Photon* hlt branches		//
//							//
//////////////////////////////////////////////////////////


#ifndef PhotonJetAnalysis_h
#define PhotonJetAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>


//Add some stuff for JEC
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorUncertainty.h"

//some general includes (i used to do this in .C files, now moved to header)
#include <iostream>
#include <fstream> 
#include <cstdio>
#include <map>
#include <string>
using namespace std;


// Header file for the classes stored in the TTree if any.

class PhotonJetAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
			     
   //int	isMC; 		// from gamjet for mc handling --> might want to use another type of flag instead? think.
   bool		isMC;		// try with boolean variable
   string	dataset;	// for the data-taking period (like 2023D for example)
   string	version;	// for the version of the code


// Fixed size dimensions of array or collections stored in the TTree if any.
// originally copied these values from gamjet GamHistosFill.h, but might have changed over time.
   static const int nJetMax=200;
   static const int nPhotonMax=200;

   static const int nElectronMax=10;
   static const int nTauMax=10;
   static const int nMuonMax=20;



   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          bunchCrossing; //nano-v12?

   Float_t         CaloMET_phi;
   Float_t         CaloMET_pt;
   Float_t         CaloMET_sumEt;
   Float_t         ChsMET_phi;		//run2?
   Float_t         ChsMET_pt;		//run2?
   Float_t         ChsMET_sumEt;	//run2?

   Int_t           nCorrT1METJet;
   Float_t         CorrT1METJet_area[12];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_eta[12];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_muonSubtrFactor[12];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_phi[12];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_rawPt[12];   //[nCorrT1METJet]
   Float_t         DeepMETResolutionTune_phi;
   Float_t         DeepMETResolutionTune_pt;
   Float_t         DeepMETResponseTune_phi;
   Float_t         DeepMETResponseTune_pt;

   Int_t           nFsrPhoton;
   Short_t         FsrPhoton_electronIdx[3];   //[nFsrPhoton]
   Short_t         FsrPhoton_muonIdx[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_dROverEt2[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_eta[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_phi[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_pt[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_relIso03[3];   //[nFsrPhoton]
 
   Int_t           nJet;
   UChar_t         Jet_jetId[nJetMax];   //[nJet], size used to be 15
   UChar_t         Jet_nConstituents[nJetMax];   //[nJet]
   UChar_t         Jet_nElectrons[nJetMax];   //[nJet]
   UChar_t         Jet_nMuons[];   //[nJet]
   UChar_t         Jet_nSVs[];   //[nJet]
   Short_t         Jet_electronIdx1[];   //[nJet]
   Short_t         Jet_electronIdx2[];   //[nJet]
   Short_t         Jet_muonIdx1[];   //[nJet]
   Short_t         Jet_muonIdx2[];   //[nJet]
   Short_t         Jet_svIdx1[];   //[nJet]
   Short_t         Jet_svIdx2[];   //[nJet]
   Int_t           Jet_hfadjacentEtaStripsSize[];   //[nJet]
   Int_t           Jet_hfcentralEtaStripSize[];   //[nJet]
   Float_t         Jet_PNetRegPtRawCorr[];   //[nJet]
   Float_t         Jet_PNetRegPtRawCorrNeutrino[];   //[nJet]
   Float_t         Jet_PNetRegPtRawRes[];   //[nJet]
   Float_t         Jet_area[];   //[nJet]
   Float_t         Jet_btagDeepFlavB[];   //[nJet]
   Float_t         Jet_btagDeepFlavCvB[];   //[nJet]
   Float_t         Jet_btagDeepFlavCvL[];   //[nJet]
   Float_t         Jet_btagDeepFlavQG[];   //[nJet]
   Float_t         Jet_btagPNetB[];   //[nJet]
   Float_t         Jet_btagPNetCvB[];   //[nJet]
   Float_t         Jet_btagPNetCvL[];   //[nJet]
   Float_t         Jet_btagPNetQvG[];   //[nJet]
   Float_t         Jet_btagPNetTauVJet[];   //[nJet]
   Float_t         Jet_btagRobustParTAK4B[];   //[nJet]
   Float_t         Jet_btagRobustParTAK4CvB[];   //[nJet]
   Float_t         Jet_btagRobustParTAK4CvL[];   //[nJet]
   Float_t         Jet_btagRobustParTAK4QG[];   //[nJet]
   Float_t         Jet_chEmEF[];   //[nJet]
   Float_t         Jet_chHEF[];   //[nJet]
   Float_t         Jet_eta[];   //[nJet]
   Float_t         Jet_hfsigmaEtaEta[];   //[nJet]
   Float_t         Jet_hfsigmaPhiPhi[];   //[nJet]
   Float_t         Jet_mass[];   //[nJet]
   Float_t         Jet_muEF[];   //[nJet]
   Float_t         Jet_muonSubtrFactor[];   //[nJet]
   Float_t         Jet_neEmEF[];   //[nJet]
   Float_t         Jet_neHEF[];   //[nJet]
   Float_t         Jet_phi[];   //[nJet]
   Float_t         Jet_pt[];   //[nJet]
   Float_t         Jet_rawFactor[];   //[nJet]
				      
   Float_t         MET_MetUnclustEnUpDeltaX;
   Float_t         MET_MetUnclustEnUpDeltaY;
   Float_t         MET_covXX;
   Float_t         MET_covXY;
   Float_t         MET_covYY;
   Float_t         MET_phi;
   Float_t         MET_pt;
   Float_t         MET_significance;
   Float_t         MET_sumEt;
   Float_t         MET_sumPtUnclustered;

   Int_t           nPhoton;
   Char_t          Photon_seediEtaOriX[10];   //[nPhoton]
   UChar_t         Photon_cutBased[10];   //[nPhoton]
   Bool_t          Photon_electronVeto[10];   //[nPhoton]
   Bool_t          Photon_hasConversionTracks[10];   //[nPhoton]
   Bool_t          Photon_isScEtaEB[10];   //[nPhoton]
   Bool_t          Photon_isScEtaEE[10];   //[nPhoton]
   Bool_t          Photon_mvaID_WP80[10];   //[nPhoton]
   Bool_t          Photon_mvaID_WP90[10];   //[nPhoton]
   Bool_t          Photon_pixelSeed[10];   //[nPhoton]
   UChar_t         Photon_seedGain[10];   //[nPhoton]
   Short_t         Photon_electronIdx[10];   //[nPhoton]
   Short_t         Photon_jetIdx[10];   //[nPhoton]
   Int_t           Photon_seediPhiOriY[10];   //[nPhoton]
   Int_t           Photon_vidNestedWPBitmap[10];   //[nPhoton]
   Float_t         Photon_ecalPFClusterIso[10];   //[nPhoton]
   Float_t         Photon_energyErr[10];   //[nPhoton]
   Float_t         Photon_energyRaw[10];   //[nPhoton]
   Float_t         Photon_esEffSigmaRR[10];   //[nPhoton]
   Float_t         Photon_esEnergyOverRawE[10];   //[nPhoton]
   Float_t         Photon_eta[10];   //[nPhoton]
   Float_t         Photon_etaWidth[10];   //[nPhoton]
   Float_t         Photon_haloTaggerMVAVal[10];   //[nPhoton]
   Float_t         Photon_hcalPFClusterIso[10];   //[nPhoton]
   Float_t         Photon_hoe[10];   //[nPhoton]
   Float_t         Photon_hoe_PUcorr[10];   //[nPhoton]
   Float_t         Photon_mvaID[10];   //[nPhoton]
   Float_t         Photon_pfChargedIso[10];   //[nPhoton]
   Float_t         Photon_pfChargedIsoPFPV[10];   //[nPhoton]
   Float_t         Photon_pfChargedIsoWorstVtx[10];   //[nPhoton]
   Float_t         Photon_pfPhoIso03[10];   //[nPhoton]
   Float_t         Photon_pfRelIso03_all_quadratic[10];   //[nPhoton]
   Float_t         Photon_pfRelIso03_chg_quadratic[10];   //[nPhoton]
   Float_t         Photon_phi[10];   //[nPhoton]
   Float_t         Photon_phiWidth[10];   //[nPhoton]
   Float_t         Photon_pt[10];   //[nPhoton]
   Float_t         Photon_r9[10];   //[nPhoton]
   Float_t         Photon_s4[10];   //[nPhoton]
   Float_t         Photon_sieie[10];   //[nPhoton]
   Float_t         Photon_sieip[10];   //[nPhoton]
   Float_t         Photon_sipip[10];   //[nPhoton]
   Float_t         Photon_trkSumPtHollowConeDR03[10];   //[nPhoton]
   Float_t         Photon_trkSumPtSolidConeDR04[10];   //[nPhoton]
   Float_t         Photon_x_calo[10];   //[nPhoton]
   Float_t         Photon_y_calo[10];   //[nPhoton]
   Float_t         Photon_z_calo[10];   //[nPhoton]
					//
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_phiJERDown;
   Float_t         PuppiMET_phiJERUp;
   Float_t         PuppiMET_phiJESDown;
   Float_t         PuppiMET_phiJESUp;
   Float_t         PuppiMET_phiUnclusteredDown;
   Float_t         PuppiMET_phiUnclusteredUp;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_ptJERDown;
   Float_t         PuppiMET_ptJERUp;
   Float_t         PuppiMET_ptJESDown;
   Float_t         PuppiMET_ptJESUp;
   Float_t         PuppiMET_ptUnclusteredDown;
   Float_t         PuppiMET_ptUnclusteredUp;
   Float_t         PuppiMET_sumEt;
   Float_t         RawMET_phi;
   Float_t         RawMET_pt;
   Float_t         RawMET_sumEt;
   Float_t         RawPuppiMET_phi;	//run3?
   Float_t         RawPuppiMET_pt; 	//run3?
   Float_t         RawPuppiMET_sumEt;	//run3?
   Float_t         Rho_fixedGridRhoAll;
   Float_t         Rho_fixedGridRhoFastjetAll;
   Float_t         Rho_fixedGridRhoFastjetCentral;
   Float_t         Rho_fixedGridRhoFastjetCentralCalo;
   Float_t         Rho_fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         Rho_fixedGridRhoFastjetCentralNeutral;

   Int_t           nSoftActivityJet; //--> remove corresponding variables?
   Float_t         SoftActivityJet_eta[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_phi[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_pt[6];   //[nSoftActivityJet]
   Int_t           SoftActivityJetNjets10;
   Int_t           SoftActivityJetNjets2;
   Int_t           SoftActivityJetNjets5;
   Float_t         SoftActivityJetHT;
   Float_t         SoftActivityJetHT10;
   Float_t         SoftActivityJetHT2;
   Float_t         SoftActivityJetHT5;

   Int_t           nSubJet; //--> remove corresponding variables?
   Float_t         SubJet_btagDeepB[10];   //[nSubJet]
   Float_t         SubJet_eta[10];   //[nSubJet]
   Float_t         SubJet_mass[10];   //[nSubJet]
   Float_t         SubJet_n2b1[10];   //[nSubJet]
   Float_t         SubJet_n3b1[10];   //[nSubJet]
   Float_t         SubJet_phi[10];   //[nSubJet]
   Float_t         SubJet_pt[10];   //[nSubJet]
   Float_t         SubJet_rawFactor[10];   //[nSubJet]
   Float_t         SubJet_tau1[10];   //[nSubJet]
   Float_t         SubJet_tau2[10];   //[nSubJet]
   Float_t         SubJet_tau3[10];   //[nSubJet]
   Float_t         SubJet_tau4[10];   //[nSubJet]
				      //
   Float_t         TkMET_phi;
   Float_t         TkMET_pt;
   Float_t         TkMET_sumEt;

   Int_t           nTrigObj;
   Short_t         TrigObj_l1charge[67];   //[nTrigObj]
   UShort_t        TrigObj_id[67];   //[nTrigObj]
   Int_t           TrigObj_l1iso[67];   //[nTrigObj]
   Int_t           TrigObj_filterBits[67];   //[nTrigObj]
   Float_t         TrigObj_pt[67];   //[nTrigObj]
   Float_t         TrigObj_eta[67];   //[nTrigObj]
   Float_t         TrigObj_phi[67];   //[nTrigObj]
   Float_t         TrigObj_l1pt[67];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[67];   //[nTrigObj]
   Float_t         TrigObj_l2pt[67];   //[nTrigObj]
				       //
   Int_t           nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         OtherPV_score[3];   //[nOtherPV]
   UChar_t         PV_npvs;
   UChar_t         PV_npvsGood;
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           nSV;
   Short_t         SV_charge[15];   //[nSV]
   Float_t         SV_dlen[15];   //[nSV]
   Float_t         SV_dlenSig[15];   //[nSV]
   Float_t         SV_dxy[15];   //[nSV]
   Float_t         SV_dxySig[15];   //[nSV]
   Float_t         SV_pAngle[15];   //[nSV]
   UChar_t         SV_ntracks[15];   //[nSV]
   Float_t         SV_chi2[15];   //[nSV]
   Float_t         SV_eta[15];   //[nSV]
   Float_t         SV_mass[15];   //[nSV]
   Float_t         SV_ndof[15];   //[nSV]
   Float_t         SV_phi[15];   //[nSV]
   Float_t         SV_pt[15];   //[nSV]
   Float_t         SV_x[15];   //[nSV]
   Float_t         SV_y[15];   //[nSV]
   Float_t         SV_z[15];   //[nSV]
			       
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadPFMuonDzFilter;
   Bool_t          Flag_hfNoisyHitsFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;

   Bool_t          Flag_HBHENoiseFilter_pRECO;
   Bool_t          Flag_HBHENoiseIsoFilter_pRECO;
   Bool_t          Flag_CSCTightHaloFilter_pRECO;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO;
   Bool_t          Flag_CSCTightHalo2015Filter_pRECO;
   Bool_t          Flag_globalTightHalo2016Filter_pRECO;
   Bool_t          Flag_globalSuperTightHalo2016Filter_pRECO;
   Bool_t          Flag_HcalStripHaloFilter_pRECO;
   Bool_t          Flag_hcalLaserEventFilter_pRECO;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter_pRECO;
   Bool_t          Flag_ecalBadCalibFilter_pRECO;
   Bool_t          Flag_goodVertices_pRECO;
   Bool_t          Flag_eeBadScFilter_pRECO;
   Bool_t          Flag_ecalLaserCorrFilter_pRECO;
   Bool_t          Flag_trkPOGFilters_pRECO;
   Bool_t          Flag_chargedHadronTrackResolutionFilter_pRECO;
   Bool_t          Flag_muonBadTrackFilter_pRECO;
   Bool_t          Flag_BadChargedCandidateFilter_pRECO;
   Bool_t          Flag_BadPFMuonFilter_pRECO;
   Bool_t          Flag_BadPFMuonDzFilter_pRECO;
   Bool_t          Flag_hfNoisyHitsFilter_pRECO;
   Bool_t          Flag_BadChargedCandidateSummer16Filter_pRECO;
   Bool_t          Flag_BadPFMuonSummer16Filter_pRECO;
   Bool_t          Flag_trkPOG_manystripclus53X_pRECO;
   Bool_t          Flag_trkPOG_toomanystripclus53X_pRECO;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters_pRECO;
   Bool_t          Flag_METFilters_pRECO;

   Bool_t          HLTriggerFirstPath;
   Bool_t          HLT_EcalCalibration;
   Bool_t          HLT_HcalCalibration;
   Bool_t          HLT_HcalNZS;
   Bool_t          HLT_HcalPhiSym;
   Bool_t          HLT_Random;
   Bool_t          HLT_Physics;
   Bool_t          HLT_ZeroBias;
   Bool_t          HLT_ZeroBias_Alignment;
   Bool_t          HLT_ZeroBias_Beamspot;
   Bool_t          HLT_ZeroBias_IsolatedBunches;
   Bool_t          HLT_ZeroBias_FirstBXAfterTrain;
   Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap;
   Bool_t          HLT_ZeroBias_FirstCollisionInTrain;
   Bool_t          HLT_ZeroBias_LastCollisionInTrain;
   Bool_t          HLT_HT300_Beamspot;
   Bool_t          HLT_IsoTrackHB;
   Bool_t          HLT_IsoTrackHE;

   //comment these out for now, but could probably even remove from photon analysis! (TO DO)
   //Bool_t          HLT_PFJet40_GPUvsCPU;
   //Bool_t          HLT_AK8PFJet400_MassSD30;
   //Bool_t          HLT_AK8PFJet420_MassSD30;
   //Bool_t          HLT_AK8PFJet450_MassSD30;
   //Bool_t          HLT_AK8PFJet470_MassSD30;
   //Bool_t          HLT_AK8PFJet500_MassSD30;
   //Bool_t          HLT_CaloJet500_NoJetID;
   //Bool_t          HLT_CaloJet550_NoJetID;
   //Bool_t          HLT_DoubleEle25_CaloIdL_MW;
   //Bool_t          HLT_DoubleEle27_CaloIdL_MW;
   //Bool_t          HLT_DoubleEle33_CaloIdL_MW;
   //Bool_t          HLT_DoubleEle24_eta2p1_WPTight_Gsf;
   //Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
   //Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
   //Bool_t          HLT_Mu27_Ele37_CaloIdL_MW;
   //Bool_t          HLT_Mu37_Ele27_CaloIdL_MW;
   //Bool_t          HLT_Mu37_TkMu27;
   //Bool_t          HLT_DoubleMu4_3_Bs;
   //Bool_t          HLT_DoubleMu4_3_Jpsi;
   //Bool_t          HLT_DoubleMu4_3_LowMass;
   //Bool_t          HLT_DoubleMu4_LowMass_Displaced;
   //Bool_t          HLT_Mu0_L1DoubleMu;
   //Bool_t          HLT_Mu4_L1DoubleMu;
   //Bool_t          HLT_DoubleMu4_3_Photon4_BsToMMG;
   //Bool_t          HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG;
   //Bool_t          HLT_DoubleMu3_Trk_Tau3mu;
   //Bool_t          HLT_DoubleMu3_TkMu_DsTau3Mu;
   //Bool_t          HLT_DoubleMu4_Mass3p8_DZ_PFHT350;
   //Bool_t          HLT_DoubleMu4_MuMuTrk_Displaced;
   //Bool_t          HLT_Mu3_PFJet40;
   //Bool_t          HLT_Mu7p5_L2Mu2_Jpsi;
   //Bool_t          HLT_Mu7p5_L2Mu2_Upsilon;
   //Bool_t          HLT_Mu3_L1SingleMu5orSingleMu7;
   //Bool_t          HLT_DoublePhoton33_CaloIdL;
   //Bool_t          HLT_DoublePhoton70;
   //Bool_t          HLT_DoublePhoton85;
   //Bool_t          HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
   //Bool_t          HLT_Ele30_WPTight_Gsf;
   //Bool_t          HLT_Ele32_WPTight_Gsf;
   //Bool_t          HLT_Ele35_WPTight_Gsf;
   //Bool_t          HLT_Ele38_WPTight_Gsf;
   //Bool_t          HLT_Ele40_WPTight_Gsf;
   //Bool_t          HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   //Bool_t          HLT_IsoMu27_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;
   //Bool_t          HLT_IsoMu20;
   //Bool_t          HLT_IsoMu24;
   //Bool_t          HLT_IsoMu24_eta2p1;
   //Bool_t          HLT_IsoMu27;
   
   //Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   //Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   //Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   //Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   //Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   //Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   //Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   //Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   //Bool_t          HLT_Mu30_TkMu0_Psi;
   //Bool_t          HLT_Mu30_TkMu0_Upsilon;
   //Bool_t          HLT_Mu25_TkMu0_Phi;
   //Bool_t          HLT_Mu15;
   //Bool_t          HLT_Mu20;
   //Bool_t          HLT_Mu27;
   //Bool_t          HLT_Mu50;
   //Bool_t          HLT_Mu55;
   //Bool_t          HLT_CascadeMu100;
   //Bool_t          HLT_HighPtTkMu100;
   
   //Bool_t          HLT_AK8PFJet40;
   //Bool_t          HLT_AK8PFJet60;
   //Bool_t          HLT_AK8PFJet80;
   //Bool_t          HLT_AK8PFJet140;
   //Bool_t          HLT_AK8PFJet200;
   //Bool_t          HLT_AK8PFJet260;
   //Bool_t          HLT_AK8PFJet320;
   //Bool_t          HLT_AK8PFJet400;
   //Bool_t          HLT_AK8PFJet450;
   //Bool_t          HLT_AK8PFJet500;
   //Bool_t          HLT_AK8PFJet550;
   //Bool_t          HLT_PFJet40;
   //Bool_t          HLT_PFJet60;
   //Bool_t          HLT_PFJet80;
   //Bool_t          HLT_PFJet110;
   //Bool_t          HLT_PFJet140;
   //Bool_t          HLT_PFJet200;
   //Bool_t          HLT_PFJet260;
   //Bool_t          HLT_PFJet320;
   //Bool_t          HLT_PFJet400;
   //Bool_t          HLT_PFJet450;
   //Bool_t          HLT_PFJet500;
   //Bool_t          HLT_PFJet550;
   //Bool_t          HLT_PFJetFwd40;
   //Bool_t          HLT_PFJetFwd60;
   //Bool_t          HLT_PFJetFwd80;
   //Bool_t          HLT_PFJetFwd140;
   //Bool_t          HLT_PFJetFwd200;
   //Bool_t          HLT_PFJetFwd260;
   //Bool_t          HLT_PFJetFwd320;
   //Bool_t          HLT_PFJetFwd400;
   //Bool_t          HLT_PFJetFwd450;
   //Bool_t          HLT_PFJetFwd500;
   //Bool_t          HLT_AK8PFJetFwd15;
   //Bool_t          HLT_AK8PFJetFwd25;
   //Bool_t          HLT_AK8PFJetFwd40;
   //Bool_t          HLT_AK8PFJetFwd60;
   //Bool_t          HLT_AK8PFJetFwd80;
   //Bool_t          HLT_AK8PFJetFwd140;
   //Bool_t          HLT_AK8PFJetFwd200;
   //Bool_t          HLT_AK8PFJetFwd260;
   //Bool_t          HLT_AK8PFJetFwd320;
   //Bool_t          HLT_AK8PFJetFwd400;
   //Bool_t          HLT_AK8PFJetFwd450;
   //Bool_t          HLT_AK8PFJetFwd500;
   
   //Bool_t          HLT_PFHT180;
   //Bool_t          HLT_PFHT250;
   //Bool_t          HLT_PFHT370;
   //Bool_t          HLT_PFHT430;
   //Bool_t          HLT_PFHT510;
   //Bool_t          HLT_PFHT590;
   //Bool_t          HLT_PFHT680;
   //Bool_t          HLT_PFHT780;
   //Bool_t          HLT_PFHT890;
   //Bool_t          HLT_PFHT1050;
   //Bool_t          HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   //Bool_t          HLT_PFHT500_PFMET110_PFMHT110_IDTight;
   //Bool_t          HLT_PFHT700_PFMET85_PFMHT85_IDTight;
   //Bool_t          HLT_PFHT800_PFMET75_PFMHT75_IDTight;
   //Bool_t          HLT_PFMET120_PFMHT120_IDTight;
   //Bool_t          HLT_PFMET130_PFMHT130_IDTight;
   //Bool_t          HLT_PFMET140_PFMHT140_IDTight;
   //Bool_t          HLT_PFMET120_PFMHT120_IDTight_PFHT60;
   

   //TODO: continue to remove branches below this (of course only the ones that are not needed, keep photon stuff)
   
   //Bool_t          HLT_L1ETMHadSeeds;
   //Bool_t          HLT_CaloMHT90;
   //Bool_t          HLT_CaloMET90_NotCleaned;
   //Bool_t          HLT_CaloMET350_NotCleaned;
   //Bool_t          HLT_PFMET200_NotCleaned;
   //Bool_t          HLT_PFMET250_NotCleaned;
   //Bool_t          HLT_PFMET300_NotCleaned;
   //Bool_t          HLT_PFMET200_BeamHaloCleaned;
   //Bool_t          HLT_PFMETTypeOne200_BeamHaloCleaned;
   //Bool_t          HLT_MET105_IsoTrk50;
   //Bool_t          HLT_MET120_IsoTrk50;

   Bool_t          HLT_Photon300_NoHE;

   Bool_t          HLT_Photon33;
   Bool_t          HLT_Photon50;
   Bool_t          HLT_Photon75;
   Bool_t          HLT_Photon90;
   Bool_t          HLT_Photon120;
   Bool_t          HLT_Photon150;
   Bool_t          HLT_Photon175;
   Bool_t          HLT_Photon200;
   Bool_t          HLT_Photon30EB_TightID_TightIso;
   Bool_t          HLT_Photon50EB_TightID_TightIso;
   Bool_t          HLT_Photon75EB_TightID_TightIso;
   Bool_t          HLT_Photon90EB_TightID_TightIso;
   Bool_t          HLT_Photon110EB_TightID_TightIso;
   Bool_t          HLT_Photon130EB_TightID_TightIso;
   Bool_t          HLT_Photon150EB_TightID_TightIso;
   Bool_t          HLT_Photon175EB_TightID_TightIso;
   Bool_t          HLT_Photon200EB_TightID_TightIso;
   Bool_t          HLT_Photon100EBHE10;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon90_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon120_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon165_R9Id90_HE10_IsoM;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;
   Bool_t          HLT_Photon35_TwoProngs35;


   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400;

   Bool_t          HLT_ECALHT800;

   Bool_t          HLT_DiSC30_18_EIso_AND_HE_Mass70;
   Bool_t          HLT_Photon20_HoverELoose;
   Bool_t          HLT_Photon30_HoverELoose;
   
   Bool_t          HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;

   Bool_t          HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;
   Bool_t          HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;


   //Bool_t          HLT_HT350;
   //Bool_t          HLT_HT425;
   //Bool_t          HLT_CaloMET60_DTCluster50;
   //Bool_t          HLT_CaloMET60_DTClusterNoMB1S50;
   //Bool_t          HLT_L1MET_DTCluster50;
   //Bool_t          HLT_L1MET_DTClusterNoMB1S50;
   //Bool_t          HLT_CscCluster_Loose;
   //Bool_t          HLT_CscCluster_Medium;
   //Bool_t          HLT_CscCluster_Tight;
   //Bool_t          HLT_DoubleCscCluster75;
   //Bool_t          HLT_DoubleCscCluster100;
   //Bool_t          HLT_L1CSCShower_DTCluster50;
   //Bool_t          HLT_L1CSCShower_DTCluster75;
   //Bool_t          HLT_PFMET105_IsoTrk50;

   /*
   Bool_t          HLT_DiPhoton10Time1ns;
   Bool_t          HLT_DiPhoton10Time1p2ns;
   Bool_t          HLT_DiPhoton10Time1p4ns;
   Bool_t          HLT_DiPhoton10Time1p6ns;
   Bool_t          HLT_DiPhoton10Time1p8ns;
   Bool_t          HLT_DiPhoton10Time2ns;
   Bool_t          HLT_DiPhoton10sminlt0p1;
   Bool_t          HLT_DiPhoton10sminlt0p12;
   Bool_t          HLT_DiPhoton10_CaloIdL;
   */

   /*
   Bool_t          HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   */

   //Bool_t          HLT_Photon32_OneProng32_M50To105;
   //Bool_t          HLT_Photon50_TimeLtNeg2p5ns;
   //Bool_t          HLT_Photon50_TimeGt2p5ns;
   //Bool_t          HLTriggerFinalPath;



   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_BeamSpot_type;   //!
   TBranch        *b_BeamSpot_sigmaZ;   //!
   TBranch        *b_BeamSpot_sigmaZError;   //!
   TBranch        *b_BeamSpot_z;   //!
   TBranch        *b_BeamSpot_zError;   //!
					
   TBranch        *b_CaloMET_phi;   //!
   TBranch        *b_CaloMET_pt;   //!
   TBranch        *b_CaloMET_sumEt;   //!
   TBranch        *b_ChsMET_phi;   //!
   TBranch        *b_ChsMET_pt;   //!
   TBranch        *b_ChsMET_sumEt;   //!
   TBranch        *b_nCorrT1METJet;   //!
   TBranch        *b_CorrT1METJet_area;   //!
   TBranch        *b_CorrT1METJet_eta;   //!
   TBranch        *b_CorrT1METJet_muonSubtrFactor;   //!
   TBranch        *b_CorrT1METJet_phi;   //!
   TBranch        *b_CorrT1METJet_rawPt;   //!
   TBranch        *b_DeepMETResolutionTune_phi;   //!
   TBranch        *b_DeepMETResolutionTune_pt;   //!
   TBranch        *b_DeepMETResponseTune_phi;   //!
   TBranch        *b_DeepMETResponseTune_pt;   //!
					       
   TBranch        *b_nFsrPhoton;   //!
   TBranch        *b_FsrPhoton_electronIdx;   //!
   TBranch        *b_FsrPhoton_muonIdx;   //!
   TBranch        *b_FsrPhoton_dROverEt2;   //!
   TBranch        *b_FsrPhoton_eta;   //!
   TBranch        *b_FsrPhoton_phi;   //!
   TBranch        *b_FsrPhoton_pt;   //!
   TBranch        *b_FsrPhoton_relIso03;   //!
					   
						  
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
   TBranch        *b_Jet_nSVs;   //!
   TBranch        *b_Jet_electronIdx1;   //!
   TBranch        *b_Jet_electronIdx2;   //!
   TBranch        *b_Jet_muonIdx1;   //!
   TBranch        *b_Jet_muonIdx2;   //!
   TBranch        *b_Jet_svIdx1;   //!
   TBranch        *b_Jet_svIdx2;   //!
   TBranch        *b_Jet_hfadjacentEtaStripsSize;   //!
   TBranch        *b_Jet_hfcentralEtaStripSize;   //!
   TBranch        *b_Jet_PNetRegPtRawCorr;   //!
   TBranch        *b_Jet_PNetRegPtRawCorrNeutrino;   //!
   TBranch        *b_Jet_PNetRegPtRawRes;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
   TBranch        *b_Jet_btagDeepFlavCvB;   //!
   TBranch        *b_Jet_btagDeepFlavCvL;   //!
   TBranch        *b_Jet_btagDeepFlavQG;   //!
   TBranch        *b_Jet_btagPNetB;   //!
   TBranch        *b_Jet_btagPNetCvB;   //!
   TBranch        *b_Jet_btagPNetCvL;   //!
   TBranch        *b_Jet_btagPNetQvG;   //!
   TBranch        *b_Jet_btagPNetTauVJet;   //!
   TBranch        *b_Jet_btagRobustParTAK4B;   //!
   TBranch        *b_Jet_btagRobustParTAK4CvB;   //!
   TBranch        *b_Jet_btagRobustParTAK4CvL;   //!
   TBranch        *b_Jet_btagRobustParTAK4QG;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_hfsigmaEtaEta;   //!
   TBranch        *b_Jet_hfsigmaPhiPhi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_muonSubtrFactor;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_rawFactor;   //!
				      

   TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_MET_covXX;   //!
   TBranch        *b_MET_covXY;   //!
   TBranch        *b_MET_covYY;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_MET_sumPtUnclustered;   //!
					  
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_seediEtaOriX;   //!
   TBranch        *b_Photon_cutBased;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_hasConversionTracks;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_mvaID_WP80;   //!
   TBranch        *b_Photon_mvaID_WP90;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_seedGain;   //!
   TBranch        *b_Photon_electronIdx;   //!
   TBranch        *b_Photon_jetIdx;   //!
   TBranch        *b_Photon_seediPhiOriY;   //!
   TBranch        *b_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_Photon_ecalPFClusterIso;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_energyRaw;   //!
   TBranch        *b_Photon_esEffSigmaRR;   //!
   TBranch        *b_Photon_esEnergyOverRawE;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_etaWidth;   //!
   TBranch        *b_Photon_haloTaggerMVAVal;   //!
   TBranch        *b_Photon_hcalPFClusterIso;   //!
   TBranch        *b_Photon_hoe;   //!
   TBranch        *b_Photon_hoe_PUcorr;   //!
   TBranch        *b_Photon_mvaID;   //!
   TBranch        *b_Photon_pfChargedIso;   //!
   TBranch        *b_Photon_pfChargedIsoPFPV;   //!
   TBranch        *b_Photon_pfChargedIsoWorstVtx;   //!
   TBranch        *b_Photon_pfPhoIso03;   //!
   TBranch        *b_Photon_pfRelIso03_all_quadratic;   //!
   TBranch        *b_Photon_pfRelIso03_chg_quadratic;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_phiWidth;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_s4;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_sieip;   //!
   TBranch        *b_Photon_sipip;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_x_calo;   //!
   TBranch        *b_Photon_y_calo;   //!
   TBranch        *b_Photon_z_calo;   //!
				      
   TBranch        *b_nPPSLocalTrack;   //!
   TBranch        *b_PPSLocalTrack_multiRPProtonIdx;   //!
   TBranch        *b_PPSLocalTrack_singleRPProtonIdx;   //!
   TBranch        *b_PPSLocalTrack_decRPId;   //!
   TBranch        *b_PPSLocalTrack_rpType;   //!
   TBranch        *b_PPSLocalTrack_x;   //!
   TBranch        *b_PPSLocalTrack_y;   //!
   TBranch        *b_PPSLocalTrack_time;   //!
   TBranch        *b_PPSLocalTrack_timeUnc;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_phiJERDown;   //!
   TBranch        *b_PuppiMET_phiJERUp;   //!
   TBranch        *b_PuppiMET_phiJESDown;   //!
   TBranch        *b_PuppiMET_phiJESUp;   //!
   TBranch        *b_PuppiMET_phiUnclusteredDown;   //!
   TBranch        *b_PuppiMET_phiUnclusteredUp;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_ptJERDown;   //!
   TBranch        *b_PuppiMET_ptJERUp;   //!
   TBranch        *b_PuppiMET_ptJESDown;   //!
   TBranch        *b_PuppiMET_ptJESUp;   //!
   TBranch        *b_PuppiMET_ptUnclusteredDown;   //!
   TBranch        *b_PuppiMET_ptUnclusteredUp;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_RawMET_phi;   //!
   TBranch        *b_RawMET_pt;   //!
   TBranch        *b_RawMET_sumEt;   //!
   TBranch        *b_RawPuppiMET_phi;   //!
   TBranch        *b_RawPuppiMET_pt;   //!
   TBranch        *b_RawPuppiMET_sumEt;   //!
   TBranch        *b_Rho_fixedGridRhoAll;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetAll;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralNeutral;   //!
							      
   //TBranch        *b_nSoftActivityJet;   //!

   TBranch        *b_TkMET_phi;   //!
   TBranch        *b_TkMET_pt;   //!
   TBranch        *b_TkMET_sumEt;   //!
				    
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
				     //
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_OtherPV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_charge;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_dxy;   //!
   TBranch        *b_SV_dxySig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_SV_ntracks;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!

   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!
					
   //TBranch        *b_L1_SingleEG10er2p5;   //! --> are those of interest to us? removed all similar ones...
   TBranch        *b_Flag_HBHENoiseFilter_pRECO;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter_pRECO;   //!
   TBranch        *b_Flag_CSCTightHaloFilter_pRECO;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter_pRECO;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter_pRECO;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter_pRECO;   //!
   TBranch        *b_Flag_HcalStripHaloFilter_pRECO;   //!
   TBranch        *b_Flag_hcalLaserEventFilter_pRECO;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter_pRECO;   //!
   TBranch        *b_Flag_ecalBadCalibFilter_pRECO;   //!
   TBranch        *b_Flag_goodVertices_pRECO;   //!
   TBranch        *b_Flag_eeBadScFilter_pRECO;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter_pRECO;   //!
   TBranch        *b_Flag_trkPOGFilters_pRECO;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter_pRECO;   //!
   TBranch        *b_Flag_muonBadTrackFilter_pRECO;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter_pRECO;   //!
   TBranch        *b_Flag_BadPFMuonFilter_pRECO;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter_pRECO;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter_pRECO;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter_pRECO;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter_pRECO;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X_pRECO;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X_pRECO;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters_pRECO;   //!
   TBranch        *b_Flag_METFilters_pRECO;   //!
   TBranch        *b_L1_AlwaysTrue_pRECO;   //!
   TBranch        *b_L1_MinimumBiasHF0_pRECO;   //!
   TBranch        *b_L1_MinimumBiasHF0_AND_BptxAND_pRECO;   //!
   TBranch        *b_L1_NotBptxOR_pRECO;   //!
   TBranch        *b_L1_ZeroBias_pRECO;   //!
   TBranch        *b_L1_ZeroBias_copy_pRECO;   //!
   TBranch        *b_L1_UnprefireableEvent_pRECO;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLT_EphemeralPhysics;   //!
   TBranch        *b_HLT_EphemeralZeroBias;   //!
   TBranch        *b_HLT_EcalCalibration;   //!
   TBranch        *b_HLT_HcalCalibration;   //!
   TBranch        *b_HLT_HcalNZS;   //!
   TBranch        *b_HLT_HcalPhiSym;   //!
   TBranch        *b_HLT_Random;   //!
   TBranch        *b_HLT_Physics;   //!
   TBranch        *b_HLT_ZeroBias;   //!
   TBranch        *b_HLT_ZeroBias_Alignment;   //!
   TBranch        *b_HLT_ZeroBias_Beamspot;   //!
   TBranch        *b_HLT_ZeroBias_IsolatedBunches;   //!
   TBranch        *b_HLT_ZeroBias_FirstBXAfterTrain;   //!
   TBranch        *b_HLT_ZeroBias_FirstCollisionAfterAbortGap;   //!
   TBranch        *b_HLT_ZeroBias_FirstCollisionInTrain;   //!
   TBranch        *b_HLT_ZeroBias_LastCollisionInTrain;   //!
   TBranch        *b_HLT_HT300_Beamspot;   //!
   TBranch        *b_HLT_IsoTrackHB;   //!
   TBranch        *b_HLT_IsoTrackHE;   //!
				       
   //TBranch        *b_HLT_CaloJet500_NoJetID;   //!
   //TBranch        *b_HLT_CaloJet550_NoJetID;   //!
   //TBranch        *b_HLT_DoublePhoton33_CaloIdL;   //!
   //TBranch        *b_HLT_DoublePhoton70;   //!
   //TBranch        *b_HLT_DoublePhoton85;   //!
   
   TBranch        *b_HLT_UncorrectedJetE30_NoBPTX;   //!
   TBranch        *b_HLT_UncorrectedJetE30_NoBPTX3BX;   //!
   TBranch        *b_HLT_UncorrectedJetE60_NoBPTX3BX;   //!
   TBranch        *b_HLT_UncorrectedJetE70_NoBPTX3BX;   //!
							
   /*
   TBranch        *b_HLT_AK8PFJet40;   //!
   TBranch        *b_HLT_AK8PFJet60;   //!
   TBranch        *b_HLT_AK8PFJet80;   //!
   TBranch        *b_HLT_AK8PFJet140;   //!
   TBranch        *b_HLT_AK8PFJet200;   //!
   TBranch        *b_HLT_AK8PFJet260;   //!
   TBranch        *b_HLT_AK8PFJet320;   //!
   TBranch        *b_HLT_AK8PFJet400;   //!
   TBranch        *b_HLT_AK8PFJet450;   //!
   TBranch        *b_HLT_AK8PFJet500;   //!
   TBranch        *b_HLT_AK8PFJet550;   //!
   TBranch        *b_HLT_PFJet40;   //!
   TBranch        *b_HLT_PFJet60;   //!
   TBranch        *b_HLT_PFJet80;   //!
   TBranch        *b_HLT_PFJet110;   //!
   TBranch        *b_HLT_PFJet140;   //!
   TBranch        *b_HLT_PFJet200;   //!
   TBranch        *b_HLT_PFJet260;   //!
   TBranch        *b_HLT_PFJet320;   //!
   TBranch        *b_HLT_PFJet400;   //!
   TBranch        *b_HLT_PFJet450;   //!
   TBranch        *b_HLT_PFJet500;   //!
   TBranch        *b_HLT_PFJet550;   //!
   TBranch        *b_HLT_PFJetFwd40;   //!
   TBranch        *b_HLT_PFJetFwd60;   //!
   TBranch        *b_HLT_PFJetFwd80;   //!
   TBranch        *b_HLT_PFJetFwd140;   //!
   TBranch        *b_HLT_PFJetFwd200;   //!
   TBranch        *b_HLT_PFJetFwd260;   //!
   TBranch        *b_HLT_PFJetFwd320;   //!
   TBranch        *b_HLT_PFJetFwd400;   //!
   TBranch        *b_HLT_PFJetFwd450;   //!
   TBranch        *b_HLT_PFJetFwd500;   //!
   TBranch        *b_HLT_AK8PFJetFwd15;   //!
   TBranch        *b_HLT_AK8PFJetFwd25;   //!
   TBranch        *b_HLT_AK8PFJetFwd40;   //!
   TBranch        *b_HLT_AK8PFJetFwd60;   //!
   TBranch        *b_HLT_AK8PFJetFwd80;   //!
   TBranch        *b_HLT_AK8PFJetFwd140;   //!
   TBranch        *b_HLT_AK8PFJetFwd200;   //!
   TBranch        *b_HLT_AK8PFJetFwd260;   //!
   TBranch        *b_HLT_AK8PFJetFwd320;   //!
   TBranch        *b_HLT_AK8PFJetFwd400;   //!
   TBranch        *b_HLT_AK8PFJetFwd450;   //!
   TBranch        *b_HLT_AK8PFJetFwd500;   //!
   TBranch        *b_HLT_PFHT180;   //!
   TBranch        *b_HLT_PFHT250;   //!
   TBranch        *b_HLT_PFHT370;   //!
   TBranch        *b_HLT_PFHT430;   //!
   TBranch        *b_HLT_PFHT510;   //!
   TBranch        *b_HLT_PFHT590;   //!
   TBranch        *b_HLT_PFHT680;   //!
   TBranch        *b_HLT_PFHT780;   //!
   TBranch        *b_HLT_PFHT890;   //!
   TBranch        *b_HLT_PFHT1050;   //!
   */

   TBranch        *b_HLT_PFHT500_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_HLT_PFHT500_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFHT700_PFMET85_PFMHT85_IDTight;   //!
   TBranch        *b_HLT_PFHT800_PFMET75_PFMHT75_IDTight;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_PFMET130_PFMHT130_IDTight;   //!
   TBranch        *b_HLT_PFMET140_PFMHT140_IDTight;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETTypeOne140_PFMHT140_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;   //!
   TBranch        *b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF;   //!
   TBranch        *b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF;   //!
   TBranch        *b_HLT_L1ETMHadSeeds;   //!
   TBranch        *b_HLT_CaloMHT90;   //!
   TBranch        *b_HLT_CaloMET90_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET350_NotCleaned;   //!
   TBranch        *b_HLT_PFMET200_NotCleaned;   //!
   TBranch        *b_HLT_PFMET250_NotCleaned;   //!
   TBranch        *b_HLT_PFMET300_NotCleaned;   //!
   TBranch        *b_HLT_PFMET200_BeamHaloCleaned;   //!
   TBranch        *b_HLT_PFMETTypeOne200_BeamHaloCleaned;   //!
   TBranch        *b_HLT_MET105_IsoTrk50;   //!
   TBranch        *b_HLT_MET120_IsoTrk50;   //!
					       //
   TBranch        *b_HLT_DoublePFJets40_PFBTagDeepJet_p71;   //!
   TBranch        *b_HLT_DoublePFJets100_PFBTagDeepJet_p71;   //!
   TBranch        *b_HLT_DoublePFJets200_PFBTagDeepJet_p71;   //!
   TBranch        *b_HLT_DoublePFJets350_PFBTagDeepJet_p71;   //!
   TBranch        *b_HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71;   //!
   TBranch        *b_HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71;   //!
									      
   TBranch        *b_HLT_Photon300_NoHE;   //!

   TBranch        *b_HLT_Photon33;   //!
   TBranch        *b_HLT_Photon50;   //!
   TBranch        *b_HLT_Photon75;   //!
   TBranch        *b_HLT_Photon90;   //!
   TBranch        *b_HLT_Photon120;   //!
   TBranch        *b_HLT_Photon150;   //!
   TBranch        *b_HLT_Photon175;   //!
   TBranch        *b_HLT_Photon200;   //!
   TBranch        *b_HLT_Photon30EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon50EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon75EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon90EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon110EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon130EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon150EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon175EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon200EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon100EBHE10;   //!
   TBranch        *b_HLT_Photon50_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon90_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon120_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon165_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;   //!
   TBranch        *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;   //!
   TBranch        *b_HLT_Photon35_TwoProngs35;   //!
						 
   //TBranch        *b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40;   //!
   //TBranch        *b_HLT_PFHT400_SixPFJet32;   //!
   //TBranch        *b_HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50;   //!
   //TBranch        *b_HLT_PFHT450_SixPFJet36;   //!
   //TBranch        *b_HLT_PFHT450_SixPFJet36_PNetBTag0p35;   //!
   //TBranch        *b_HLT_PFHT400_FivePFJet_100_100_60_30_30;   //!
   //TBranch        *b_HLT_PFHT350;   //!
				    
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400;   //!
									      
   //TBranch        *b_HLT_ECALHT800;   //!
   //TBranch        *b_HLT_DiSC30_18_EIso_AND_HE_Mass70;   //!
   TBranch        *b_HLT_Photon20_HoverELoose;   //!
   TBranch        *b_HLT_Photon30_HoverELoose;   //!
						 
   TBranch        *b_HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;   //!
									       
   TBranch        *b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;   //!
   TBranch        *b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;   //!
									    
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05;   //!
   TBranch        *b_HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10;   //!
   TBranch        *b_HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03;   //!
   TBranch        *b_HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05;   //!
   TBranch        *b_HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10;   //!
   TBranch        *b_HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03;   //!
   TBranch        *b_HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05;   //!
   TBranch        *b_HLT_AK8PFJet425_SoftDropMass40;   //!
   TBranch        *b_HLT_AK8PFJet450_SoftDropMass40;   //!
						       //
   //TBranch        *b_HLT_HT350;   //!
   //TBranch        *b_HLT_HT425;   //!
   //TBranch        *b_HLT_CaloMET60_DTCluster50;   //!
   //TBranch        *b_HLT_CaloMET60_DTClusterNoMB1S50;   //!
   //TBranch        *b_HLT_L1MET_DTCluster50;   //!
   //TBranch        *b_HLT_L1MET_DTClusterNoMB1S50;   //!
   TBranch        *b_HLT_CscCluster_Loose;   //!
   TBranch        *b_HLT_CscCluster_Medium;   //!
   TBranch        *b_HLT_CscCluster_Tight;   //!
   TBranch        *b_HLT_DoubleCscCluster75;   //!
   TBranch        *b_HLT_DoubleCscCluster100;   //!
   TBranch        *b_HLT_L1CSCShower_DTCluster50;   //!
   TBranch        *b_HLT_L1CSCShower_DTCluster75;   //!
   TBranch        *b_HLT_PFMET105_IsoTrk50;   //!
   TBranch        *b_HLT_L1SingleLLPJet;   //!

										      
   TBranch        *b_HLT_DiPhoton10Time1ns;   //!
   TBranch        *b_HLT_DiPhoton10Time1p2ns;   //!
   TBranch        *b_HLT_DiPhoton10Time1p4ns;   //!
   TBranch        *b_HLT_DiPhoton10Time1p6ns;   //!
   TBranch        *b_HLT_DiPhoton10Time1p8ns;   //!
   TBranch        *b_HLT_DiPhoton10Time2ns;   //!
   TBranch        *b_HLT_DiPhoton10sminlt0p1;   //!
   TBranch        *b_HLT_DiPhoton10sminlt0p12;   //!
   TBranch        *b_HLT_DiPhoton10_CaloIdL;   //!
					       
   TBranch        *b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;   //!
									       
   TBranch        *b_HLT_Photon32_OneProng32_M50To105;   //!
							 
   /*
   TBranch        *b_HLT_PFJet200_TimeLtNeg2p5ns;   //!
   TBranch        *b_HLT_PFJet200_TimeGt2p5ns;   //!
   TBranch        *b_HLT_Photon50_TimeLtNeg2p5ns;   //!
   TBranch        *b_HLT_Photon50_TimeGt2p5ns;   //!
   TBranch        *b_HLT_ExpressMuons;   //!
   TBranch        *b_HLT_PPSMaxTracksPerArm1;   //!
   TBranch        *b_HLT_PPSMaxTracksPerRP4;   //!
   TBranch        *b_HLT_PPSRandom;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   */



   PhotonJetAnalysis(TTree *tree=0);
   virtual ~PhotonJetAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PhotonJetAnalysis_cxx
PhotonJetAnalysis::PhotonJetAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("Events",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("Events","");
      chain->Add("data/2ac63590-ec9d-4d95-a37e-13f8819e0503.root/Events");
      chain->Add("data/0c71c68e-416e-42a7-a6cf-d7d168dad678.root/Events");
      chain->Add("data/18c12337-03a0-4e82-a10a-307ece47505a.root/Events");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

PhotonJetAnalysis::~PhotonJetAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PhotonJetAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PhotonJetAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PhotonJetAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("BeamSpot_type", &BeamSpot_type, &b_BeamSpot_type);
   fChain->SetBranchAddress("BeamSpot_sigmaZ", &BeamSpot_sigmaZ, &b_BeamSpot_sigmaZ);
   fChain->SetBranchAddress("BeamSpot_sigmaZError", &BeamSpot_sigmaZError, &b_BeamSpot_sigmaZError);
   fChain->SetBranchAddress("BeamSpot_z", &BeamSpot_z, &b_BeamSpot_z);
   fChain->SetBranchAddress("BeamSpot_zError", &BeamSpot_zError, &b_BeamSpot_zError);

   fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   fChain->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi);
   fChain->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt);
   fChain->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt);
   fChain->SetBranchAddress("nCorrT1METJet", &nCorrT1METJet, &b_nCorrT1METJet);
   fChain->SetBranchAddress("CorrT1METJet_area", CorrT1METJet_area, &b_CorrT1METJet_area);
   fChain->SetBranchAddress("CorrT1METJet_eta", CorrT1METJet_eta, &b_CorrT1METJet_eta);
   fChain->SetBranchAddress("CorrT1METJet_muonSubtrFactor", CorrT1METJet_muonSubtrFactor, &b_CorrT1METJet_muonSubtrFactor);
   fChain->SetBranchAddress("CorrT1METJet_phi", CorrT1METJet_phi, &b_CorrT1METJet_phi);
   fChain->SetBranchAddress("CorrT1METJet_rawPt", CorrT1METJet_rawPt, &b_CorrT1METJet_rawPt);
   fChain->SetBranchAddress("DeepMETResolutionTune_phi", &DeepMETResolutionTune_phi, &b_DeepMETResolutionTune_phi);
   fChain->SetBranchAddress("DeepMETResolutionTune_pt", &DeepMETResolutionTune_pt, &b_DeepMETResolutionTune_pt);
   fChain->SetBranchAddress("DeepMETResponseTune_phi", &DeepMETResponseTune_phi, &b_DeepMETResponseTune_phi);
   fChain->SetBranchAddress("DeepMETResponseTune_pt", &DeepMETResponseTune_pt, &b_DeepMETResponseTune_pt);

   fChain->SetBranchAddress("nFsrPhoton", &nFsrPhoton, &b_nFsrPhoton);
   fChain->SetBranchAddress("FsrPhoton_electronIdx", FsrPhoton_electronIdx, &b_FsrPhoton_electronIdx);
   fChain->SetBranchAddress("FsrPhoton_muonIdx", FsrPhoton_muonIdx, &b_FsrPhoton_muonIdx);
   fChain->SetBranchAddress("FsrPhoton_dROverEt2", FsrPhoton_dROverEt2, &b_FsrPhoton_dROverEt2);
   fChain->SetBranchAddress("FsrPhoton_eta", FsrPhoton_eta, &b_FsrPhoton_eta);
   fChain->SetBranchAddress("FsrPhoton_phi", FsrPhoton_phi, &b_FsrPhoton_phi);
   fChain->SetBranchAddress("FsrPhoton_pt", FsrPhoton_pt, &b_FsrPhoton_pt);
   fChain->SetBranchAddress("FsrPhoton_relIso03", FsrPhoton_relIso03, &b_FsrPhoton_relIso03);

   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_nSVs", Jet_nSVs, &b_Jet_nSVs);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_svIdx1", Jet_svIdx1, &b_Jet_svIdx1);
   fChain->SetBranchAddress("Jet_svIdx2", Jet_svIdx2, &b_Jet_svIdx2);
   fChain->SetBranchAddress("Jet_hfadjacentEtaStripsSize", Jet_hfadjacentEtaStripsSize, &b_Jet_hfadjacentEtaStripsSize);
   fChain->SetBranchAddress("Jet_hfcentralEtaStripSize", Jet_hfcentralEtaStripSize, &b_Jet_hfcentralEtaStripSize);
   fChain->SetBranchAddress("Jet_PNetRegPtRawCorr", Jet_PNetRegPtRawCorr, &b_Jet_PNetRegPtRawCorr);
   fChain->SetBranchAddress("Jet_PNetRegPtRawCorrNeutrino", Jet_PNetRegPtRawCorrNeutrino, &b_Jet_PNetRegPtRawCorrNeutrino);
   fChain->SetBranchAddress("Jet_PNetRegPtRawRes", Jet_PNetRegPtRawRes, &b_Jet_PNetRegPtRawRes);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvB", Jet_btagDeepFlavCvB, &b_Jet_btagDeepFlavCvB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvL", Jet_btagDeepFlavCvL, &b_Jet_btagDeepFlavCvL);
   fChain->SetBranchAddress("Jet_btagDeepFlavQG", Jet_btagDeepFlavQG, &b_Jet_btagDeepFlavQG);
   fChain->SetBranchAddress("Jet_btagPNetB", Jet_btagPNetB, &b_Jet_btagPNetB);
   fChain->SetBranchAddress("Jet_btagPNetCvB", Jet_btagPNetCvB, &b_Jet_btagPNetCvB);
   fChain->SetBranchAddress("Jet_btagPNetCvL", Jet_btagPNetCvL, &b_Jet_btagPNetCvL);
   fChain->SetBranchAddress("Jet_btagPNetQvG", Jet_btagPNetQvG, &b_Jet_btagPNetQvG);
   fChain->SetBranchAddress("Jet_btagPNetTauVJet", Jet_btagPNetTauVJet, &b_Jet_btagPNetTauVJet);
   fChain->SetBranchAddress("Jet_btagRobustParTAK4B", Jet_btagRobustParTAK4B, &b_Jet_btagRobustParTAK4B);
   fChain->SetBranchAddress("Jet_btagRobustParTAK4CvB", Jet_btagRobustParTAK4CvB, &b_Jet_btagRobustParTAK4CvB);
   fChain->SetBranchAddress("Jet_btagRobustParTAK4CvL", Jet_btagRobustParTAK4CvL, &b_Jet_btagRobustParTAK4CvL);
   fChain->SetBranchAddress("Jet_btagRobustParTAK4QG", Jet_btagRobustParTAK4QG, &b_Jet_btagRobustParTAK4QG);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_hfsigmaEtaEta", Jet_hfsigmaEtaEta, &b_Jet_hfsigmaEtaEta);
   fChain->SetBranchAddress("Jet_hfsigmaPhiPhi", Jet_hfsigmaPhiPhi, &b_Jet_hfsigmaPhiPhi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);

   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
   fChain->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
   fChain->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("MET_sumPtUnclustered", &MET_sumPtUnclustered, &b_MET_sumPtUnclustered);

   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_seediEtaOriX", Photon_seediEtaOriX, &b_Photon_seediEtaOriX);
   fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased);
   fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_hasConversionTracks", Photon_hasConversionTracks, &b_Photon_hasConversionTracks);
   fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);
   fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   fChain->SetBranchAddress("Photon_seediPhiOriY", Photon_seediPhiOriY, &b_Photon_seediPhiOriY);
   fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("Photon_ecalPFClusterIso", Photon_ecalPFClusterIso, &b_Photon_ecalPFClusterIso);
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_energyRaw", Photon_energyRaw, &b_Photon_energyRaw);
   fChain->SetBranchAddress("Photon_esEffSigmaRR", Photon_esEffSigmaRR, &b_Photon_esEffSigmaRR);
   fChain->SetBranchAddress("Photon_esEnergyOverRawE", Photon_esEnergyOverRawE, &b_Photon_esEnergyOverRawE);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth, &b_Photon_etaWidth);
   fChain->SetBranchAddress("Photon_haloTaggerMVAVal", Photon_haloTaggerMVAVal, &b_Photon_haloTaggerMVAVal);
   fChain->SetBranchAddress("Photon_hcalPFClusterIso", Photon_hcalPFClusterIso, &b_Photon_hcalPFClusterIso);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   fChain->SetBranchAddress("Photon_hoe_PUcorr", Photon_hoe_PUcorr, &b_Photon_hoe_PUcorr);
   fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   fChain->SetBranchAddress("Photon_pfChargedIso", Photon_pfChargedIso, &b_Photon_pfChargedIso);
   fChain->SetBranchAddress("Photon_pfChargedIsoPFPV", Photon_pfChargedIsoPFPV, &b_Photon_pfChargedIsoPFPV);
   fChain->SetBranchAddress("Photon_pfChargedIsoWorstVtx", Photon_pfChargedIsoWorstVtx, &b_Photon_pfChargedIsoWorstVtx);
   fChain->SetBranchAddress("Photon_pfPhoIso03", Photon_pfPhoIso03, &b_Photon_pfPhoIso03);
   fChain->SetBranchAddress("Photon_pfRelIso03_all_quadratic", Photon_pfRelIso03_all_quadratic, &b_Photon_pfRelIso03_all_quadratic);
   fChain->SetBranchAddress("Photon_pfRelIso03_chg_quadratic", Photon_pfRelIso03_chg_quadratic, &b_Photon_pfRelIso03_chg_quadratic);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth, &b_Photon_phiWidth);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_s4", Photon_s4, &b_Photon_s4);
   fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_sieip", Photon_sieip, &b_Photon_sieip);
   fChain->SetBranchAddress("Photon_sipip", Photon_sipip, &b_Photon_sipip);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_x_calo", Photon_x_calo, &b_Photon_x_calo);
   fChain->SetBranchAddress("Photon_y_calo", Photon_y_calo, &b_Photon_y_calo);
   fChain->SetBranchAddress("Photon_z_calo", Photon_z_calo, &b_Photon_z_calo);

   fChain->SetBranchAddress("nPPSLocalTrack", &nPPSLocalTrack, &b_nPPSLocalTrack);
   fChain->SetBranchAddress("PPSLocalTrack_multiRPProtonIdx", &PPSLocalTrack_multiRPProtonIdx, &b_PPSLocalTrack_multiRPProtonIdx);
   fChain->SetBranchAddress("PPSLocalTrack_singleRPProtonIdx", &PPSLocalTrack_singleRPProtonIdx, &b_PPSLocalTrack_singleRPProtonIdx);
   fChain->SetBranchAddress("PPSLocalTrack_decRPId", &PPSLocalTrack_decRPId, &b_PPSLocalTrack_decRPId);
   fChain->SetBranchAddress("PPSLocalTrack_rpType", &PPSLocalTrack_rpType, &b_PPSLocalTrack_rpType);
   fChain->SetBranchAddress("PPSLocalTrack_x", &PPSLocalTrack_x, &b_PPSLocalTrack_x);
   fChain->SetBranchAddress("PPSLocalTrack_y", &PPSLocalTrack_y, &b_PPSLocalTrack_y);
   fChain->SetBranchAddress("PPSLocalTrack_time", &PPSLocalTrack_time, &b_PPSLocalTrack_time);
   fChain->SetBranchAddress("PPSLocalTrack_timeUnc", &PPSLocalTrack_timeUnc, &b_PPSLocalTrack_timeUnc);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_phiJERDown", &PuppiMET_phiJERDown, &b_PuppiMET_phiJERDown);
   fChain->SetBranchAddress("PuppiMET_phiJERUp", &PuppiMET_phiJERUp, &b_PuppiMET_phiJERUp);
   fChain->SetBranchAddress("PuppiMET_phiJESDown", &PuppiMET_phiJESDown, &b_PuppiMET_phiJESDown);
   fChain->SetBranchAddress("PuppiMET_phiJESUp", &PuppiMET_phiJESUp, &b_PuppiMET_phiJESUp);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredDown", &PuppiMET_phiUnclusteredDown, &b_PuppiMET_phiUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredUp", &PuppiMET_phiUnclusteredUp, &b_PuppiMET_phiUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_ptJERDown", &PuppiMET_ptJERDown, &b_PuppiMET_ptJERDown);
   fChain->SetBranchAddress("PuppiMET_ptJERUp", &PuppiMET_ptJERUp, &b_PuppiMET_ptJERUp);
   fChain->SetBranchAddress("PuppiMET_ptJESDown", &PuppiMET_ptJESDown, &b_PuppiMET_ptJESDown);
   fChain->SetBranchAddress("PuppiMET_ptJESUp", &PuppiMET_ptJESUp, &b_PuppiMET_ptJESUp);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredDown", &PuppiMET_ptUnclusteredDown, &b_PuppiMET_ptUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredUp", &PuppiMET_ptUnclusteredUp, &b_PuppiMET_ptUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   fChain->SetBranchAddress("RawPuppiMET_phi", &RawPuppiMET_phi, &b_RawPuppiMET_phi);
   fChain->SetBranchAddress("RawPuppiMET_pt", &RawPuppiMET_pt, &b_RawPuppiMET_pt);
   fChain->SetBranchAddress("RawPuppiMET_sumEt", &RawPuppiMET_sumEt, &b_RawPuppiMET_sumEt);
   fChain->SetBranchAddress("Rho_fixedGridRhoAll", &Rho_fixedGridRhoAll, &b_Rho_fixedGridRhoAll);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetAll", &Rho_fixedGridRhoFastjetAll, &b_Rho_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentral", &Rho_fixedGridRhoFastjetCentral, &b_Rho_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralCalo", &Rho_fixedGridRhoFastjetCentralCalo, &b_Rho_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralChargedPileUp", &Rho_fixedGridRhoFastjetCentralChargedPileUp, &b_Rho_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralNeutral", &Rho_fixedGridRhoFastjetCentralNeutral, &b_Rho_fixedGridRhoFastjetCentralNeutral);

   //fChain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   fChain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   fChain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   fChain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);

   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);

   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("OtherPV_score", OtherPV_score, &b_OtherPV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_charge", SV_charge, &b_SV_charge);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_dxy", SV_dxy, &b_SV_dxy);
   fChain->SetBranchAddress("SV_dxySig", SV_dxySig, &b_SV_dxySig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("SV_ntracks", SV_ntracks, &b_SV_ntracks);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);

   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);


   //fChain->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   
   //fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   //fChain->SetBranchAddress("L1_ZeroBias_copy", &L1_ZeroBias_copy, &b_L1_ZeroBias_copy);
   //fChain->SetBranchAddress("L1_UnprefireableEvent", &L1_UnprefireableEvent, &b_L1_UnprefireableEvent);
   //fChain->SetBranchAddress("L1Reco_step", &L1Reco_step, &b_L1Reco_step);

   fChain->SetBranchAddress("Flag_HBHENoiseFilter_pRECO", &Flag_HBHENoiseFilter_pRECO, &b_Flag_HBHENoiseFilter_pRECO);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter_pRECO", &Flag_HBHENoiseIsoFilter_pRECO, &b_Flag_HBHENoiseIsoFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter_pRECO", &Flag_CSCTightHaloFilter_pRECO, &b_Flag_CSCTightHaloFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO", &Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO, &b_Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter_pRECO", &Flag_CSCTightHalo2015Filter_pRECO, &b_Flag_CSCTightHalo2015Filter_pRECO);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter_pRECO", &Flag_globalTightHalo2016Filter_pRECO, &b_Flag_globalTightHalo2016Filter_pRECO);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter_pRECO", &Flag_globalSuperTightHalo2016Filter_pRECO, &b_Flag_globalSuperTightHalo2016Filter_pRECO);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter_pRECO", &Flag_HcalStripHaloFilter_pRECO, &b_Flag_HcalStripHaloFilter_pRECO);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter_pRECO", &Flag_hcalLaserEventFilter_pRECO, &b_Flag_hcalLaserEventFilter_pRECO);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO", &Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO, &b_Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter_pRECO", &Flag_EcalDeadCellBoundaryEnergyFilter_pRECO, &b_Flag_EcalDeadCellBoundaryEnergyFilter_pRECO);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter_pRECO", &Flag_ecalBadCalibFilter_pRECO, &b_Flag_ecalBadCalibFilter_pRECO);
   fChain->SetBranchAddress("Flag_goodVertices_pRECO", &Flag_goodVertices_pRECO, &b_Flag_goodVertices_pRECO);
   fChain->SetBranchAddress("Flag_eeBadScFilter_pRECO", &Flag_eeBadScFilter_pRECO, &b_Flag_eeBadScFilter_pRECO);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter_pRECO", &Flag_ecalLaserCorrFilter_pRECO, &b_Flag_ecalLaserCorrFilter_pRECO);
   fChain->SetBranchAddress("Flag_trkPOGFilters_pRECO", &Flag_trkPOGFilters_pRECO, &b_Flag_trkPOGFilters_pRECO);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter_pRECO", &Flag_chargedHadronTrackResolutionFilter_pRECO, &b_Flag_chargedHadronTrackResolutionFilter_pRECO);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter_pRECO", &Flag_muonBadTrackFilter_pRECO, &b_Flag_muonBadTrackFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter_pRECO", &Flag_BadChargedCandidateFilter_pRECO, &b_Flag_BadChargedCandidateFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter_pRECO", &Flag_BadPFMuonFilter_pRECO, &b_Flag_BadPFMuonFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter_pRECO", &Flag_BadPFMuonDzFilter_pRECO, &b_Flag_BadPFMuonDzFilter_pRECO);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter_pRECO", &Flag_hfNoisyHitsFilter_pRECO, &b_Flag_hfNoisyHitsFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter_pRECO", &Flag_BadChargedCandidateSummer16Filter_pRECO, &b_Flag_BadChargedCandidateSummer16Filter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter_pRECO", &Flag_BadPFMuonSummer16Filter_pRECO, &b_Flag_BadPFMuonSummer16Filter_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X_pRECO", &Flag_trkPOG_manystripclus53X_pRECO, &b_Flag_trkPOG_manystripclus53X_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X_pRECO", &Flag_trkPOG_toomanystripclus53X_pRECO, &b_Flag_trkPOG_toomanystripclus53X_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters_pRECO", &Flag_trkPOG_logErrorTooManyClusters_pRECO, &b_Flag_trkPOG_logErrorTooManyClusters_pRECO);
   fChain->SetBranchAddress("Flag_METFilters_pRECO", &Flag_METFilters_pRECO, &b_Flag_METFilters_pRECO);

   //fChain->SetBranchAddress("L1_AlwaysTrue_pRECO", &L1_AlwaysTrue_pRECO, &b_L1_AlwaysTrue_pRECO);

   //fChain->SetBranchAddress("L1_SingleEG10er2p5_pRECO", &L1_SingleEG10er2p5_pRECO, &b_L1_SingleEG10er2p5_pRECO);
   
   //fChain->SetBranchAddress("L1_UnpairedBunchBptxPlus_pRECO", &L1_UnpairedBunchBptxPlus_pRECO, &b_L1_UnpairedBunchBptxPlus_pRECO);
   //fChain->SetBranchAddress("L1_ZeroBias_pRECO", &L1_ZeroBias_pRECO, &b_L1_ZeroBias_pRECO);
   //fChain->SetBranchAddress("L1_ZeroBias_copy_pRECO", &L1_ZeroBias_copy_pRECO, &b_L1_ZeroBias_copy_pRECO);
   //fChain->SetBranchAddress("L1_UnprefireableEvent_pRECO", &L1_UnprefireableEvent_pRECO, &b_L1_UnprefireableEvent_pRECO);
   
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLT_EphemeralPhysics", &HLT_EphemeralPhysics, &b_HLT_EphemeralPhysics);
   fChain->SetBranchAddress("HLT_EphemeralZeroBias", &HLT_EphemeralZeroBias, &b_HLT_EphemeralZeroBias);
   fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
   fChain->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
   fChain->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
   fChain->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
   fChain->SetBranchAddress("HLT_Physics", &HLT_Physics, &b_HLT_Physics);
   fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   fChain->SetBranchAddress("HLT_ZeroBias_Alignment", &HLT_ZeroBias_Alignment, &b_HLT_ZeroBias_Alignment);
   fChain->SetBranchAddress("HLT_ZeroBias_Beamspot", &HLT_ZeroBias_Beamspot, &b_HLT_ZeroBias_Beamspot);
   fChain->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches, &b_HLT_ZeroBias_IsolatedBunches);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain, &b_HLT_ZeroBias_FirstBXAfterTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain, &b_HLT_ZeroBias_FirstCollisionInTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_LastCollisionInTrain", &HLT_ZeroBias_LastCollisionInTrain, &b_HLT_ZeroBias_LastCollisionInTrain);
   fChain->SetBranchAddress("HLT_HT300_Beamspot", &HLT_HT300_Beamspot, &b_HLT_HT300_Beamspot);
   fChain->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB, &b_HLT_IsoTrackHB);
   fChain->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE, &b_HLT_IsoTrackHE);

   //fChain->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID, &b_HLT_CaloJet500_NoJetID);
   //fChain->SetBranchAddress("HLT_CaloJet550_NoJetID", &HLT_CaloJet550_NoJetID, &b_HLT_CaloJet550_NoJetID);
   //fChain->SetBranchAddress("HLT_DoublePhoton33_CaloIdL", &HLT_DoublePhoton33_CaloIdL, &b_HLT_DoublePhoton33_CaloIdL);
   //fChain->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70, &b_HLT_DoublePhoton70);
   //fChain->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85, &b_HLT_DoublePhoton85);
 

   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX", &HLT_UncorrectedJetE30_NoBPTX, &b_HLT_UncorrectedJetE30_NoBPTX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX3BX", &HLT_UncorrectedJetE30_NoBPTX3BX, &b_HLT_UncorrectedJetE30_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE60_NoBPTX3BX", &HLT_UncorrectedJetE60_NoBPTX3BX, &b_HLT_UncorrectedJetE60_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE70_NoBPTX3BX", &HLT_UncorrectedJetE70_NoBPTX3BX, &b_HLT_UncorrectedJetE70_NoBPTX3BX);



   /*
   fChain->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40, &b_HLT_AK8PFJet40);
   fChain->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60, &b_HLT_AK8PFJet60);
   fChain->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80, &b_HLT_AK8PFJet80);
   fChain->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140, &b_HLT_AK8PFJet140);
   fChain->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200, &b_HLT_AK8PFJet200);
   fChain->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260, &b_HLT_AK8PFJet260);
   fChain->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320, &b_HLT_AK8PFJet320);
   fChain->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400, &b_HLT_AK8PFJet400);
   fChain->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
   fChain->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
   fChain->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550, &b_HLT_AK8PFJet550);
   fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   fChain->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
   fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   fChain->SetBranchAddress("HLT_PFJet110", &HLT_PFJet110, &b_HLT_PFJet110);
   fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550, &b_HLT_PFJet550);
   fChain->SetBranchAddress("HLT_PFJetFwd40", &HLT_PFJetFwd40, &b_HLT_PFJetFwd40);
   fChain->SetBranchAddress("HLT_PFJetFwd60", &HLT_PFJetFwd60, &b_HLT_PFJetFwd60);
   fChain->SetBranchAddress("HLT_PFJetFwd80", &HLT_PFJetFwd80, &b_HLT_PFJetFwd80);
   fChain->SetBranchAddress("HLT_PFJetFwd140", &HLT_PFJetFwd140, &b_HLT_PFJetFwd140);
   fChain->SetBranchAddress("HLT_PFJetFwd200", &HLT_PFJetFwd200, &b_HLT_PFJetFwd200);
   fChain->SetBranchAddress("HLT_PFJetFwd260", &HLT_PFJetFwd260, &b_HLT_PFJetFwd260);
   fChain->SetBranchAddress("HLT_PFJetFwd320", &HLT_PFJetFwd320, &b_HLT_PFJetFwd320);
   fChain->SetBranchAddress("HLT_PFJetFwd400", &HLT_PFJetFwd400, &b_HLT_PFJetFwd400);
   fChain->SetBranchAddress("HLT_PFJetFwd450", &HLT_PFJetFwd450, &b_HLT_PFJetFwd450);
   fChain->SetBranchAddress("HLT_PFJetFwd500", &HLT_PFJetFwd500, &b_HLT_PFJetFwd500);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd15", &HLT_AK8PFJetFwd15, &b_HLT_AK8PFJetFwd15);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd25", &HLT_AK8PFJetFwd25, &b_HLT_AK8PFJetFwd25);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd40", &HLT_AK8PFJetFwd40, &b_HLT_AK8PFJetFwd40);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd60", &HLT_AK8PFJetFwd60, &b_HLT_AK8PFJetFwd60);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd80", &HLT_AK8PFJetFwd80, &b_HLT_AK8PFJetFwd80);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd140", &HLT_AK8PFJetFwd140, &b_HLT_AK8PFJetFwd140);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd200", &HLT_AK8PFJetFwd200, &b_HLT_AK8PFJetFwd200);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd260", &HLT_AK8PFJetFwd260, &b_HLT_AK8PFJetFwd260);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd320", &HLT_AK8PFJetFwd320, &b_HLT_AK8PFJetFwd320);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd400", &HLT_AK8PFJetFwd400, &b_HLT_AK8PFJetFwd400);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd450", &HLT_AK8PFJetFwd450, &b_HLT_AK8PFJetFwd450);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd500", &HLT_AK8PFJetFwd500, &b_HLT_AK8PFJetFwd500);
   fChain->SetBranchAddress("HLT_PFHT180", &HLT_PFHT180, &b_HLT_PFHT180);
   fChain->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
   fChain->SetBranchAddress("HLT_PFHT370", &HLT_PFHT370, &b_HLT_PFHT370);
   fChain->SetBranchAddress("HLT_PFHT430", &HLT_PFHT430, &b_HLT_PFHT430);
   fChain->SetBranchAddress("HLT_PFHT510", &HLT_PFHT510, &b_HLT_PFHT510);
   fChain->SetBranchAddress("HLT_PFHT590", &HLT_PFHT590, &b_HLT_PFHT590);
   fChain->SetBranchAddress("HLT_PFHT680", &HLT_PFHT680, &b_HLT_PFHT680);
   fChain->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780, &b_HLT_PFHT780);
   fChain->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890, &b_HLT_PFHT890);
   fChain->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050, &b_HLT_PFHT1050);
   */

   fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight", &HLT_PFHT500_PFMET110_PFMHT110_IDTight, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight", &HLT_PFHT700_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight", &HLT_PFHT800_PFMET75_PFMHT75_IDTight, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight, &b_HLT_PFMET130_PFMHT130_IDTight);
   fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight, &b_HLT_PFMET140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight, &b_HLT_PFMETTypeOne140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_L1ETMHadSeeds", &HLT_L1ETMHadSeeds, &b_HLT_L1ETMHadSeeds);
   fChain->SetBranchAddress("HLT_CaloMHT90", &HLT_CaloMHT90, &b_HLT_CaloMHT90);
   fChain->SetBranchAddress("HLT_CaloMET90_NotCleaned", &HLT_CaloMET90_NotCleaned, &b_HLT_CaloMET90_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET350_NotCleaned", &HLT_CaloMET350_NotCleaned, &b_HLT_CaloMET350_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned, &b_HLT_PFMET200_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET250_NotCleaned", &HLT_PFMET250_NotCleaned, &b_HLT_PFMET250_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET300_NotCleaned", &HLT_PFMET300_NotCleaned, &b_HLT_PFMET300_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET200_BeamHaloCleaned", &HLT_PFMET200_BeamHaloCleaned, &b_HLT_PFMET200_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMETTypeOne200_BeamHaloCleaned", &HLT_PFMETTypeOne200_BeamHaloCleaned, &b_HLT_PFMETTypeOne200_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50, &b_HLT_MET105_IsoTrk50);
   fChain->SetBranchAddress("HLT_MET120_IsoTrk50", &HLT_MET120_IsoTrk50, &b_HLT_MET120_IsoTrk50);

  fChain->SetBranchAddress("HLT_DoublePFJets40_PFBTagDeepJet_p71", &HLT_DoublePFJets40_PFBTagDeepJet_p71, &b_HLT_DoublePFJets40_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets100_PFBTagDeepJet_p71", &HLT_DoublePFJets100_PFBTagDeepJet_p71, &b_HLT_DoublePFJets100_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets200_PFBTagDeepJet_p71", &HLT_DoublePFJets200_PFBTagDeepJet_p71, &b_HLT_DoublePFJets200_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets350_PFBTagDeepJet_p71", &HLT_DoublePFJets350_PFBTagDeepJet_p71, &b_HLT_DoublePFJets350_PFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71", &HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71, &b_HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71);
   fChain->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71", &HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71, &b_HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71);


   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);


   fChain->SetBranchAddress("HLT_Photon33", &HLT_Photon33, &b_HLT_Photon33);
   fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_Photon200", &HLT_Photon200, &b_HLT_Photon200);
   fChain->SetBranchAddress("HLT_Photon30EB_TightID_TightIso", &HLT_Photon30EB_TightID_TightIso, &b_HLT_Photon30EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon50EB_TightID_TightIso", &HLT_Photon50EB_TightID_TightIso, &b_HLT_Photon50EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon75EB_TightID_TightIso", &HLT_Photon75EB_TightID_TightIso, &b_HLT_Photon75EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon90EB_TightID_TightIso", &HLT_Photon90EB_TightID_TightIso, &b_HLT_Photon90EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso, &b_HLT_Photon110EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon130EB_TightID_TightIso", &HLT_Photon130EB_TightID_TightIso, &b_HLT_Photon130EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon150EB_TightID_TightIso", &HLT_Photon150EB_TightID_TightIso, &b_HLT_Photon150EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon175EB_TightID_TightIso", &HLT_Photon175EB_TightID_TightIso, &b_HLT_Photon175EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon200EB_TightID_TightIso", &HLT_Photon200EB_TightID_TightIso, &b_HLT_Photon200EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon100EBHE10", &HLT_Photon100EBHE10, &b_HLT_Photon100EBHE10);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95);
   fChain->SetBranchAddress("HLT_Photon35_TwoProngs35", &HLT_Photon35_TwoProngs35, &b_HLT_Photon35_TwoProngs35);


   //fChain->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40, &b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40);
   //fChain->SetBranchAddress("HLT_PFHT400_SixPFJet32", &HLT_PFHT400_SixPFJet32, &b_HLT_PFHT400_SixPFJet32);
   //fChain->SetBranchAddress("HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50", &HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50, &b_HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50);
   //fChain->SetBranchAddress("HLT_PFHT450_SixPFJet36", &HLT_PFHT450_SixPFJet36, &b_HLT_PFHT450_SixPFJet36);
   //fChain->SetBranchAddress("HLT_PFHT450_SixPFJet36_PNetBTag0p35", &HLT_PFHT450_SixPFJet36_PNetBTag0p35, &b_HLT_PFHT450_SixPFJet36_PNetBTag0p35);
   //fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_100_100_60_30_30", &HLT_PFHT400_FivePFJet_100_100_60_30_30, &b_HLT_PFHT400_FivePFJet_100_100_60_30_30);
   //fChain->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);

   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400);

   //fChain->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
   //fChain->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_HLT_DiSC30_18_EIso_AND_HE_Mass70);
   fChain->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose, &b_HLT_Photon20_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose, &b_HLT_Photon30_HoverELoose);

   fChain->SetBranchAddress("HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);



   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40", &HLT_AK8PFJet220_SoftDropMass40, &b_HLT_AK8PFJet220_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50", &HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50, &b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53", &HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53, &b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55", &HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55, &b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60", &HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60, &b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40, &b_HLT_AK8PFJet230_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06, &b_HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10", &HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10, &b_HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03", &HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03, &b_HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05", &HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05, &b_HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06, &b_HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10", &HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10, &b_HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03", &HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03, &b_HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05", &HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05, &b_HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06, &b_HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10", &HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10, &b_HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03", &HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03, &b_HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05", &HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05, &b_HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05);
   fChain->SetBranchAddress("HLT_AK8PFJet425_SoftDropMass40", &HLT_AK8PFJet425_SoftDropMass40, &b_HLT_AK8PFJet425_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet450_SoftDropMass40", &HLT_AK8PFJet450_SoftDropMass40, &b_HLT_AK8PFJet450_SoftDropMass40);



   //fChain->SetBranchAddress("HLT_HT350", &HLT_HT350, &b_HLT_HT350);
   //fChain->SetBranchAddress("HLT_HT425", &HLT_HT425, &b_HLT_HT425);
   //fChain->SetBranchAddress("HLT_CaloMET60_DTCluster50", &HLT_CaloMET60_DTCluster50, &b_HLT_CaloMET60_DTCluster50);
   //fChain->SetBranchAddress("HLT_CaloMET60_DTClusterNoMB1S50", &HLT_CaloMET60_DTClusterNoMB1S50, &b_HLT_CaloMET60_DTClusterNoMB1S50);
   //fChain->SetBranchAddress("HLT_L1MET_DTCluster50", &HLT_L1MET_DTCluster50, &b_HLT_L1MET_DTCluster50);
   //fChain->SetBranchAddress("HLT_L1MET_DTClusterNoMB1S50", &HLT_L1MET_DTClusterNoMB1S50, &b_HLT_L1MET_DTClusterNoMB1S50);
   fChain->SetBranchAddress("HLT_CscCluster_Loose", &HLT_CscCluster_Loose, &b_HLT_CscCluster_Loose);
   fChain->SetBranchAddress("HLT_CscCluster_Medium", &HLT_CscCluster_Medium, &b_HLT_CscCluster_Medium);
   fChain->SetBranchAddress("HLT_CscCluster_Tight", &HLT_CscCluster_Tight, &b_HLT_CscCluster_Tight);
   fChain->SetBranchAddress("HLT_DoubleCscCluster75", &HLT_DoubleCscCluster75, &b_HLT_DoubleCscCluster75);
   fChain->SetBranchAddress("HLT_DoubleCscCluster100", &HLT_DoubleCscCluster100, &b_HLT_DoubleCscCluster100);
   fChain->SetBranchAddress("HLT_L1CSCShower_DTCluster50", &HLT_L1CSCShower_DTCluster50, &b_HLT_L1CSCShower_DTCluster50);
   fChain->SetBranchAddress("HLT_L1CSCShower_DTCluster75", &HLT_L1CSCShower_DTCluster75, &b_HLT_L1CSCShower_DTCluster75);
   fChain->SetBranchAddress("HLT_PFMET105_IsoTrk50", &HLT_PFMET105_IsoTrk50, &b_HLT_PFMET105_IsoTrk50);
   fChain->SetBranchAddress("HLT_L1SingleLLPJet", &HLT_L1SingleLLPJet, &b_HLT_L1SingleLLPJet);


   fChain->SetBranchAddress("HLT_DiPhoton10Time1ns", &HLT_DiPhoton10Time1ns, &b_HLT_DiPhoton10Time1ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p2ns", &HLT_DiPhoton10Time1p2ns, &b_HLT_DiPhoton10Time1p2ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns, &b_HLT_DiPhoton10Time1p4ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p6ns", &HLT_DiPhoton10Time1p6ns, &b_HLT_DiPhoton10Time1p6ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p8ns", &HLT_DiPhoton10Time1p8ns, &b_HLT_DiPhoton10Time1p8ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time2ns", &HLT_DiPhoton10Time2ns, &b_HLT_DiPhoton10Time2ns);
   fChain->SetBranchAddress("HLT_DiPhoton10sminlt0p1", &HLT_DiPhoton10sminlt0p1, &b_HLT_DiPhoton10sminlt0p1);
   fChain->SetBranchAddress("HLT_DiPhoton10sminlt0p12", &HLT_DiPhoton10sminlt0p12, &b_HLT_DiPhoton10sminlt0p12);
   fChain->SetBranchAddress("HLT_DiPhoton10_CaloIdL", &HLT_DiPhoton10_CaloIdL, &b_HLT_DiPhoton10_CaloIdL);

   fChain->SetBranchAddress("HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT", &HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT, &b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT", &HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT, &b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);

   fChain->SetBranchAddress("HLT_IsoMu24_OneProng32", &HLT_IsoMu24_OneProng32, &b_HLT_IsoMu24_OneProng32);
   fChain->SetBranchAddress("HLT_Photon32_OneProng32_M50To105", &HLT_Photon32_OneProng32_M50To105, &b_HLT_Photon32_OneProng32_M50To105);


   /*
   fChain->SetBranchAddress("HLT_PFJet200_TimeLtNeg2p5ns", &HLT_PFJet200_TimeLtNeg2p5ns, &b_HLT_PFJet200_TimeLtNeg2p5ns);
   fChain->SetBranchAddress("HLT_PFJet200_TimeGt2p5ns", &HLT_PFJet200_TimeGt2p5ns, &b_HLT_PFJet200_TimeGt2p5ns);
   fChain->SetBranchAddress("HLT_Photon50_TimeLtNeg2p5ns", &HLT_Photon50_TimeLtNeg2p5ns, &b_HLT_Photon50_TimeLtNeg2p5ns);
   fChain->SetBranchAddress("HLT_Photon50_TimeGt2p5ns", &HLT_Photon50_TimeGt2p5ns, &b_HLT_Photon50_TimeGt2p5ns);
   fChain->SetBranchAddress("HLT_ExpressMuons", &HLT_ExpressMuons, &b_HLT_ExpressMuons);
   fChain->SetBranchAddress("HLT_PPSMaxTracksPerArm1", &HLT_PPSMaxTracksPerArm1, &b_HLT_PPSMaxTracksPerArm1);
   fChain->SetBranchAddress("HLT_PPSMaxTracksPerRP4", &HLT_PPSMaxTracksPerRP4, &b_HLT_PPSMaxTracksPerRP4);
   fChain->SetBranchAddress("HLT_PPSRandom", &HLT_PPSRandom, &b_HLT_PPSRandom);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   */
   Notify();
}

Bool_t PhotonJetAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PhotonJetAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PhotonJetAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PhotonJetAnalysis_cxx
