#define PhotonJetAnalysis_cxx
#include "PhotonJetAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


//first histogram for testing purposes
TH1D *hn;

void PhotonJetAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PhotonJetAnalysis.C
//      root> PhotonJetAnalysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


	fChain->SetBranchStatus("*",0);	//switch off all branches
	fChain->SetBranchStatus("nJet",1);	//include branch in analysis
	fChain->SetBranchStatus("Jet_pt",1);	//include branch in analysis
	fChain->SetBranchStatus("Jet_eta",1);	//include branch in analysis
						
	fChain->SetBranchStatus("nPhoton",1);
	fChain->SetBranchStatus("Photon_pt",1);
	fChain->SetBranchStatus("Photon_eta",1);


	//create a first testing output file
	TFile *testfile;

	testfile = new TFile("photonjet_testfile.root", "RECREATE");




	//actual loop through events
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	//starting event loop
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

	if(jentry%1000==0){
		cout << "testing: jentry=" << jentry << endl;
		cout << "Jet pT: " << Jet_pt[jentry];
	}
		
		//hn->Fill(Jet_pt[jentry]);
      // if (Cut(ientry) < 0) continue;
   }
}
