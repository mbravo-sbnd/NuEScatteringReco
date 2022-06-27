{
	/********************************************************************************************
	 * Macro that computes the reconstructed direction of the shower according to maximum 
	 * eigenvalue cut.
	 * Input:SegmentAndPortions.root ROOT file. 
	 * Output: ROOT histogram where angles between reco and true electron direction are stored
	 * (thetaRT variable, in degrees);
	 ********************************************************************************************/

	TFile *f = new TFile("SegmentAndPortions.root");
	TTree *t = (TTree*)f->Get("Tree");

	int Run, SubRun, ID;

	double VertX;
	double VertY;
	double VertZ;
	double NeutE;
	double ElecE;
	double TrueEDirX, TrueEDirY, TrueEDirZ;
	double TrueNeutDirX, TrueNeutDirY, TrueNeutDirZ;

	std::vector<double> *PointsX;
	std::vector<double> *PointsY;
	std::vector<double> *PointsZ;

	std::vector<int> *nseg;
	std::vector<double> *segXDir;
	std::vector<double> *segYDir;
	std::vector<double> *segZDir;
	std::vector<int> *nLastPoint;
	std::vector<double> *ZposLastPoint;

	std::vector<double> *cumulEValue;
	std::vector<double> *cumulThetaRT;
	std::vector<double> *cumulXDir;
	std::vector<double> *cumulYDir;
	std::vector<double> *cumulZDir;

	std::vector<double> *thetaWithPrim;
	std::vector<double> *meanThetaWithPrim;
	std::vector<double> *StdDevThetaWithPrim;

	std::vector<double> *PandTracksDirX;
	std::vector<double> *PandTracksDirY;
	std::vector<double> *PandTracksDirZ;

	std::vector<double> *PandShowersDirX;
	std::vector<double> *PandShowersDirY;
	std::vector<double> *PandShowersDirZ;

	t->SetBranchAddress("ID", &ID);
	t->SetBranchAddress("SubRun", &SubRun);
	t->SetBranchAddress("Run", &Run);
	t->SetBranchAddress("VertX", &VertX);
	t->SetBranchAddress("VertY", &VertY);
	t->SetBranchAddress("VertZ", &VertZ);
	t->SetBranchAddress("NeutE", &NeutE);
	t->SetBranchAddress("ElecE", &ElecE);
	t->SetBranchAddress("TrueEDirX", &TrueEDirX);
	t->SetBranchAddress("TrueEDirY", &TrueEDirY);
	t->SetBranchAddress("TrueEDirZ", &TrueEDirZ);
	t->SetBranchAddress("TrueNeutDirX", &TrueNeutDirX);
	t->SetBranchAddress("TrueNeutDirY", &TrueNeutDirY);
	t->SetBranchAddress("TrueNeutDirZ", &TrueNeutDirZ);

	t->SetBranchAddress("PointsX", &PointsX);
	t->SetBranchAddress("PointsY", &PointsY);
	t->SetBranchAddress("PointsZ", &PointsZ);

	t->SetBranchAddress("nSeg", &nseg);
	t->SetBranchAddress("segXDir", &segXDir);
	t->SetBranchAddress("segYDir", &segYDir);
	t->SetBranchAddress("segZDir", &segZDir);
	t->SetBranchAddress("nLastPoint", &nLastPoint);
	t->SetBranchAddress("ZposLastPoint", &ZposLastPoint);

	t->SetBranchAddress("cumulEValue", &cumulEValue);
	t->SetBranchAddress("cumulThetaRT", &cumulThetaRT);
	t->SetBranchAddress("cumulXDir", &cumulXDir);
	t->SetBranchAddress("cumulYDir", &cumulYDir);
	t->SetBranchAddress("cumulZDir", &cumulZDir);

	t->SetBranchAddress("thetaWithPrim", &thetaWithPrim);
	t->SetBranchAddress("meanThetaWithPrim", &meanThetaWithPrim);
	t->SetBranchAddress("StdDevThetaWithPrim", &StdDevThetaWithPrim);

	t->SetBranchAddress("PandTracksDirX", &PandTracksDirX);
	t->SetBranchAddress("PandTracksDirY", &PandTracksDirY);
	t->SetBranchAddress("PandTracksDirZ", &PandTracksDirZ);
	t->SetBranchAddress("PandShowersDirX", &PandShowersDirX);
	t->SetBranchAddress("PandShowersDirY", &PandShowersDirY);
	t->SetBranchAddress("PandShowersDirZ", &PandShowersDirZ);

	int i, j, nlastsegment;
	double thetaRT;

	TH1D *h = new TH1D("max #lambda_{n} cut", "#theta_{#hat{s}, e} with reco e direction (using reco points);#theta_{#hat{s}, e}(^{0});entries", 180, 0, 90);
	for (i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);

		//Looking for the shower portion to establish the cut
		nlastsegment = cumulEValue->size() - 1;
		for (j=0; j<cumulEValue->size(); j++){
			if (cumulEValue->at(j) > maxEValue){
				maxEValue = cumulEValue->at(j);
				nlastsegment = j;
			}
		}

		//Calculate theta (true, reco), with the direction reconstructed for the corresponding portion
		thetaRT = (180./TMath::Pi())*TMath::ACos(cumulXDir->at(nlastsegment)*ElecDirX + cumulYDir->at(nlastsegment)*ElecDirY + cumulZDir->at(nlastsegment)*ElecDirZ); 
		
		h->Fill(thetaRT);

	}//end entries loop

	//Obtaining mode value and FWHM of the histogram, to dsplay them on the screen
	double Mode = h->GetBinCenter(h->GetMaximumBin());
	double Lim1 = h->GetBinLowEdge(h->FindFirstBinAbove(h->GetMaximum()/2));
	double Lim2 = h->GetBinLowEdge(h->FindLastBinAbove(h->GetMaximum()/2)+1);
	double FWHM = Lim2 - Lim1;

	std::cout<<"ThetaRMS cut reco info:"<<std::endl;
	std::cout<<"Number of entries = "<<h->GetEntries()<<std::endl;
	std::cout<<"Mode value = "<<Mode<<"  FWHM = "<<FWHM<<std::endl;
	std::cout<<"Mean value = "<<h->GetMean()<<"  RMS = "<<h->GetRMS()<<std::endl;
	std::cout<<std::endl;

	TCanvas *c = new TCanvas("c", "c", 3000, 3000);
	h->Draw();




}
