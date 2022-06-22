{	
	/***************************************************************************************************************
	 * Macro that creates a ROOT output file with a TTree that stores the created 10-space-point
	 * segments and shower portions and their information, using analysisOutput.root as ROOT 
	 * input file with the events information.
	 * The output Tree has the following branches:
	 * - Run (int): Run number of the event
	 * - SubRun (int): Subun number of the event
	 * - ID (int): ID number of the event
	 * - VertX (double): x position of the true interaction vertex (in cm , SBND coordinate system)
	 * - VertY (double): y position of the true interaction vertex (in cm , SBND coordinate system)
	 * - VertZ (double): z position of the true interaction vertex (in cm , SBND coordinate system)
	 * - NeutE (double): true incoming neutrino energy (in GeV)
	 * - ElecE (double): true outgoing electron energy (in GeV)
	 * - TrueEDirX (double): x component of the true outgoing electron direction
	 * - TrueEDirY (double): y component of the true outgoing electron direction
	 * - TrueEDirZ (double): z component of the true outgoing electron direction
	 * - TrueNeutDirX (double): x component of the true incoming neutrino direction
	 * - TrueNeutDirY (double): y component of the true incoming neutrino direction
	 * - TrueNeutDirZ (double): z component of the true incoming neutrino direction
	 * - PointsX (vector<double>): x position of the reconstructed 3D space-points (in cm , SBND coordinate system)
	 * - PointsY (vector<double>): y position of the reconstructed 3D space-points (in cm , SBND coordinate system)
	 * - PointsZ (vector<double>): z position of the reconstructed 3D space-points (in cm , SBND coordinate system)
	 * - nseg (vector<int>): number each of the created 10-space-points segment
	 * - segXDir (vector<double>): x component of the direction of the created segments
	 * - segYDir (vector<double>): y component of the direction of the created segments
	 * - segZDir (vector<double>): z component of the direction of the created segments
	 * - nLastPoint (vector<int>): number (in Points vector) of the last point included in each segment
	 * - ZposLastPoint (vector<double>): z position (in cm , SBND coordinate system) of the last point included in 
	 * 	each segment
	 * - thetaWithPrim (vector<double>): angle between first segment and current segment direction (in degrees)
	 * - meanThetaWithPrim (vector<double>): mean angle(1st segment, n segment) up to current segment (in degrees)
	 * - StdDevThetaWithPrim (vector<double>): RMS angle(1st segment, n segment) up to current segment (in degrees^2)
	 * - cumulEValue (vector<double>): eigenvalue of the shower portion direction (which is an eigenvector)
	 * - cumulThetaRT (vector<double>): angle between true e direction and shower portion direction (degrees)
	 * - cumulXDir (vector<double>): x component of the shower portion direction
	 * - cumulYDir (vector<double>): y component of the shower portion direction
	 * - cumulZDir (vector<double>): z component of the shower portion direction
	 * - PandTracksDirX (vector<double>): x direction of the Pandora reconstructed tracks
	 * - PandTracksDirY (vector<double>): y direction of the Pandora reconstructed tracks
	 * - PandTracksDirZ (vector<double>): z direction of the Pandora reconstructed tracks
	 * - PandShowersDirX (vector<double>): x direction of the Pandora reconstructed showers
	 * - PandShowersDirY (vector<double>): y direction of the Pandora reconstructed showers
	 * - PandShowersDirZ (vector<double>): z direction of the Pandora reconstructed showers
	 * 
	 * All directions are calculated using PCA (Principal Component Analysis) tool, implemente in ROOT in TPrincipal
	 * class.
	 * *************************************************************************************************************/

	TFile *f=new TFile("NuE_SimEDepo.root");
	TTree *t=(TTree*)f->Get("ana/analysisOutputTree");

	int EventID, EventRun, EventSubRun;
	double TrueIntVerX, TrueIntVerY, TrueIntVerZ;

	std::vector<int> *trueParts;
	std::vector<double> *truePartsX;
	std::vector<double> *truePartsY;
	std::vector<double> *truePartsZ;
	std::vector<double> *truePartsPx;
	std::vector<double> *truePartsPy;
	std::vector<double> *truePartsPz;
	std::vector<double> *truePartsE;
	std::vector<int> *truePartsStatus;

	std::vector<double> *spX;
	std::vector<double> *spY;
	std::vector<double> *spZ;
	std::vector<double> *spHitsPTime;
	std::vector<double> *spHitsInteg;
	std::vector<double> *spHitsPlane;
	std::vector<double> *spHitsWire;
	std::vector<double> *spHitsTPC;

	std::vector<double> *G4EDepos;
	std::vector<double> *G4EDeposT;
	std::vector<double> *G4EDeposPDG;
	std::vector<double> *G4EDeposMidX;
	std::vector<double> *G4EDeposMidY;
	std::vector<double> *G4EDeposMidZ;

	std::vector<double> *showersXStartPoint;
	std::vector<double> *showersYStartPoint;
	std::vector<double> *showersZStartPoint;
	std::vector<double> *showersXDirection;
	std::vector<double> *showersYDirection;
	std::vector<double> *showersZDirection;

	std::vector<double> *tracksXDirection;
	std::vector<double> *tracksYDirection;
	std::vector<double> *tracksZDirection;

	
	t->SetBranchAddress("EventID", &EventID);
	t->SetBranchAddress("EventRun", &EventRun);
	t->SetBranchAddress("EventSubRun", &EventSubRun);
	t->SetBranchAddress("TrueIntVerX", &TrueIntVerX);
	t->SetBranchAddress("TrueIntVerY", &TrueIntVerY);
	t->SetBranchAddress("TrueIntVerZ", &TrueIntVerZ);
	t->SetBranchAddress("trueParts", &trueParts);
	t->SetBranchAddress("truePartsX", &truePartsX);
	t->SetBranchAddress("truePartsY", &truePartsY);
	t->SetBranchAddress("truePartsZ", &truePartsZ);
	t->SetBranchAddress("truePartsPx", &truePartsPx);
	t->SetBranchAddress("truePartsPy", &truePartsPy);
	t->SetBranchAddress("truePartsPz", &truePartsPz);
	t->SetBranchAddress("truePartsE", &truePartsE);
	t->SetBranchAddress("truePartsStatus", &truePartsStatus);
	t->SetBranchAddress("spX", &spX);
	t->SetBranchAddress("spY", &spY);
	t->SetBranchAddress("spZ", &spZ);
	t->SetBranchAddress("spHitsPTime", &spHitsPTime);
	t->SetBranchAddress("spHitsInteg", &spHitsInteg);
	t->SetBranchAddress("spHitsPlane", &spHitsPlane);
	t->SetBranchAddress("spHitsWire", &spHitsWire);
	t->SetBranchAddress("spHitsTPC", &spHitsTPC);

	t->SetBranchAddress("showersXStartPoint", &showersXStartPoint);
	t->SetBranchAddress("showersYStartPoint", &showersYStartPoint);
	t->SetBranchAddress("showersZStartPoint", &showersZStartPoint);
	t->SetBranchAddress("showersXDirection", &showersXDirection);
	t->SetBranchAddress("showersYDirection", &showersYDirection);
	t->SetBranchAddress("showersZDirection", &showersZDirection);

	t->SetBranchAddress("tracksXDirection", &tracksXDirection);
	t->SetBranchAddress("tracksYDirection", &tracksYDirection);
	t->SetBranchAddress("tracksZDirection", &tracksZDirection);
	
	

	TFile *fSave = new TFile("SegmentAndPortions.root", "RECREATE");
	TTree *tSave = new TTree("Tree", "Tree");

	int Run, SubRun, ID;

	double VertX;
	double VertY;
	double VertZ;
	double NeutE;
	double ElecE;
	double TrueEDirX, TrueEDirY, TrueEDirZ;
	double TrueNeutDirX, TrueNeutDirY, TrueNeutDirZ;


	std::vector<double> PointsX;
	std::vector<double> PointsY;
	std::vector<double> PointsZ;

	std::vector<int> nseg;
	std::vector<double> segXDir;
	std::vector<double> segYDir;
	std::vector<double> segZDir;
	std::vector<int> nLastPoint;
	std::vector<double> ZposLastPoint;

	std::vector<double> cumulEValue;
	std::vector<double> cumulThetaRT;
	std::vector<double> cumulXDir;
	std::vector<double> cumulYDir;
	std::vector<double> cumulZDir;

	std::vector<double> thetaWithPrim;
	std::vector<double> meanThetaWithPrim;
	std::vector<double> StdDevThetaWithPrim;

	std::vector<double> PandTracksDirX;
	std::vector<double> PandTracksDirY;
	std::vector<double> PandTracksDirZ;

	std::vector<double> PandShowersDirX;
	std::vector<double> PandShowersDirY;
	std::vector<double> PandShowersDirZ;

	tSave->Branch("ID", &ID);
	tSave->Branch("SubRun", &SubRun);
	tSave->Branch("Run", &Run);
	tSave->Branch("VertX", &VertX);
	tSave->Branch("VertY", &VertY);
	tSave->Branch("VertZ", &VertZ);
	tSave->Branch("NeutE", &NeutE);
	tSave->Branch("ElecE", &ElecE);
	tSave->Branch("TrueEDirX", &TrueEDirX);
	tSave->Branch("TrueEDirY", &TrueEDirY);
	tSave->Branch("TrueEDirZ", &TrueEDirZ);
	tSave->Branch("TrueNeutDirX", &TrueNeutDirX);
	tSave->Branch("TrueNeutDirY", &TrueNeutDirY);
	tSave->Branch("TrueNeutDirZ", &TrueNeutDirZ);

	tSave->Branch("PointsX", &PointsX);
	tSave->Branch("PointsY", &PointsY);
	tSave->Branch("PointsZ", &PointsZ);
	tSave->Branch("hPlane", &hPlane);
	tSave->Branch("hPTime", &hPTime);
	tSave->Branch("hInteg", &hInteg);

	tSave->Branch("nSeg", &nseg);
	tSave->Branch("segXDir", &segXDir);
	tSave->Branch("segYDir", &segYDir);
	tSave->Branch("segZDir", &segZDir);
	tSave->Branch("nLastPoint", &nLastPoint);
	tSave->Branch("ZposLastPoint", &ZposLastPoint);

	tSave->Branch("cumulEValue", &cumulEValue);
	tSave->Branch("cumulThetaRT", &cumulThetaRT);
	tSave->Branch("cumulXDir", &cumulXDir);
	tSave->Branch("cumulYDir", &cumulYDir);
	tSave->Branch("cumulZDir", &cumulZDir);

	tSave->Branch("thetaWithPrim", &thetaWithPrim);
	tSave->Branch("meanThetaWithPrim", &meanThetaWithPrim);
	tSave->Branch("StdDevThetaWithPrim", &StdDevThetaWithPrim);

	tSave->Branch("PandTracksDirX", &PandTracksDirX);
	tSave->Branch("PandTracksDirY", &PandTracksDirY);
	tSave->Branch("PandTracksDirZ", &PandTracksDirZ);
	tSave->Branch("PandShowersDirX", &PandShowersDirX);
	tSave->Branch("PandShowersDirY", &PandShowersDirY);
	tSave->Branch("PandShowersDirZ", &PandShowersDirZ);

	int i, j, k, nSegsMax, nHitMax; 
	double thetai, meanthetai, DirX, DirY, DirZ, DirXPrim, DirYPrim, DirZPrim;
	double sumThetaiSquare;
	double rX, rY, rZ, cDirX, cDirY, cDirZ, cThetaRT;

	for (i=0; i<t->GetEntries(); i++){

		t->GetEntry(i);

		if (spX->size() < 20) continue;

		PointsX.clear();
		PointsY.clear();
		PointsZ.clear();
		hPlane.clear();
		hPTime.clear();
		hInteg.clear();

		nseg.clear();
		segXDir.clear();
		segYDir.clear();
		segZDir.clear();
		nLastPoint.clear();
		ZposLastPoint.clear();

		cumulEValue.clear();
		cumulThetaRT.clear();
		cumulXDir.clear();
		cumulYDir.clear();
		cumulZDir.clear();

		thetaWithPrim.clear();
		meanThetaWithPrim.clear();
		StdDevThetaWithPrim.clear();

		PandTracksDirX.clear();
		PandTracksDirY.clear();
		PandTracksDirZ.clear();
		PandShowersDirX.clear();
		PandShowersDirY.clear();
		PandShowersDirZ.clear();

		sumThetaiSquare = 0.;

		//Storing variables directly from input TTree

		ID=EventID;
		SubRun=EventSubRun;
		Run=EventRun;
		

		for (j=0; j<truePartsStatus->size(); j++){

			if (truePartsStatus->at(j)==1 && trueParts->at(j)==11){
				VertX=truePartsX->at(j);
				VertY=truePartsY->at(j);
				VertZ=truePartsZ->at(j);

				double peX=truePartsPx->at(j);
				double peY=truePartsPy->at(j);
				double peZ=truePartsPz->at(j);
				double pe=TMath::Sqrt(peX*peX + peY*peY + peZ*peZ);
				TrueEDirX=peX/pe;
				TrueEDirY=peY/pe;
				TrueEDirZ=peZ/pe;

				ElecE=truePartsE->at(j);
			}


			if (truePartsStatus->at(j)==0 && (TMath::Abs(trueParts->at(j))==12 || TMath::Abs(trueParts->at(j))==14)){
				
				double pnX=truePartsPx->at(j);
				double pnY=truePartsPy->at(j);
				double pnZ=truePartsPz->at(j);
				double pn=TMath::Sqrt(pnX*pnX + pnY*pnY + pnZ*pnZ);
				TrueNeutDirX=pnX/pn;
				TrueNeutDirY=pnY/pn;
				TrueNeutDirZ=pnZ/pn;

				NeutE=truePartsE->at(j);
			}
		}

		for (j=0; j<showersXDirection->size(); j++){
			PandShowersDirX.push_back(showersXDirection->at(j));
			PandShowersDirY.push_back(showersYDirection->at(j));
			PandShowersDirZ.push_back(showersZDirection->at(j));
		}

		for (j=0; j<tracksXDirection->size(); j++){
			PandTracksDirX.push_back(tracksXDirection->at(j));
			PandTracksDirY.push_back(tracksYDirection->at(j));
			PandTracksDirZ.push_back(tracksZDirection->at(j));
		}

		for (j=0; j<spX->size(); j++){
			PointsX.push_back(spX->at(j));
			PointsY.push_back(spY->at(j));
			PointsZ.push_back(spZ->at(j));
		}

		//Calculating number of shower segments and portions
		nSegsMax=int(spZ->size()/10.);


		for (j=0; j<nSegsMax; j++){
			
			int nHitMin = 10*j;
			nHitMax=10*(j+1);
			nseg.push_back(j);

			//Variables for j-th segment
			double pointsXYZ[3];
        	TPrincipal *princ=new TPrincipal(3, "");

			for (k=nHitMin; k<nHitMax; k++){
				if (k > spX->size() - 1) break;
				pointsXYZ[0]=spX->at(k);
				pointsXYZ[1]=spY->at(k);
				pointsXYZ[2]=spZ->at(k);

				princ->AddRow(pointsXYZ);

				if (k == nHitMax-1){
					nLastPoint.push_back(k);
					ZposLastPoint.push_back(spZ->at(k));
				}

			}//end hits-in-the-segment loop

			princ->MakePrincipals();
			const TMatrixD *eigenMatrix = princ->GetEigenVectors();
			DirX=eigenMatrix->operator()(0,0);
			DirY=eigenMatrix->operator()(1,0);
			DirZ=eigenMatrix->operator()(2,0);

			double mod = TMath::Sqrt(DirX*DirX + DirY*DirY + DirZ*DirZ);

			DirX = DirX/mod;
			DirY = DirY/mod;
			DirZ = DirZ/mod;

			if (DirX*DirXPrim<0) DirX=-DirX;
			if (DirY*DirYPrim<0) DirY=-DirY;
			if (DirZ*DirZPrim<0) DirZ=-DirZ;

			segXDir.push_back(DirX);
			segYDir.push_back(DirY);
			segZDir.push_back(DirZ);

			//Segments angles obtention *************************************************************
			if (j==0){

				thetaWithPrim.push_back(-1);
				meanThetaWithPrim.push_back(-1);
				StdDevThetaWithPrim.push_back(0);

				DirXPrim = DirX;
				DirYPrim = DirY;
				DirZPrim = DirZ;
	

			}//end is-first-segment condition

			else if (j==1){

				thetai = (180./TMath::Pi())*TMath::ACos(DirX*DirXPrim + DirY*DirYPrim + DirZ*DirZPrim);
				sumThetaiSquare = sumThetaiSquare + thetai*thetai;

				thetaWithPrim.push_back(thetai);
				meanThetaWithPrim.push_back(thetai);
				StdDevThetaWithPrim.push_back(0);

			}//end is-second-segment condition

			else if (j>1){

				thetai = (180./TMath::Pi())*TMath::ACos(DirX*DirXPrim + DirY*DirYPrim + DirZ*DirZPrim);
				sumThetaiSquare = sumThetaiSquare + thetai*thetai;
				thetaWithPrim.push_back(thetai);

				meanthetai = ((j-1)*meanThetaWithPrim.at(meanThetaWithPrim.size()-1) + thetai)/j;
				meanThetaWithPrim.push_back(meanthetai);

				StdDevThetaWithPrim.push_back(sumThetaiSquare/j - meanthetai*meanthetai);

			}//end is-more-than-second-segment condition

			//***************************************************************************************+**


			//Variables for j-th shower portion
			if (nHitMax < spX->size()-1){
				rX = spX->at(nHitMax - 1) - spX->at(0);
				rY = spY->at(nHitMax - 1) - spY->at(0);
				rZ = spZ->at(nHitMax - 1) - spZ->at(0);
			}
			else if (nHitMax > spX->size()-1){
				rX = spX->back() - spX->at(0);
				rY = spY->back() - spY->at(0);
				rZ = spZ->back() - spZ->at(0);
			}
			

			TPrincipal *princ2=new TPrincipal(3, "");
			
			double pointsXYZ2[3];
			for (k=0; k<nHitMax; k++){
				if (k > spX->size() - 1) break;
				pointsXYZ2[0]=spX->at(k);
				pointsXYZ2[1]=spY->at(k);
				pointsXYZ2[2]=spZ->at(k);

				princ2->AddRow(pointsXYZ2);
			}

			princ2->MakePrincipals();
			const TMatrixD *ceigenMatrix = princ2->GetEigenVectors();
			cDirX=ceigenMatrix->operator()(0,0);
			cDirY=ceigenMatrix->operator()(1,0);
			cDirZ=ceigenMatrix->operator()(2,0);

			const TVectorD *ceValues = princ2->GetEigenValues();
			cumulEValue.push_back(ceValues->operator()(0));

			princ2->Clear();

			double cmod = TMath::Sqrt(cDirX*cDirX + cDirY*cDirY + cDirZ*cDirZ);

			cDirX = cDirX/cmod;
			cDirY = cDirY/cmod;
			cDirZ = cDirZ/cmod;

			if (cDirX*rX<0) cDirX=-cDirX;
			if (cDirY*rY<0) cDirY=-cDirY;
			if (cDirZ*rZ<0) cDirZ=-cDirZ;

			cumulXDir.push_back(cDirX);
			cumulYDir.push_back(cDirY);
			cumulZDir.push_back(cDirZ);

			cThetaRT = (180./TMath::Pi())*TMath::ACos(cDirX*TrueEDirX + cDirY*TrueEDirY + cDirZ*TrueEDirZ);
			cumulThetaRT.push_back(cThetaRT);

		}//end loop-over-possible-segments loop

		std::cout<<"Event "<<Run<<"."<<SubRun<<"."<<ID<<" done."<<std::endl;
		tSave->Fill();


	}//end entries loop

	tSave->Write();
	fSave->Close();
	f->Close();














































/*
 IN CASE YOU HAVE TO CHECK IF HITS ARE CORRECTLY ODERED
 int nHits1, nHits2, nHits3, nHits4, nHits5, nHits6;

	TGraph *gXZ1 = new TGraph();
	gXZ1->SetLineWidth(0);
	gXZ1->SetMarkerStyle(8);
	gXZ1->SetMarkerSize(2);
	gXZ1->SetMarkerColor(kRed+3);

	TGraph *gYZ1 = new TGraph();
	gYZ1->SetLineWidth(0);
	gYZ1->SetMarkerStyle(8);
	gYZ1->SetMarkerSize(2);
	gYZ1->SetMarkerColor(kRed+3);

	TGraph *gXZ2 = new TGraph();
	gXZ2->SetLineWidth(0);
	gXZ2->SetMarkerStyle(8);
	gXZ2->SetMarkerSize(2);
	gXZ2->SetMarkerColor(kRed);

	TGraph *gYZ2 = new TGraph();
	gYZ2->SetLineWidth(0);
	gYZ2->SetMarkerStyle(8);
	gYZ2->SetMarkerSize(2);
	gYZ2->SetMarkerColor(kRed);

	TGraph *gXZ3 = new TGraph();
	gXZ3->SetLineWidth(0);
	gXZ3->SetMarkerStyle(8);
	gXZ3->SetMarkerSize(2);
	gXZ3->SetMarkerColor(kOrange);

	TGraph *gYZ3 = new TGraph();
	gYZ3->SetLineWidth(0);
	gYZ3->SetMarkerStyle(8);
	gYZ3->SetMarkerSize(2);
	gYZ3->SetMarkerColor(kOrange);

	TGraph *gXZ4 = new TGraph();
	gXZ4->SetLineWidth(0);
	gXZ4->SetMarkerStyle(8);
	gXZ4->SetMarkerSize(2);
	gXZ4->SetMarkerColor(kYellow);

	TGraph *gYZ4 = new TGraph();
	gYZ4->SetLineWidth(0);
	gYZ4->SetMarkerStyle(8);
	gYZ4->SetMarkerSize(2);
	gYZ4->SetMarkerColor(kYellow);

	TGraph *gXZ5 = new TGraph();
	gXZ5->SetLineWidth(0);
	gXZ5->SetMarkerStyle(8);
	gXZ5->SetMarkerSize(2);
	gXZ5->SetMarkerColor(kGreen);

	TGraph *gYZ5 = new TGraph();
	gYZ5->SetLineWidth(0);
	gYZ5->SetMarkerStyle(8);
	gYZ5->SetMarkerSize(2);
	gYZ5->SetMarkerColor(kGreen);

	TGraph *gXZ6 = new TGraph();
	gXZ6->SetLineWidth(0);
	gXZ6->SetMarkerStyle(8);
	gXZ6->SetMarkerSize(2);
	gXZ6->SetMarkerColor(kBlue);

	TGraph *gYZ6 = new TGraph();
	gYZ6->SetLineWidth(0);
	gYZ6->SetMarkerStyle(8);
	gYZ6->SetMarkerSize(2);
	gYZ6->SetMarkerColor(kBlue);

	t->GetEntry(2);

	std::cout<<"number of hits for first segment (n total hits = "<<PointsX->size()<<"): ";
	std::cin>>nHits1;
	std::cout<<std::endl;

	std::cout<<"number of hits for second segment (n total hits = "<<PointsX->size()<<"): ";
	std::cin>>nHits2;
	std::cout<<std::endl;

	std::cout<<"number of hits for third segment (n total hits = "<<PointsX->size()<<"): ";
	std::cin>>nHits3;
	std::cout<<std::endl;

	std::cout<<"number of hits for fourth segment (n total hits = "<<PointsX->size()<<"): ";
	std::cin>>nHits4;
	std::cout<<std::endl;

	std::cout<<"number of hits for fifth segment (n total hits = "<<PointsX->size()<<"): ";
	std::cin>>nHits5;
	std::cout<<std::endl;

	std::cout<<"number of hits for sixth segment (n total hits = "<<PointsX->size()<<"): ";
	std::cin>>nHits6;
	std::cout<<std::endl;
	
	for (j=0; j<nHits1; j++){

		gXZ1->SetPoint(gXZ1->GetN(), PointsZ->at(j), PointsX->at(j));
		gYZ1->SetPoint(gYZ1->GetN(), PointsZ->at(j), PointsY->at(j));

	}

	for (j=nHits1; j<nHits2; j++){

		gXZ2->SetPoint(gXZ2->GetN(), PointsZ->at(j), PointsX->at(j));
		gYZ2->SetPoint(gYZ2->GetN(), PointsZ->at(j), PointsY->at(j));

	}

	for (j=nHits2; j<nHits3; j++){

		gXZ3->SetPoint(gXZ3->GetN(), PointsZ->at(j), PointsX->at(j));
		gYZ3->SetPoint(gYZ3->GetN(), PointsZ->at(j), PointsY->at(j));

	}

	for (j=nHits3; j<nHits4; j++){

		gXZ4->SetPoint(gXZ4->GetN(), PointsZ->at(j), PointsX->at(j));
		gYZ4->SetPoint(gYZ4->GetN(), PointsZ->at(j), PointsY->at(j));

	}

	for (j=nHits4; j<nHits5; j++){

		gXZ5->SetPoint(gXZ5->GetN(), PointsZ->at(j), PointsX->at(j));
		gYZ5->SetPoint(gYZ5->GetN(), PointsZ->at(j), PointsY->at(j));

	}

	for (j=nHits5; j<nHits6; j++){

		gXZ6->SetPoint(gXZ6->GetN(), PointsZ->at(j), PointsX->at(j));
		gYZ6->SetPoint(gYZ6->GetN(), PointsZ->at(j), PointsY->at(j));

	}

	TMultiGraph *mgXZ = new TMultiGraph();
	mgXZ->Add(gXZ1);
	mgXZ->Add(gXZ2);
	mgXZ->Add(gXZ3);
	mgXZ->Add(gXZ4);
	mgXZ->Add(gXZ5);
	mgXZ->Add(gXZ6);


	TMultiGraph *mgYZ = new TMultiGraph();
	mgYZ->Add(gYZ1);
	mgYZ->Add(gYZ2);
	mgYZ->Add(gYZ3);
	mgYZ->Add(gYZ4);
	mgYZ->Add(gYZ5);
	mgYZ->Add(gYZ6);

	TCanvas *cXZ = new TCanvas("cXZ", "cXZ", 3000, 3000);
	mgXZ->Draw("APL");

	TCanvas *cYZ = new TCanvas("cYZ", "cYZ", 3000, 3000);
	mgYZ->Draw("APL");

*/





}