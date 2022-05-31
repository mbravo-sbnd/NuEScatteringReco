/**************************************************************************************************
 * Macro that creates a ROOT file where the segments of the tracks of the muon 
 * sample are stored (for using in Kalman Filter reconstruction).
 * Input analysisOutput.root with the muon sample + length for the segment (lmax),
 * asked to the user.
 * Output: SegmentedMuon.root file, with the TTree called Tree. It has the following branches:
 * 
 * 	- Run (int): Run number of the event
	- SubRun (int): Sub-Run number of the event
	- ID (int): ID number of the event
	- MuE (double): MC initial energy of the muon (MeV)
	- MuP (double): MC initial momentum of the muon (MeV/c)
	- MuDepoE (double): MC muon total deposited energy (MeV)
	- MuDirX (double): MC X component of the muon inital direction 
	- MuDirY (double): MC Y component of the muon inital direction 
	- MuDirZ (double): MC Z component of the muon inital direction  
	- MuStartX (double): MC X component of the muon start point (cm, SBND coordinate system) 
	- MuStartY (double): MC Y component of the muon start point (cm, SBND coordinate system)
	- MuStartZ (double): MC Z component of the muon start point (cm, SBND coordinate system) 
	- IsEscaping (int): 0 = it is a contained muon, 1 = it is an escaping muon
	- spX (double): X position of each reconstructed space-point (cm, SBND coordinate system) 
	- spY (double): Y position of each reconstructed space-point (cm, SBND coordinate system) 
	- spZ (double): Z position of each reconstructed space-point (cm, SBND coordinate system) 
	- spInteg (double): integral of the hit associated to each space-point (ADC x tick)
	- spPlane (double): number of the plane if the hit associated to each spae-point
	- EdepoX (double): X position of each MC energy deposition (cm, SBND coordinate system) 
	- EdepoY (double): Y position of each MC energy deposition (cm, SBND coordinate system) 
	- EdepoZ (double): Z position of each MC energy deposition (cm, SBND coordinate system) 
	- Edepo (double): MC energy depositions (MeV) 
	- segX (double): X position of the origin of each created segment (cm) 
	- segY (double): Y position of the origin of each created segment (cm) 
	- segZ (double): Z position of the origin of each created segment (cm) 
	- segdxdz (double): dx/dz of each created segment 
	- segdydz (double): dy/dz of each created segment 
	- segIncz (double): (z_final - z_initial) of each created segment (cm)
	- segXdir (double): X component of the direction of each created segment 
	- segYdir (double): Y component of the direction of each created segment 
	- segZdir (double): Z component of the direction of each created segment 
	- segXdirErr (double): error in X component of the direction of each created segment 
	- segYdirErr (double): error in Y component of the direction of each created segment
	- segZdirErr (double): error in Z component of the direction of each created segment
	- segLength (double): (precise) length of each segment (cm)
	- segDepoE (double): deposited reconstructed energy in each segment, with calorimetry
		relation (MeV) 
	- segTrueDepoE (double): deposited MC energy in each segment (MeV)
 * 
 * For direction reconstruction, PCA method is used. ROOT::TPrincipal is the object used
 * for this purpose.
 * *******************************************************************************************/


{

	TFile *f = new TFile("analysisOutput.root");
	TTree *t=(TTree*)f->Get("ana/analysisOutputTree");

	int Run, SubRun, ID;
	std::vector<double> *MCPartsPx;
	std::vector<double> *MCPartsPy;
	std::vector<double> *MCPartsPz;
	std::vector<double> *MCPartsFinalPx;
	std::vector<double> *MCPartsFinalPy;
	std::vector<double> *MCPartsFinalPz;
	std::vector<double> *MCPartsX;
	std::vector<double> *MCPartsY;
	std::vector<double> *MCPartsZ;
	std::vector<double> *MCPartsE;
	std::vector<double> *G4EDepos;
	std::vector<int> *G4EDeposPDG;
	std::vector<double> *spX;
	std::vector<double> *spY;
	std::vector<double> *spZ;
	std::vector<double> *spHitsInteg;
	std::vector<int> *spHitsPlane;
	std::vector<double> *G4EDeposMidX;
	std::vector<double> *G4EDeposMidY;
	std::vector<double> *G4EDeposMidZ;

	t->SetBranchAddress("Run", &Run);
	t->SetBranchAddress("SubRun", &SubRun);
	t->SetBranchAddress("ID", &ID);
	t->SetBranchAddress("MCPartsX", &MCPartsX);
	t->SetBranchAddress("MCPartsY", &MCPartsY);
	t->SetBranchAddress("MCPartsZ", &MCPartsZ);
	t->SetBranchAddress("MCPartsPx", &MCPartsPx);
	t->SetBranchAddress("MCPartsPy", &MCPartsPy);
	t->SetBranchAddress("MCPartsPz", &MCPartsPz);
	t->SetBranchAddress("MCPartsFinalPx", &MCPartsFinalPx);
	t->SetBranchAddress("MCPartsFinalPy", &MCPartsFinalPy);
	t->SetBranchAddress("MCPartsFinalPz", &MCPartsFinalPz);
	t->SetBranchAddress("MCPartsE", &MCPartsE);
	t->SetBranchAddress("G4EDepos", &G4EDepos);
	t->SetBranchAddress("G4EDeposPDG", &G4EDeposPDG);
	t->SetBranchAddress("spX", &spX);
	t->SetBranchAddress("spY", &spY);
	t->SetBranchAddress("spZ", &spZ);
	t->SetBranchAddress("spHitsInteg", &spHitsInteg);
	t->SetBranchAddress("spHitsPlane", &spHitsPlane);
	t->SetBranchAddress("G4EDeposMidX", &G4EDeposMidX);
	t->SetBranchAddress("G4EDeposMidY", &G4EDeposMidY);
	t->SetBranchAddress("G4EDeposMidZ", &G4EDeposMidZ);

	//Creating a new .root file with a new tree with the muon information

	TFile *fSave = new TFile("SegmentedMuon.root", "RECREATE");
	TTree *tSave=new TTree("Tree", "Tree");

	int Run1, SubRun1, ID1;
	int IsEscaping;
	double MuE, MuP, MuDepoE, MuDirX, MuDirY, MuDirZ, MuStartX, MuStartY, MuStartZ;
	std::vector<double> spX1, spY1, spZ1, spInteg1, spPlane1;
	std::vector<double>  segX, segY, segZ, segdxdz, segdydz, segIncz, segXdir, segYdir, segZdir;
	std::vector<double>  segDepoE, segLength, segTrueDepoE;
	std::vector<double> segXdirErr, segYdirErr, segZdirErr;
	std::vector<double> EdepoX, EdepoY, EdepoZ, Edepo;

	tSave->Branch("Run", &Run1, "Run/I");
	tSave->Branch("SubRun", &SubRun1, "SubRun/I");
	tSave->Branch("ID", &ID1, "ID/I");
	tSave->Branch("MuE", &MuE, "MuE/D");
	tSave->Branch("MuP", &MuP, "MuP/D");
	tSave->Branch("MuDepoE", &MuDepoE, "MuDepoE/D");
	tSave->Branch("MuDirX", &MuDirX, "MuDirX/D");
	tSave->Branch("MuDirY", &MuDirY, "MuDirY/D");
	tSave->Branch("MuDirZ", &MuDirZ, "MuDirZ/D");
	tSave->Branch("MuStartX", &MuStartX, "MuStartX/D");
	tSave->Branch("MuStartY", &MuStartY, "MuStartY/D");
	tSave->Branch("MuStartZ", &MuStartZ, "MuStartZ/D");
	tSave->Branch("IsEscaping", &IsEscaping, "IsEscaping/I");
	tSave->Branch("spX", &spX1);
	tSave->Branch("spY", &spY1);
	tSave->Branch("spZ", &spZ1);
	tSave->Branch("spInteg", &spInteg1);
	tSave->Branch("spPlane", &spPlane1);
	tSave->Branch("EdepoX", &EdepoX);
	tSave->Branch("EdepoY", &EdepoY);
	tSave->Branch("EdepoZ", &EdepoZ);
	tSave->Branch("Edepo", &Edepo);
	tSave->Branch("segX", &segX);
	tSave->Branch("segY", &segY);
	tSave->Branch("segZ", &segZ);
	tSave->Branch("segdxdz", &segdxdz);
	tSave->Branch("segdydz", &segdydz);
	tSave->Branch("segIncz", &segIncz);
	tSave->Branch("segXdir", &segXdir);
	tSave->Branch("segYdir", &segYdir);
	tSave->Branch("segZdir", &segZdir);
	tSave->Branch("segXdirErr", &segXdirErr);
	tSave->Branch("segYdirErr", &segYdirErr);
	tSave->Branch("segZdirErr", &segZdirErr);
	tSave->Branch("segLength", &segLength);
	tSave->Branch("segDepoE", &segDepoE);
	tSave->Branch("segTrueDepoE", &segTrueDepoE);


	int i, j, k, m, a, b, It;
	double E, xdirAux, ydirAux, zdirAux, truedepoE;
	double seglim1x, seglim2x, seglim1y, seglim2y, seglim1z, seglim2z;

	bool bFound;
	double distab, lmax, seglen, EdeposInseg;

	std::cout<<"Segments length (in cm): ";
	std::cin>>lmax;
	std::cout<<std::endl;

	for (i=0; i<t->GetEntries(); i++){
		spX1.clear();
		spY1.clear();
		spZ1.clear();
		spPlane1.clear();
		spInteg1.clear();
		segX.clear();
		segY.clear();
		segZ.clear();
		segdxdz.clear();
		segdydz.clear();
		segIncz.clear();
		segXdir.clear();
		segYdir.clear();
		segZdir.clear();
		segDepoE.clear();
		segLength.clear();
		segTrueDepoE.clear();
		EdepoX.clear();
		EdepoY.clear();
		EdepoZ.clear();
		Edepo.clear();
		segXdirErr.clear();
		segYdirErr.clear();
		segZdirErr.clear();
		
		double MuPx, MuPy, MuPz;

		t->GetEntry(i);

		if (spX->size()<30) continue;

		Run1=Run;
		SubRun1=SubRun;
		ID1=ID;

		MuStartX = MCPartsX->at(0);
		MuStartY = MCPartsY->at(0);
		MuStartZ = MCPartsZ->at(0);

		MuE = MCPartsE->at(0)*1000.;
		MuPx = MCPartsPx->at(0)*1000.;
		MuPy = MCPartsPy->at(0)*1000.;
		MuPz = MCPartsPz->at(0)*1000.;
		MuP = TMath::Sqrt( MuPx* MuPx +  MuPy* MuPy +  MuPz* MuPz);

		MuDirX = MuPx/MuP;
		MuDirY = MuPy/MuP;
		MuDirZ = MuPz/MuP;

		MuDepoE = 0.;
		for (j=0; j<G4EDepos->size(); j++){
			if (G4EDeposPDG->at(j)==13) MuDepoE = MuDepoE + G4EDepos->at(j);
			Edepo.push_back(G4EDepos->at(j));
			EdepoX.push_back(G4EDeposMidX->at(j));
			EdepoY.push_back(G4EDeposMidY->at(j));
			EdepoZ.push_back(G4EDeposMidZ->at(j));
		}

		IsEscaping = 0;

		for(j=0; j<spX->size(); j++){
			spX1.push_back(spX->at(j));
			spY1.push_back(spY->at(j));
			spZ1.push_back(spZ->at(j));
			spPlane1.push_back(spHitsPlane->at(j));
			spInteg1.push_back(spHitsInteg->at(j));

			if (TMath::Abs(spX->at(j))>195 || TMath::Abs(spY->at(j))>195 || spZ->at(j)<5 || spZ->at(j)>495) IsEscaping=1;
		}


		a=0;
    	It=0;
    	do{
      		It++;
      		bFound=false;
      		k=a+1;
          
      		do{

    			distab=TMath::Sqrt((spX->at(a) - spX->at(k))*(spX->at(a) - spX->at(k)) + 
    				(spY->at(a) - spY->at(k))*(spY->at(a) - spY->at(k)) + 
    				(spZ->at(a) - spZ->at(k))*(spZ->at(a) - spZ->at(k)));

        		if (distab>=lmax){
          			b = k;
          			bFound = true;		
    	   		}
        	 		          	
        		k++;

      		}while (bFound==false && k<spX->size());
     

	      	if (bFound==true){

	      		E=0.;
	      		double incrX, incrY,incrZ;
	      		incrX = spX->at(b) - spX->at(a);
	      		incrY = spY->at(b) - spY->at(a);
	      		incrZ = spZ->at(b) - spZ->at(a);

	      		seglen = TMath::Sqrt(incrX*incrX + incrY*incrY + incrZ*incrZ);
	      		double errseglen = 2*TMath::Sqrt(0.01*incrX*incrX + 0.09*(incrY*incrY + incrZ*incrZ))/seglen;

	      		double errDirX, errDirY, errDirZ;
	      		errDirX = TMath::Sqrt(0.01 + incrX*incrX*errseglen*errseglen/(seglen*seglen))/seglen;
	      		errDirY = TMath::Sqrt(0.09 + incrY*incrY*errseglen*errseglen/(seglen*seglen))/seglen;
	      		errDirZ = TMath::Sqrt(0.09 + incrZ*incrZ*errseglen*errseglen/(seglen*seglen))/seglen;

	      		double pointsXYZ[3];
				TPrincipal *princ=new TPrincipal(3, "");

				for (m=a; m<=b; m++){

					pointsXYZ[0]=spX->at(m)-spX->at(a);
					pointsXYZ[1]=spY->at(m)-spY->at(a);
					pointsXYZ[2]=spZ->at(m)-spZ->at(a);

					princ->AddRow(pointsXYZ);

					if (spHitsPlane->at(m)==2) E+=0.00191*spHitsInteg->at(m);
				}

				princ->MakePrincipals();

				const TMatrixD *eigenMatrix = princ->GetEigenVectors();

				xdirAux=eigenMatrix->operator()(0,0);
				ydirAux=eigenMatrix->operator()(1,0);
				zdirAux=eigenMatrix->operator()(2,0);

				princ->Clear();

				double norm = TMath::Sqrt(xdirAux*xdirAux + ydirAux*ydirAux + zdirAux*zdirAux);

				xdirAux = xdirAux/norm;
				ydirAux = ydirAux/norm;
				zdirAux = zdirAux/norm;
				
				if (xdirAux*incrX<0) xdirAux=(-1)*xdirAux;
				if (ydirAux*incrY<0) ydirAux=(-1)*ydirAux;
				if (zdirAux*incrZ<0) zdirAux=(-1)*zdirAux;

				segX.push_back(spX->at(a));
				segY.push_back(spY->at(a));
				segZ.push_back(spZ->at(a));
				segXdir.push_back(xdirAux);
				segYdir.push_back(ydirAux);
				segZdir.push_back(zdirAux);
				segXdirErr.push_back(errDirX);
				segYdirErr.push_back(errDirY);
				segZdirErr.push_back(errDirZ);
				segdxdz.push_back(xdirAux/zdirAux);
				segdydz.push_back(ydirAux/zdirAux);
				segIncz.push_back(spZ->at(b)-spZ->at(a));
				segDepoE.push_back(E);
				segLength.push_back(seglen);



	      		//TRUE ENERGY DEPO IN THIS SEGMENT·····································
	      		
	      		if (spX->at(a) > spX->at(b)){
	      			seglim1x=spX->at(b);
	      			seglim2x=spX->at(a);
	      		}
	      		if (spX->at(b) > spX->at(a)){
	      			seglim1x=spX->at(a);
	      			seglim2x=spX->at(b);
	      		}

	      		if (spY->at(a) > spY->at(b)){
	      			seglim1y=spY->at(b);
	      			seglim2y=spY->at(a);
	      		}
	      		if (spY->at(b) > spY->at(a)){
	      			seglim1y=spY->at(a);
	      			seglim2y=spY->at(b);
	      		}

	      		if (spZ->at(a) > spZ->at(b)){
	      			seglim1z=spZ->at(b);
	      			seglim2z=spZ->at(a);
	      		}
	      		if (spZ->at(b) > spZ->at(a)){
	      			seglim1z=spZ->at(a);
	      			seglim2z=spZ->at(b);
	      		}

	      		EdeposInseg = 0.;
	      		for (m=0; m<G4EDepos->size(); m++){
	      			if (G4EDeposMidX->at(m)>=seglim1x && G4EDeposMidX->at(m)<seglim2x &&
	      				G4EDeposMidY->at(m)>=seglim1y && G4EDeposMidY->at(m)<seglim2y &&
	      				G4EDeposMidZ->at(m)>=seglim1z && G4EDeposMidZ->at(m)<seglim2z &&
	      				(G4EDeposPDG->at(m)==13 || TMath::Abs(G4EDeposPDG->at(m))==11)){

	      				EdeposInseg = EdeposInseg + G4EDepos->at(m);
	      			}
	      		}

	      		segTrueDepoE.push_back(EdeposInseg);
	      		
	      		//..................................................................

	       	}//end segment-found condition

	       	a=b;

	    }while (b<spX->size()-1 && k<spX->size());//end loop for segments definition 

	    tSave->Fill();	    
	}//end entries loop


	tSave->Write();
	fSave->Close();
	f->Close();






}