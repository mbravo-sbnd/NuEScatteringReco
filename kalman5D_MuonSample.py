""" 
Macro that computes the Kalman Filter process using 5-D state vector. Using the information in  SegmentedMuons.root 
ROOT file.
See Friday-meetings slides for more informaion about this Kalman Filter application.

Objects that are involved in Kalman Filter process:
- xkk_list = list where all the state vectors (stVect) made in a filter step are stored
- xkk_1_list = list where all the state vectors (stVect) made in a prediction step are stored
- Ckk_list = list where all the covariance matrices (CovMat) made in a filter step are stored
- Ckk_1_list = list where all the covariance matrices (CovMat) made in a prediction step are stored
- mk_list = python list object where all measurement vectors are stored
- Vk_list = python list object where all covariance matrices of measurement vectors (Vk) are stored
- Fk_list = python list object where all prediction matrices (FAux) are stored
- FkT_list = python list object where all transposed prediction matrices (FAuxT) are stored
- Qk_list = python list object where all covariance matrices of prediction matrices (FAux) (Qk) are stored
- K = Kalman gain matrix
- KT = transposed Kalman gain matrix
- H = H matrix 
- A = smoothing matrix

Output: FinalStateVectors.txt text file with the following information for each event:
- Run (int) = Run number of the event 
- SubRun (int) = SubRun number of the event
- ID (int) = ID number of the event
- IsEscaping (int) => 0 = it is a contained muon; 1=>it is an escaping muon
- MuP (double) = true initial muon momentum (MeV)
- MuStartX  (double)= true X component of the initial muon start point (cm, SBND coordinate system)
- MuStartY (double) = true Y component of the initial muon start point (cm, SBND coordinate system)
- MuStartZ (double) = true Z component of the initial muon start point (cm, SBND coordinate system)
- MuDirX/MuDirZ (double) = true initial dx/dz of the muon direction
- MuDirY/MuDirZ (double) = true initial dy/dz of the muon direction
- pestimation = reconstructed initial mometum obtained with first component of the final state vector
- 2nd, 3rd, 4th and 5th component of the state vector (stVect)

"""

from ROOT import *
import numpy as np
import math
import random

f = TFile("SegmentedMuons.root")
t = f.Get("Tree")

text1=open("FinalStateVectors.txt", "w")

mk_list=[]
Vk_list=[]
xkk_list=[]
xkk_1_list=[]
Ckk_list=[]
Ckk_1_list=[]
Fk_list=[]
FkT_list=[]
Qk_list=[]

for event in t:
	
	mk_list.clear()
	Vk_list.clear()
	xkk_list.clear()
	xkk_1_list.clear()
	Ckk_list.clear()
	Ckk_1_list.clear()
	Fk_list.clear()
	FkT_list.clear()
	Qk_list.clear()


	Run = event.__getattr__("Run")
	SubRun = event.__getattr__("SubRun")
	ID = event.__getattr__("ID")

	IsEscaping = event.__getattr__("IsEscaping")

	MuDirX = event.__getattr__("MuDirX")
	MuDirY = event.__getattr__("MuDirY")
	MuDirZ = event.__getattr__("MuDirZ")
	MuE = event.__getattr__("MuE")
	MuP = event.__getattr__("MuP")
	
	MuDepoE = event.__getattr__("MuDepoE")

	MuStartX = event.__getattr__("MuStartX")
	MuStartY = event.__getattr__("MuStartY")
	MuStartZ = event.__getattr__("MuStartZ")

	segX = event.__getattr__("segX")
	segY = event.__getattr__("segY")
	segZ = event.__getattr__("segZ")
	segdxdz = event.__getattr__("segdxdz")
	segdydz = event.__getattr__("segdydz")
	segIncz = event.__getattr__("segIncz")
	segXdir = event.__getattr__("segXdir")
	segYdir = event.__getattr__("segYdir")
	segZdir = event.__getattr__("segZdir")
	segDepoE = event.__getattr__("segDepoE")
	segLength = event.__getattr__("segLength")
	segXdirErr = event.__getattr__("segXdirErr")
	segYdirErr = event.__getattr__("segYdirErr")
	segZdirErr = event.__getattr__("segZdirErr")
	segTrueDepoE = event.__getattr__("segTrueDepoE")


	nsegs = len(segX)

	truep = MuP
	trueE = MuE

	text1.write(str(Run))
	text1.write("   ")
	text1.write(str(SubRun))
	text1.write("   ")
	text1.write(str(ID))
	text1.write("   ")
	text1.write(str(IsEscaping))
	text1.write("\n")
	text1.write(str(MuP))
	text1.write("   ")
	text1.write(str(MuStartX))
	text1.write("   ")
	text1.write(str(MuStartY))
	text1.write("   ")
	text1.write(str(MuStartZ))
	text1.write("   ")
	text1.write(str(1.0*MuDirX/MuDirZ))
	text1.write("   ")
	text1.write(str(1.0*MuDirY/MuDirZ))
	text1.write("\n")

	m0 = np.zeros(shape=(5,1))
	mk_list = [m0]

	m0 = np.zeros(shape=(5,5))
	Vk_list = [m0]

	#Obtaining all measurement vectors of the event and storing them in their list

	for j in range(1, nsegs):

		costheta = segXdir[j-1]*segXdir[j] + segYdir[j-1]*segYdir[j] + segZdir[j-1]*segZdir[j]

		theta = math.acos(costheta)

		errthetaAux = segXdir[j]*segXdir[j]*segXdirErr[j-1]*segXdirErr[j-1] + segXdir[j-1]*segXdir[j-1]*segXdirErr[j]*segXdirErr[j]
		errthetaAux = errthetaAux + segYdir[j]*segYdir[j]*segYdirErr[j-1]*segYdirErr[j-1] + segYdir[j-1]*segYdir[j-1]*segYdirErr[j]*segYdirErr[j]
		errthetaAux = errthetaAux + segZdir[j]*segZdir[j]*segZdirErr[j-1]*segZdirErr[j-1] + segZdir[j-1]*segZdir[j-1]*segZdirErr[j]*segZdirErr[j]

		errtheta = np.sqrt(errthetaAux/(1. - costheta*costheta))

		m0 = np.matrix([[theta], [segX[j]], [segY[j]], [segdxdz[j]], [segdydz[j]]])
		mk_list.append(m0)

		
		Vk = np.matrix([[errtheta, 0, 0, 0, 0],[0, 0.1, 0, 0, 0],[0, 0, 0.3, 0, 0.],[0, 0, 0, 0.4, 0],[0, 0, 0, 0, 0.4]])

		Vk_list.append(Vk)



	iden = np.matrix([[1,0, 0, 0, 0],[0,1, 0, 0, 0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1]])

	#Initializing state vector and covariance matrix
	#Due to difference in number of entries, some lists are initialized with 0 values

	totalDepoE = 0.
	for j in range(nsegs):
		totalDepoE+=segDepoE[j]

	estE = (totalDepoE + 106.)
	estP = np.sqrt(totalDepoE*totalDepoE + 212*totalDepoE)

	stVect = np.matrix([[estE/(estP*estP)], [0.], [0.], [0.], [0.]])
	xkk_list.append(stVect)


	m0 = np.zeros(shape=(5,1))
	xkk_1_list.append(m0)
	
	CovMat = np.matrix([[100000., 0, 0, 0, 0],[0, 100000., 0, 0, 0],[0, 0, 100000., 0, 0],[0, 0, 0, 100000., 0.],[0, 0, 0, 0, 100000.]])
	Ckk_list.append(CovMat)

	m0 = np.ones(shape=(5,5))
	Ckk_1_list.append(m0)

	#Begining of kalman filter process

	itcounter = 0

	while itcounter<25:

		#Variables used to store true information about the energy and momentum in each begining of the 
		#segment, in case it is wanted to be compared

		truep = MuP
		trueE = MuE

		#FORWARD PROCESS*************************************

		for j in range(1, nsegs):

			#True values update
			trueE = trueE - segTrueDepoE[j-1]
			truep = np.sqrt(trueE*trueE - 11236.)

			#PREDICTION PROCESS*****************************

			#Estimate of muon momentum and E in this step according to qst element of the state vector
			pestimation = (0.7071/(stVect[0,0]))*np.sqrt(1 + np.sqrt(1 + 44944*stVect[0,0]*stVect[0,0]))
			Eestimation = np.sqrt(pestimation*pestimation + 11236.)

			#Predicion matrix obtention
			alpha = segDepoE[j-1]/Eestimation
			if segDepoE[j-1]<0.1:
				alpha = 0.1/Eestimation

			alphaP = alpha*Eestimation/pestimation
			erralpha = 0.25*segDepoE[j-1]/(pestimation*(1-alphaP)*(1-alphaP))

			FAux = np.matrix([[(1-alpha)/((1-alphaP)*(1-alphaP)), 0, 0, 0, 0], [0,1,0,segIncz[j-1],0], [0,0,1,0,segIncz[j-1]], [0,0,0,1,0], [0,0,0,0,1]])
			Qk = np.matrix([[erralpha, 0, 0, 0, 0],[0, np.absolute(0.3*segdxdz[j-1]), 0, 0, 0],[0, 0, np.absolute(0.3*segdydz[j-1]), 0, 0],[0, 0, 0, 0.01, 0],[0, 0, 0, 0, 0.01]])
			
			Fk_list.append(FAux)
			Qk_list.append(Qk)

			FAuxT = np.transpose(FAux)
			FkT_list.append(FAuxT)

			#Extrapolation/prediction of state vector and covariance matrix
			stVect = np.dot(FAux, stVect)
			xkk_1_list.append(stVect)

			CovMat = np.dot(CovMat, FAuxT)
			CovMat = np.dot(FAux, CovMat)
			CovMat = CovMat + Qk
			Ckk_1_list.append(CovMat)

			#FILTERING PROCESS*****************************

			#H matrix
			Const = 13.6*np.sqrt(segLength[j]/14.)*(1 + 0.038*np.log(segLength[j]/14.))
			H = np.matrix([[Const,0, 0, 0, 0],[0,1, 0, 0, 0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1]])

			#K and KT matrix creation process
			HCH = np.dot(CovMat, H)
			HCH = np.dot(H, HCH)
			VHCH = Vk_list[j] + HCH
			VHCHInv = np.linalg.inv(VHCH)			
			K = np.dot(H, VHCHInv)
			K = np.dot(CovMat, K)
			KT = np.transpose(K)

			#Filtering of state vector and covariance matrix
			mHx = np.dot(H, stVect)
			mHx = mk_list[j] - mHx
			KmHx = np.dot(K, mHx)
			stVect = stVect + KmHx
			xkk_list.append(stVect)

			idKH = np.dot(K, H)
			idKH = iden - idKH
			idKHT = np.transpose(idKH)
			prodAux = np.dot(CovMat, idKHT)
			prodAux = np.dot(idKH, prodAux)
			KVKT = np.dot(Vk_list[j], KT)
			KVKT = np.dot(K, KVKT)
			CovMat = prodAux + KVKT
			Ckk_list.append(CovMat)

		#Due to indices need, filling some list with 0's
		m0 = np.zeros(shape=(5,5))
		Fk_list.append(m0)
		FkT_list.append(m0)
		Qk_list.append(m0)

		

		#BACKWARDS PROCESS*******************************************
		for j in range(nsegs - 2, -1, -1):

			#Smoother gain matrix creation process
			Ckk_1Inv = np.linalg.inv(Ckk_1_list[j+1])
			A = np.dot(FkT_list[j], Ckk_1Inv)
			A = np.dot(Ckk_list[j], A)
			AT = np.transpose(A)

			#Smoothing of state vector and covariance matrix
			stVectAux = stVect - xkk_1_list[j+1]
			stVectAux = np.dot(A, stVectAux)
			stVect = xkk_list[j] + stVectAux

			CovMinusCkk1 = CovMat - Ckk_1_list[j+1]
			CovMinusCkk1 = np.dot(CovMinusCkk1, AT)
			CovMinusCkk1 = np.dot(A, CovMinusCkk1)
			CovMat = Ckk_list[j] + CovMinusCkk1


		#Recontructed initial muon momentum in this forward + backwards iteration (itcounter)
		pestimation = (0.7071/stVect[0,0])*np.sqrt(1 + np.sqrt(1 + 44944*stVect[0,0]*stVect[0,0]))

		#Restoring initial values for the following forward + backwards iteration (itcounter)

		Ckk_list.clear()
		xkk_1_list.clear()
		Ckk_1_list.clear()
		xkk_list.clear()
		Fk_list.clear();
		FkT_list.clear();

		xkk_list.append(stVect)
		Ckk_list.append(CovMat)

		m0 = np.zeros(shape=(5,1))
		xkk_1_list.append(m0)

		m0 = np.zeros(shape=(5,5))
		Ckk_1_list.append(m0)

		itcounter+=1


	#The resulting reconstructed initial muon momentum 
	pestimation = (0.7071/(stVect[0,0]))*np.sqrt(1 + np.sqrt(1 + 44944*stVect[0,0]*stVect[0,0]))

	#text file filling
	text1.write(str(pestimation))
	text1.write("   ")
	text1.write(str(stVect[1,0]))
	text1.write("   ")
	text1.write(str(stVect[2,0]))
	text1.write("   ")
	text1.write(str(stVect[3,0]))
	text1.write("   ")
	text1.write(str(stVect[4,0]))
	text1.write("\n")



text1.close()



