""" 
Macro that computes the reconstructed direction using max Theta-RMS cut and Kalman Filter for 
direction reconstruction step.
Input:SegmentAndPortions.root ROOT file and number of points per segment for Kalman Filter 
algorithm (nPointsPerSeg), entered by the user. 
Output: ROOT histogram where angles between reco and true electron direction are stored
(thetaRT variable, in degrees)

"""



from ROOT import *
import numpy as np
import math
import pandas as pd
from sklearn.decomposition import PCA
from sklearn import preprocessing

#classes definition to translate info from one function to another
class SelecPoints:
	def __init__(self, selecPointsX, selecPointsY, selecPointsZ):
		#Class made of 3 vectors, each one with the X, Y and Z position of each one of the selected points (in cm, SBND coordinate system)
		self.X = selecPointsX
		self.Y = selecPointsY
		self.Z = selecPointsZ

class Segments:
	def __init__(self, validSet, segX, segY, segZ, segIncZ, segXdir, segYdir, segZdir, segXdirErr, segYdirErr, segZdirErr, segdxdz, segdydz):
		#Class containing the following information:
		#validSet: integer that tells if the set of segments is valid (if so, it returns 1, otherwise, it returns 0)
		#segX, segY, segZ: vectors with the corresponding X, Y, and Z positions of the start point of each segment (in cm, SBND coordinate system)
		#segIncZ: vector with z_finalPoint - z_InitialPoint of each segment (in cm)
		#segXdir, segYdir, segZdir: vectors with the corresponding X, Y, and Z components of the direction of each segment
		#segXdirErr, segYdirErr, segZdirErr: vectors with the corresponding direction estimated errors in X, Y and Z component of each segment
		#segdxdz, segdydz: X/Z and Y/Z direction components of each segment
		self.validSet = validSet
		self.X = segX
		self.Y = segY
		self.Z = segZ
		self.incZ = segIncZ
		self.Xdir = segXdir
		self.Ydir = segYdir
		self.Zdir = segZdir
		self.XdirErr = segXdirErr
		self.YdirErr = segYdirErr
		self.ZdirErr = segZdirErr
		self.dxdz = segdxdz
		self.dydz = segdydz
		
class FinalStVect:
	def __init__(self, iniRecoX, iniRecoY, iniRecodxdz, iniRecodydz):
		#Class containing the following information about the final reconstructed state vector.
		#iniRecoX, iniRecoY : reconstructed X and Y initial positions after all Kalman Filter iterations (in cm, SBND coorrdinate system)
		#iniRecodxdz, iniRecodydz: reconstructed X/Z and Y/Z initial direction components after all Kalman Filter iterations
		self.X = iniRecoX 
		self.Y = iniRecoY 
		self.dxdz = iniRecodxdz 
		self.dydz = iniRecodydz

#...............................................................................



#Function definitions
def HitSelection(PointsX, PointsY, PointsZ, StdDevThetaWithPrim, nLastPoint):
	#function that returns the selected points according to theta-RMS cut (in a SelecPoints object)
	maxDeriv = 0.
	nlastsegment = len(StdDevThetaWithPrim) - 1

	for j in range (2, len(StdDevThetaWithPrim)):
		Deriv = (StdDevThetaWithPrim[j+1] - StdDevThetaWithPrim[j])*2./(StdDevThetaWithPrim[j]+StdDevThetaWithPrim[j+1])

		if Deriv > MaxDeriv:
			MaxDeriv = Deriv
			nlastsegment = j

	nPointMax = nLastPoint[nlastsegment]

	#Creating and filling the vectors of the selected points
	selecPointsX = []
	selecPointsY = []
	selecPointsZ = []

	for j in range(nPointMax+1):
		selecPointsX.append(PointsX[j])
		selecPointsY.append(PointsY[j])
		selecPointsZ.append(PointsZ[j])

	#Creating a SelecPoints object with the filled vectors
	sPoints = SelecPoints(selecPointsX, selecPointsY, selecPointsZ)

	return sPoints

#**************************************************************

def KalFiltSegmentation (sPoints, nsegsKF):
	#function that, using the selected points and the required number of segments, returns the resulting segments for Kalman Filter implementation (in a Segments object)
	selecPointsX = sPoints.X
	selecPointsY = sPoints.Y
	selecPointsZ = sPoints.Z

	segX = []
	segY = []
	segZ = []
	segIncZ = []
	segXdir = []
	segYdir = []
	segZdir = []
	segdxdz = []
	segdydz = []
	segXdirErr = []
	segYdirErr = []
	segZdirErr = []

	#List where the points of each segmen are stored in order to apply PCA method to get segment directions
	pointList = []

	nPointsPerSeg = int(len(selecPointsX)/nsegsKF)

	validSet = 1

	for j in range(nsegsKF):
		#using selected points information, each segment is defined and its required information is stored in the corresponding vectors
		a = j * nPointsPerSeg
		b = (j+1) * nPointsPerSeg - 1

		incrX = selecPointsX[b] - selecPointsX[a]
		incrY = selecPointsY[b] - selecPointsY[a]
		incrZ = selecPointsZ[b] - selecPointsZ[a]

		seglength = np.sqrt(incrX*incrX + incrY*incrY + incrZ*incrZ)
		errseglen = 2*np.sqrt(0.01*incrX*incrX + 0.09*(incrY*incrY + incrZ*incrZ))/seglength

		errDirX = np.sqrt(0.01 + incrX*incrX*errseglen*errseglen/(seglength*seglength))
		errDirX = errDirX/seglength
		errDirY = np.sqrt(0.09 + incrY*incrY*errseglen*errseglen/(seglength*seglength))
		errDirY = errDirY/seglength
		errDirZ = np.sqrt(0.09 + incrZ*incrZ*errseglen*errseglen/(seglength*seglength))
		errDirZ = errDirZ/seglength

		for k in range(a, b+1):
			pointList.append([selecPointsX[k], selecPointsY[k], selecPointsZ[k]])

		if len(pointList) < 1:
			print("WARNING!!!!!! NO SELECTED HITS IN ", j, "SEGMENT \n")

		#Defininf PCA object and filling it with the points of the segment
		pca = PCA(n_components = 3)
		pca.fit(pointList)

		#Directionn obtaining
		DirX = pca.components_[2,0]
		DirY = pca.components_[2,1]
		DirZ = pca.components_[2,2]

		pointList.clear()

		#if the segment has null Z direction component, the set is not valid
		if DirZ == 0:
			validSet = 0
			break

		Dirnorm = np.sqrt(DirX*DirX + DirY*DirY + DirZ*DirZ)

		DirX = DirX/Dirnorm
		DirY = DirX/Dirnorm
		DirX = DirX/Dirnorm

		if incrX*DirX < 0:
			DirX = (-1.)*DirX

		if incrY*DirY < 0:
			DirY = (-1.)*DirY

		if incrZ*DirZ < 0:
			DirZ = (-1.)*DirZ

		segX.append(selecPointsX[a])
		segY.append(selecPointsY[a])
		segZ.append(selecPointsZ[a])
		segIncZ.append(incrZ)
		segXdir.append(DirX)
		segYdir.append(DirY)
		segZdir.append(DirZ)
		segdxdz.append(DirX/DirZ)
		segdydz.append(DirY/DirZ)
		segXdirErr.append(errDirX)
		segYdirErr.append(errDirY)
		segZdirErr.append(errDirZ)

	#Creating a Segments object with the filled vectors
	KFSegments = Segments(validSet, segX, segY, segZ, segIncZ, segXdir, segYdir, segZdir, segdxdz, segdydz, segXdirErr, segYdirErr, segZdirErr)

	return KFSegments

#**************************************************************	    

def KalmanFilter(KFSegments):

	#Function that, using the information of the created segments with the selected points, runs the Kalman Filter algorith and returns the reconstructed initial state vector information.
	#See kalman4D_MuonSample.py for more information. Here, elements of Qk matrices are otained with more detail
	segX = KFSegments.X
	segY = KFSegments.Y
	segZ = KFSegments.Z 
	segIncz = KFSegments.incZ   
	segXdir = KFSegments.Xdir
	segYdir = KFSegments.Ydir
	segZdir = KFSegments.Zdir
	segdxdz = KFSegments.dxdz
	segdydz = KFSegments.dydz
	segXdirErr = KFSegments.XdirErr
	segYdirErr = KFSegments.YdirErr
	segZdirErr = KFSegments.ZdirErr

	mk_list=[]
	Vk_list=[]
	xkk_list=[]
	xkk_1_list=[]
	Ckk_list=[]
	Ckk_1_list=[]
	Fk_list=[]
	FkT_list=[]
	Qk_list=[]

	nsegs = len(segX)

	m0 = np.zeros(shape=(4,1))
	mk_list.append(m0)

	m0 = np.zeros(shape=(4,4))
	Vk_list.append(m0)

	for j in range(1, nsegs):

		m0 = np.matrix([[segX[j]], [segY[j]], [segdxdz[j]], [segdydz[j]]])
		mk_list.append(m0)
		
		Vk = np.matrix([[0.1, 0, 0, 0],[0, 0.3, 0, 0.],[0, 0, 0.4, 0],[0, 0, 0, 0.4]])

		Vk_list.append(Vk)
		

	iden = np.matrix([[1, 0, 0, 0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])

	stVect = np.matrix([[1.], [1.], [1.], [1.]])
	xkk_list.append(stVect)

	m0 = np.zeros(shape=(4,1))
	xkk_1_list.append(m0)
	
	CovMat = np.matrix([[100000., 0, 0, 0],[0, 100000., 0, 0],[0, 0, 100000., 0.],[0, 0, 0, 100000.]])

	Ckk_list.append(CovMat)

	m0 = np.zeros(shape=(4,4))
	Ckk_1_list.append(m0)

	totalcounter = 0

	while totalcounter<25:

		#FORWARD PROCESS

		for j in range(1, nsegs):

			FAux = np.matrix([[1,0,segIncz[j-1],0], [0,1,0,segIncz[j-1]], [0,0,1,0], [0,0,0,1]])

			errdxdz = np.sqrt( np.absolute(segXdirErr[j-1]*segXdirErr[j-1]/(segZdir[j-1]*segZdir[j-1]))*np.absolute(segXdirErr[j-1]*segXdirErr[j-1]/(segZdir[j-1]*segZdir[j-1])) + np.absolute(segXdir[j-1]*segZdirErr[j-1]/(segZdir[j-1]*segZdir[j-1]))*np.absolute(segXdir[j-1]*segZdirErr[j-1]/(segZdir[j-1]*segZdir[j-1])))
			errdydz = np.sqrt( np.absolute(segYdirErr[j-1]*segYdirErr[j-1]/(segZdir[j-1]*segZdir[j-1]))*np.absolute(segYdirErr[j-1]*segYdirErr[j-1]/(segZdir[j-1]*segZdir[j-1])) + np.absolute(segYdir[j-1]*segZdirErr[j-1]/(segZdir[j-1]*segZdir[j-1]))*np.absolute(segYdir[j-1]*segZdirErr[j-1]/(segZdir[j-1]*segZdir[j-1])))
			errx = np.sqrt( np.absolute(segIncz[j-1]*errdxdz)*np.absolute(segIncz[j-1]*errdxdz) + np.absolute(segdxdz[j-1]*0.3)*np.absolute(segdxdz[j-1]*0.3) )
			erry = np.sqrt( np.absolute(segIncz[j-1]*errdydz)*np.absolute(segIncz[j-1]*errdydz) + np.absolute(segdydz[j-1]*0.3)*np.absolute(segdydz[j-1]*0.3) )

			Qk = np.matrix([[errx, 0, 0, 0],[0, erry, 0, 0],[0, 0, errdxdz, 0],[0, 0, 0, errdydz]])

			Fk_list.append(FAux)
			Qk_list.append(Qk)

			FAuxT = np.transpose(FAux)
			FkT_list.append(FAuxT)

			#Extrapolation/prediction
			stVect = np.dot(FAux, stVect)
			xkk_1_list.append(stVect)

			CovMat = np.dot(CovMat, FAuxT)
			CovMat = np.dot(FAux, CovMat)
			CovMat = CovMat + Qk
			Ckk_1_list.append(CovMat)

			#H and Kalman gain matrices
			H = iden

			HCH = np.dot(CovMat, H)
			HCH = np.dot(H, HCH)
			VHCH = Vk_list[j] + HCH
			VHCHInv = np.linalg.inv(VHCH)
			K = np.dot(H, VHCHInv)
			K = np.dot(CovMat, K)
			KT = np.transpose(K)

			#Filtering
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

		m0 = np.zeros(shape=(4,4))
		Fk_list.append(m0)
		FkT_list.append(m0)
		Qk_list.append(m0)
		
	
		#BACKWARDS PROCESS
		for j in range(nsegs - 2, -1, -1):

			#Smoother gain matrix
			Ckk_1Inv = np.linalg.inv(Ckk_1_list[j+1])
			A = np.dot(FkT_list[j], Ckk_1Inv)
			A = np.dot(Ckk_list[j], A)
			AT = np.transpose(A)

			#Smoothing
			stVectAux = stVect - xkk_1_list[j+1]
			stVectAux = np.dot(A, stVectAux)
			stVect = xkk_list[j] + stVectAux

			CovMinusCkk1 = CovMat - Ckk_1_list[j+1]
			CovMinusCkk1 = np.dot(CovMinusCkk1, AT)
			CovMinusCkk1 = np.dot(A, CovMinusCkk1)
			CovMat = Ckk_list[j] + CovMinusCkk1

		Ckk_list.clear()
		xkk_1_list.clear()
		Ckk_1_list.clear()
		xkk_list.clear()
		Fk_list.clear()
		FkT_list.clear()

		xkk_list.append(stVect)
		Ckk_list.append(CovMat)

		m0 = np.zeros(shape=(4,1))
		xkk_1_list.append(m0)

		m0 = np.zeros(shape=(4,4))
		Ckk_1_list.append(m0)

		totalcounter+=1

	#Creating the FinalStVect object with the reconstructed information of the initial shower point
	finalrecoVect = FinalStVect(stVect[0,0], stVect[1,0], stVect[2,0], stVect[3,0])

	return finalrecoVect

#***************************************************************************************************************************

def MainFunction():

	f = TFile("SegmentAndPortions.root")
	t = f.Get("Tree")

	#Asking the required points for Kalman Filter segments
	nPointsPerSeg = int(input("Number points per segment in Kalman filter direction reco: "))

	for event in t:

		ID = event.__getattr__("ID")
		SubRun = event.__getattr__("SubRun")
		Run = event.__getattr__("Run")
		ElecE = event.__getattr__("ElecE");
		TrueEDirX = event.__getattr__("TrueEDirX")
		TrueEDirY = event.__getattr__("TrueEDirY")
		TrueEDirZ = event.__getattr__("TrueEDirZ")
		PointsX = event.__getattr__("PointsX")
		PointsY = event.__getattr__("PointsY")
		PointsZ = event.__getattr__("PointsZ")
		nLastHit = event.__getattr__("nLastHit")
		StdDevThetaWithPrim = event.__getattr__("StdDevThetaWithPrim")

		ID = event.__getattr__("ID")
		SubRun = event.__getattr__("SubRun")
		Run = event.__getattr__("Run")
		VertX = event.__getattr__("VertX")
		VertY = event.__getattr__("VertY")
		VertZ = event.__getattr__("VertZ")
		NeutE = event.__getattr__("NeutE")
		ElecE = event.__getattr__("ElecE")
		TrueEDirX = event.__getattr__("TrueEDirX")
		TrueEDirY = event.__getattr__("TrueEDirY")
		TrueEDirZ = event.__getattr__("TrueEDirZ")
		TrueNeutDirX = event.__getattr__("TrueNeutDirX")
		TrueNeutDirY = event.__getattr__("TrueNeutDirY")
		TrueNeutDirZ = event.__getattr__("TrueNeutDirZ")

		PointsX = event.__getattr__("PointsX")
		PointsY = event.__getattr__("PointsY")
		PointsZ = event.__getattr__("PointsZ")

		nseg = event.__getattr__("nSeg")
		segXDir = event.__getattr__("segXDir")
		segYDir = event.__getattr__("segYDir")
		segZDir = event.__getattr__("segZDir")
		nLastPoint = event.__getattr__("nLastPoint")
		ZposLastPoint = event.__getattr__("ZposLastPoint")

		cumulEValue = event.__getattr__("cumulEValue")
		cumulThetaRT = event.__getattr__("cumulThetaRT")
		cumulXDir = event.__getattr__("cumulXDir")
		cumulYDir = event.__getattr__("cumulYDir")
		cumulZDir = event.__getattr__("cumulZDir")

		thetaWithPrim = event.__getattr__("thetaWithPrim")
		meanThetaWithPrim = event.__getattr__("meanThetaWithPrim")
		StdDevThetaWithPrim = event.__getattr__("StdDevThetaWithPrim")

		PandTracksDirX = event.__getattr__("PandTracksDirX")
		PandTracksDirY = event.__getattr__("PandTracksDirY")
		PandTracksDirZ = event.__getattr__("PandTracksDirZ")
		PandShowersDirX = event.__getattr__("PandShowersDirX")
		PandShowersDirY = event.__getattr__("PandShowersDirY")
		PandShowersDirZ = event.__getattr__("PandShowersDirZ")

		#Checking if this is a "reconstructable" event

		if ElecE < 0.05:
			continue

		if len(StdDevThetaWithPrim) < 2:
			continue

		selectedHits = HitSelection(PointsX, PointsY, PointsZ, StdDevThetaWithPrim, nLastPoint)

		#Checking if the sigma theta cut let us work with enough hits for Kalman Filter
		if len(selectedHits.X) < 2*nPointsPerSeg:
			continue

		#...........................................................................

		nKFSegs = int(len(selectedHits.X)/nPointsPerSeg)
		
		KFSegments = KalFiltSegmentation(selectedHits, nKFSegs)

		#Checking if the segment reconstruction for Kalmal Filfer is valid
		if KFSegments.validSet == 0:
			continue

		#.............................................................................

		FinalStVect = KalmanFilter(KFSegments)

		recoX = FinalStVect.X
		recoY = FinalStVect.Y
		recodxdz = FinalStVect.dxdz
		recodydz = FinalStVect.dydz

		recoDirZ = 1/np.sqrt(recodxdz*recodxdz + recodydz*recodydz + 1.)
		recoDirY = recodydz*recoDirZ
		recoDirX = recodxdz*recoDirZ

		thetaRT = math.acos(recoDirX*TrueEDirX + recoDirY*TrueEDirY + recoDirZ*TrueEDirZ)*180./3.14159

		h.Fill(theta)

		Mode = h.GetBinCenter(h.GetMaximumBin());
		Lim1 = h.GetBinLowEdge(h.FindFirstBinAbove(h.GetMaximum()/2));
		Lim2 = h.GetBinLowEdge(h.FindLastBinAbove(h.GetMaximum()/2)+1);
		FWHM = Lim2 - Lim1;

		print("Max thetaRMS derivative reco info: \n")
		print("Number of entries = ", h.GetEntries(), "\n")
		print("Mode value = ", Mode, "  FWHM = ", FWHM, "\n")
		print("Mean value = ", h.GetMean(), "  RMS = ", h.GetRMS(), "\n")


		c = TCanvas("c", "c", 3000, 3000)
		h.Draw()
		c.SaveAs("ThetaRMSCut_KalmanFilterDirReco.pdf")




MainFunction()



