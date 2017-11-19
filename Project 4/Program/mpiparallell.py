from pylab import* 
import numpy as np
import sys, math

from mpi4py import MPI


def periodic(i,limit,add): 
	if i+add == limit:
		return 0
	elif i+add == -1:
		return limit-1
	else:
		return i+add

#periodic(i, limit,add)

EAv= []
temp= []
Mag= []
CV = []
Xi = []

def MonteCarlo( latticeDim, cycles): #calculate the energy and magnetization for a given temp.
	# cycles: MonteCarlo cycles (how many times do we flip the matrix?)
	#latticeDim = dim of square matrix
	# EAverage = energy of matrix averaged over cycles, normalized to spins**2
	#MagAverage= magnetic filed of matrix, averaged over cycles, normalized to spins**2
	#EVariance = variance of energy, normalized
	#MagAbsAverage= absolute value of magnetic field, average over cycles
	#Setup spin matrix, initialize to ground state
	
	#MPI Initializations
	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()


	Tstart= 2.0
	Tend = 2.3
	dt=0.05	

	E_slutt = 0
	M_slutt = 0
	heatCapacity_slutt =0
	Susceptibility_slutt = 0

	#spinMatrix = np.zeros((latticeDim,latticeDim),np.int8) + 1 #ALL SPIN UP

	spinMatrix = np.zeros((latticeDim,latticeDim),np.int8) #RANDOM 

	for i in xrange(latticeDim):
		for j in xrange(latticeDim):
			if np.random.random() < 0.5:
				spinMatrix[i,j] = 1
			else:
				spinMatrix[i,j] = -1

	Trange = linspace(Tstart,Tend,int((Tend-Tstart)/dt))
	for temperature in Trange:
		#create and initialize variables
		StartTime = MPI.Wtime()
		E= M= 0
		EAverage= E2Average =MagAverage= Mag2Average = MagAbsAverage= 0
		k= 1.0
		J = 1.0
		beta = 1/float(k*temperature)
		NumberOfAcceptedStates= 0
		#Possible energy changes, -8J, -4J, 0J, 4J, 8J 
		w= np.zeros(17,np.float64) #17=16 +1, 
		for degeneration in xrange(-8,9,4):
			w[degeneration +8]= math.exp(-degeneration*J*beta) #legger til 8 

		#print w

		#Calculate initial magnetization
		M= spinMatrix.sum()
		for j in xrange(latticeDim):
			for i in xrange(latticeDim):
				E -= spinMatrix.item(i,j)*(spinMatrix.item(periodic(i,latticeDim,-1),j) + spinMatrix.item(i,periodic(j,latticeDim,1))) #initial energy 

		
		NumberOfAcceptedStates = 0
		#start metropolis MonteCarlo Computation
		for i in xrange(cycles): #monte carlo cycle
			#loop over all spins, pick a random spin each time
			for s in xrange(latticeDim**(2)):
				x= int(np.random.random()*latticeDim)
				y= int(np.random.random()*latticeDim)

				spinUp = spinMatrix.item(x,periodic(y,latticeDim,1)) 
				spinDown = spinMatrix.item(x,periodic(y,latticeDim,-1)) 
				spinRight = spinMatrix.item(periodic(x,latticeDim,1),y)
				spinLeft = spinMatrix.item(periodic(x,latticeDim,-1),y)

				deltaE= 2*spinMatrix.item(x,y)*( spinLeft+ spinRight + spinUp + spinDown)


				if np.random.random() <= w[deltaE +8]:
					#accept
					spinMatrix[x,y] = -spinMatrix[x,y]
					M += 2*spinMatrix[x,y] #flipped spin in x and y
					E += deltaE
					NumberOfAcceptedStates +=1


			
			#updating exceptation values
			EAverage += E
			E2Average += E**2
			MagAverage +=M
			Mag2Average += M**2
			MagAbsAverage += math.fabs(M)
	

		#To get the average values
		EAverage /=float(cycles)
		E2Average /=float(cycles)
		MagAverage /=float(cycles)
		Mag2Average /=float(cycles)
		MagAbsAverage /=float(cycles)


		heatCapacity= (E2Average- EAverage**2)/float(latticeDim**(2)*temperature**(2))

		Susceptibility  = (Mag2Average- MagAbsAverage**2)/float(latticeDim**(2)*temperature)

		EAverage /=float(latticeDim**2)
		MagAverage /=float(latticeDim**2)
		MagAbsAverage /=float(latticeDim**2)


		Elist = np.array(EAverage)
		Eslist = np.array(0.)
		Mlist = np.array(MagAbsAverage)
		Mslist = np.array(0.)
		CVlist = np.array(heatCapacity)
		CVslist = np.array(0.)
		Xilist = np.array(Susceptibility)
		Xislist= np.array(0.)

		comm.Reduce(Elist,Eslist, op= MPI.SUM)
		comm.Reduce(Mlist,Mslist, op=MPI.SUM)
		comm.Reduce(CVlist,CVslist, op= MPI.SUM)
		comm.Reduce(Xilist,Xislist, op=MPI.SUM)

		EAv.append(Eslist/size)
		temp.append(temperature)
		Mag.append(Mslist/size)
		CV.append(CVslist/size)
		Xi.append(Xislist/size)

		if rank == 0:
			print Eslist/size, Mslist/size, CVslist/size, Xislist/size, temperature


		EndTime = MPI.Wtime()

		Totaltime = EndTime- StartTime

		if (rank ==0):
			print ####
			print Totaltime
			



	return  Eslist/size, CVslist/size, Mslist/size, Xislist/size, temperature

def PhaseTransisions():
	# to run: -n 2 python mpi_parallell.py 
	cycles= 10**5

	latticeDim= 40

	EAverage, EVariance, MagAbsAverage, MagVariance,temperature = MonteCarlo(latticeDim,cycles)

	EAvv = np.array(EAv)
	ttemp= np.array(temp)
	MM= np.array(Mag)
	CVV= np.array(CV)
	Xii = np.array(Xi)

	figure()
	plt.plot(ttemp, EAvv*(1./float(latticeDim)/float(latticeDim)),'r') #eaverage per spin
	hold('on')
	plt.plot(ttemp, MM*(1./float(latticeDim)/float(latticeDim)),'b') 
	hold('on')
	plt.plot(ttemp, CVV*(1./float(latticeDim)/float(latticeDim)),'g') 
	hold('on')
	plt.plot(ttemp, Xii*(1./float(latticeDim)/float(latticeDim)),'brown') 

	plt.xlabel('k_BT')
	plt.ylabel('$ <E>/J $')
	plt.title(' $%d\\times%d$ Ising model with $T = %s$ and $MCS = %s$ ' % (latticeDim,latticeDim, temperature, cycles))
	plt.legend ([r"E/J","|M|/ J","$C_V$","$\chi$"], loc = 'best')
	plt.savefig('energyising.pdf',bbox_inches='tight')
	plt.show()

PhaseTransisions()


