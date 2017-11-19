from pylab import* 
import numpy as np
import sys, math




#choose correct matrix with periodic boundary conditions
def periodic(i,limit,add): 
	return (i+limit+add) % limit

#periodic(i, limit,add)

Ehist, MagAbsHist, E2, histogram= [], [], [], []


def MonteCarlo(temperature, latticeDim, cycles): 
	#calculate the energy and magnetization for a given temp.
	# cycles: MonteCarlo cycles (how many times we flip the matrix)
	#latticeDim = dimention of square matrix
	# EAverage = energy of matrix averaged over cycles
	#MagAverage= magnetic filed of matrix, averaged over cycles
	#EVariance = variance of energy
	#MagAbsAverage= absolute value of magnetic field



	#set up spin matrix: ONE FOR RANDOM AND ONE WITH ALL SPIN UP
	#spinMatrix = np.zeros((latticeDim,latticeDim),np.int8) + 1 #ALL SPIN UP

	spinMatrix = np.zeros((latticeDim,latticeDim),np.int8) #RANDOM 

	for i in xrange(latticeDim):
		for j in xrange(latticeDim):
			if np.random.random() < 0.5:
				spinMatrix[i,j] = 1
			else:

				spinMatrix[i,j] = -1   




	#create and initialize variables
	E= M= 0.
	EAverage= E2Average =MagAverage= Mag2Average = MagAbsAverage= 0
	k= 1.0
	J = 1.0
	beta = 1/float(k*temperature)
	NumberOfAcceptedStates= 0.
	#Possible energy changes, -8J, -4J, 0J, 4J, 8J 
	w= np.zeros(17,np.float64) #17=16 +1, 
	for degeneration in xrange(-8,9,4):
		w[degeneration +8]= math.exp(-degeneration*J*beta) #add 8 

	#print w

	#Calculate initial magnetization and energy
	M= spinMatrix.sum()
	for j in xrange(latticeDim):
		for i in xrange(latticeDim):
			E -= spinMatrix.item(i,j)*(spinMatrix.item(periodic(i,latticeDim,-1),j) + spinMatrix.item(i,periodic(j,latticeDim,1))) #initial energy 

	
	NumberOfAcceptedStates = 0.
	#start metropolis MonteCarlo Computation
	for i in xrange(cycles): 
		#metropolis algorithm
		#loop over all spins, and pick a random spin each time
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
		Ehist.append(EAverage/float(i+1))
		MagAbsHist.append(MagAbsAverage/float(i+1))
		E2.append(E**2)
		histogram.append(E)


	# average values
	EAverage /=float(cycles)
	E2Average /=float(cycles)
	MagAverage /=float(cycles)
	Mag2Average /=float(cycles)
	MagAbsAverage /=float(cycles)

	#Calculate variance and normalize 
	EVariance = (E2Average- EAverage**2)/float(latticeDim**(2) * temperature**(2))
	print " heat= ", EVariance

	MagVariance = (Mag2Average- MagAbsAverage**2)/float(latticeDim**(2)*temperature)
	print "susceptibility= ", MagVariance

	#Normalize averages
	EAverage /=float(latticeDim**2)
	MagAverage /=float(latticeDim**2)
	MagAbsAverage /=float(latticeDim**2)
	print " EAverage = ", EAverage
	print "MagAbsAverage=", MagAbsAverage

	print "acceptable:", NumberOfAcceptedStates/float(latticeDim**2)



	return EAverage, EVariance, MagAverage, MagVariance, MagAbsAverage, E



def IsingModel2x2():
	"""
	Ising model for 2x2 lattice. 
	"""
	temperature= 1.0
	latticeDim= 2
	cycles= 10**(5)

	EAverage, EVariance, MagAverage, MagVariance, MagAbsAverage, E = MonteCarlo(temperature, latticeDim, cycles)

#IsingModel2x2()


def Analytic():

	"""
	Anlytical values for partition function, expectation values for energy, mean absolute value of magnetization,
	specific heat and susceptibility. All divided by latticesize**(2). 
	"""

	N= 4 

	z= 4*math.cosh(8)+ 12

	E_forventa =-32*math.sinh(8)/float(z*N)

	E2_forventa = 256*math.cosh(8)/float(z*N)


	Mabs_forventa = (8*math.exp(8) +16)/float(z*N)

	M2abs_forventa = (32*math.exp(8) +32)/float(z*N)

	Cv = (E2_forventa*N - (E_forventa*N)**(2))/float(N)

	X= (M2abs_forventa*N -(Mabs_forventa*N)**(2))/float(N)

	print ######
	print "ANALYTIC:"
	print "Z= ", z, "E_forventa:", E_forventa, "Mabs_forventa:", Mabs_forventa, "Cv:", Cv, "X:", X
	return z, E_forventa, Mabs_forventa, Cv, X

#z, E_forventa, Mabs_forventa, Cv, X= Analytic()



def IsingModel20x20ExpectationValues():
	"""
	Ising model for 20 x 20 latticesize. 
	Plotting mean energy and mean magnetisation as a functions of the number of Monte Carlo cycles for temperature 2.4 and 1.

	"""
	temperature= 2.4
	latticeDim= 20	
	cycles= 10**(4)


	EAverage, EVariance, MagAverage, MagVariance, MagAbsAverage, E= MonteCarlo(temperature,latticeDim,cycles)

	EEhist = np.array(Ehist)
	MMagAbsHist = np.array(MagAbsHist)


	
	figure()
	#plt.plot(range(cycles), EEhist*(1./float(latticeDim)/float(latticeDim)),'r')
	plt.plot(range(cycles),MMagAbsHist*(1./float(latticeDim)/float(latticeDim)),'b')

	
	#plt.axis([-0.5, 10000,-0.5, 1.5])
	plt.xlabel("Monte Carlo Cycles")
	plt.ylabel("$<|M|/N>$")
	#plt.ylabel("$<E/N>$")
	plt.title('All spin configurations up for 20 x 20 lattice with T = %s and MCS =%s ' % (temperature, cycles))
	#plt.savefig('20x20energyspinuptemp24.pdf',bbox_inches='tight')
	#plt.savefig('spinupmagtemp1.pdf')
	#plt.tight_layout()
	plt.show()


#IsingModel20x20ExpectationValues()
	

def IsingModel20x20AcceptedConfiguration():

	"""
	Ising model for 20x20- lattice.
	Studying how the number of accepted configuration behave as a funtion of temperature T.
	"""
	T = [1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4] #temperature

	#accepted config. with random spin orientation as starting config.
	acceptRandom = [12.3017, 33.0117, 91.1792, 210.5983, 424.4942, 795.6292, 1486.4667, 2697.7233] #accepted config. with random spin orientation as starting config.
	#accepted config. with ordered spin orientation as starting config.
	acceptOrdered= [11.1433, 32.9683, 91.7609, 215.555, 418.0967, 814.4617, 1527.8092, 2680.7783]

	figure()
	plt.plot(T,acceptRandom, 'b')
	plt.title('Number of accepted configurations as a function of temperature')
	plt.xlabel('$k_BT$')
	plt.ylabel('Accepted configurations')
	#plt.savefig('acceptedconfig.pdf')
	plt.show()
	
#IsingModel20x20AcceptedConfiguration()

def IsingModel20x20ProbabilityDistribution():

	"""
	Ising Model for 20x 20 lattice.
	Analyzing the probability distribution for temperature 1.0 and 2.4
	"""

	temperature= 2.4
	latticeDim= 20	
	cycles= 10**(4)


	EAverage, EVariance, MagAverage, MagVariance, MagAbsAverage, E= MonteCarlo(temperature,latticeDim,cycles)

	HHistogram= np.array(histogram)

	plt.hist(HHistogram, bins= 25, normed = 1)
	plt.title('Probability distribution for $k_BT$= %s' % temperature)
	plt.xlabel('Energy')
	plt.ylabel('Probability distribution')
	#plt.axis([-800, -200, 0.0, 0.05]) #T=1.0
	#plt.axis([-700, -100, 0.0, 0.012]) #T=2.4
	#plt.savefig('Ehist24random.pdf')
	plt.show()

	#EvarNY= sqrt(np.sum(E2)/float(len(E2)) - (np.sum(histogram)/float(len(histogram)))**2 )

	#print "Evar: ", EvarNY



#IsingModel20x20ProbabilityDistribution()




