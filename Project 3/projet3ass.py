from pylab import*
import numpy as np
from celestialbody import *
from datetime import datetime
#from bodyclass import *
#from celestialbody.py import celestialbody 


def MassPlanet():
    phys_sun = 1.998*10**(30)
    sun_mass = 1.0
    earth_mass = 6.0*10**(24)/phys_sun

    jupiter_mass = 1.9*10**(27)/phys_sun
    mars_mass = 6.6*10**(23)/float(phys_sun) 
    venus_mass = 4.9*10**(24)/float(phys_sun) 
    saturn_mass = 5.5*10**(26)/float(phys_sun)
    mercury_mass = 3.3*10**(23)/float(phys_sun)
    uranus_mass = 8.8*10**(25)/float(phys_sun) 
    neptun_mass = 1.03*10**(26)/float(phys_sun)

    M= jupiter_mass + sun_mass + earth_mass


    return sun_mass, earth_mass, jupiter_mass, mars_mass, venus_mass, saturn_mass, mercury_mass, uranus_mass, neptun_mass, M

sun_mass, earth_mass, jupiter_mass, mars_mass, venus_mass, saturn_mass, mercury_mass, uranus_mass, neptun_mass, M = MassPlanet()

def values():

	timeFinal = 10 #year #100 #T 
	stepnumber = 100 #1000 #dt #integrationpoints
	dt= (timeFinal)/float(stepnumber) #timestep

	return timeFinal, stepnumber, dt

timeFinal, stepnumber, dt = values()


def EarthSunSystem():
	
	Earth =Solarsystems(position=[1.,0.,0.], velocity=[0.,2.*np.pi,0.], mass=earth_mass)
	Sun = Solarsystems(position=[0.,0.,0.], velocity=[0.,0.,0.], mass=sun_mass)

	bodies = []

	bodies.append(Sun)  # Adds a new element at the end of the vector
	bodies.append(Earth)
	return Earth, Sun, bodies 

#Earth, Sun, bodies = EarthSunSystem()


def EscapeVelocity():
	Earth =Solarsystems(position=[1.,0.,0.], velocity=[0.,8.9,0.], mass=earth_mass) #8.39 , 8.4 

	Sun = Solarsystems(position=[0.,0.,0.], velocity=[0.,0.,0.], mass=sun_mass)
	
	bodies = []

	bodies.append(Sun)  # Adds a new element at the end of the vector
	bodies.append(Earth)
	return Earth, Sun, bodies 

#Earth, Sun , bodies = EscapeVelocity()

def ThreeBodySystem():
	Sun = Solarsystems(position=[0.,0.,0.], velocity=[0.,0.,0.], mass=sun_mass)
	Earth =Solarsystems(position=[1.,0.,0.], velocity=[0.,2.*np.pi,0.], mass=earth_mass)
	Jupiter = Solarsystems(position=[5.4,0.,0.], velocity=[0.,2.*np.pi/float(sqrt(5.4)),0.], mass=jupiter_mass*1000) #5.4 er radiusen til jupiter. 

	bodies = []

	bodies.append(Sun)  # Adds a new element at the end of the vector
	bodies.append(Earth)
	bodies.append(Jupiter)

	return Earth, Sun, Jupiter, bodies

#Earth, Sun, Jupiter, bodies = ThreeBodySystem()

def RealThreebody():

	position_Earth = 1. -(earth_mass+ 5.4*jupiter_mass )/float(M)
	position_Sun =-(earth_mass + jupiter_mass*5.4)/float(M)
	position_Jupiter= 5.4-(earth_mass + jupiter_mass*5.4)/float(M)

	#print position_Jupiter, position_Sun, position_Earth

	#print vo*sun_mass + earth_mass*2*np.pi+jupiter_mass*2*np.pi/sqrt(5.4)

	v0_earth = 2*np.pi
	v0_jupiter = 2*np.pi/sqrt(5.4)
	v0_sun = -(v0_earth*earth_mass + v0_jupiter*jupiter_mass)/float(sun_mass)
	#print v0_earth, v0_jupiter, v0_sun
	#print v0_earth*earth_mass + v0_jupiter*jupiter_mass + v0_sun*sun_mass 

	#REAL THREE BODY
	Earth =Solarsystems(position=[position_Earth,0.,0.], velocity=[0.,v0_earth,0.], mass=earth_mass)
	Sun = Solarsystems(position=[position_Sun,0.,0.], velocity=[0.,v0_sun,0.], mass=1)
	Jupiter= Solarsystems(position=[position_Jupiter,0.,0.], velocity=[0.,v0_jupiter,0.], mass=jupiter_mass)

	bodies = []

	bodies.append(Sun)  # Adds a new element at the end of the vector
	bodies.append(Earth)
	bodies.append(Jupiter)
	return Earth, Sun, Jupiter, bodies

#Earth, Sun, Jupiter , bodies = RealThreebody()

def AllPlanets():
	Sun =Solarsystems(position=[2.249742341967526E-03,5.702259020296952E-03,-1.308868418264601E-04], velocity=[-5.176331796449613E-06, 5.522815922934838E-06,1.208552996089031E-07], mass=1)
	Earth =Solarsystems(position=[9.480957405756560E-01, 3.242804441221033E-01,-1.465820905998789E-04], velocity=[-5.770012343598492E-03, 1.624711177600315E-02,-1.024894917958096E-06],mass=earth_mass)
	Jupiter =Solarsystems(position=[-4.604440597510837,-2.891007523243876,1.149787438139007E-01], velocity=[3.925073055856761E-03,-6.032031115361013E-03,-6.270510562667993E-05], mass= jupiter_mass)
	Mercury= Solarsystems(position=[-3.787765979540608E-01,-1.822820556608564E-01, 1.946364572516068E-02], velocity=[6.639513358754547E-03,-2.399984280022258E-02,-2.571010368861300E-03], mass= mercury_mass)
	Mars = Solarsystems(position=[-1.540775834917753, 6.318374030745555E-01, 5.085869013425003E-02], velocity=[-4.742251167677201E-03,-1.176588347335337E-02,-1.302988318954905E-04], mass=mars_mass)
	Venus = Solarsystems(position=[-5.895447258261073E-01,4.113399941375999E-01,3.958359273166973E-02], velocity=[-1.152142811351218E-02,-1.677859771759771E-02, 4.344523056718858E-04], mass= venus_mass )
	Saturn = Solarsystems(position=[-3.790768464407142E-01,-1.004811933118735E+01,1.897936890208218E-01], velocity=[5.268080367079689E-03,-2.277337042199548E-04,-2.059484097729237E-04],mass= saturn_mass)
	Uranus = Solarsystems(position=[1.786847909037686E+01, 8.793101880397845,-1.988312286911673E-01], velocity=[-1.765372076056564E-03, 3.345607679959963E-03, 3.534922695016412E-05], mass=uranus_mass )
	Neptun = Solarsystems(position=[ 2.860931715814051E+01,-8.836432401622716,-4.773613959400135E-01], velocity=[9.053217800158139E-04, 3.017556437109072E-03,-8.326738981989289E-05],mass=neptun_mass )
	

	bodies = []

	bodies.append(Sun)  # Adds a new element at the end of the vector
	bodies.append(Earth)
	bodies.append(Jupiter)
	bodies.append(Mercury)
	bodies.append(Mars)
	bodies.append(Venus)
	bodies.append(Saturn)
	bodies.append(Uranus)
	bodies.append(Neptun)


	for planet in (bodies):
		for i in range(0,3):
			planet.velocity[i]*= 365.25 #year

	return Sun, Earth, Jupiter, Mercury, Mars, Venus, Saturn, Uranus, Neptun, bodies

#Sun, Earth, Jupiter, Mercury, Mars, Venus, Saturn, Uranus, Neptun, bodies = AllPlanets()


xArr = np.zeros(stepnumber)
yArr = np.zeros(stepnumber)
x1Arr = np.zeros(stepnumber)
y1Arr = np.zeros(stepnumber)

x2Arr = np.zeros(stepnumber)
y2Arr = np.zeros(stepnumber)

x3Arr = np.zeros(stepnumber)
y3Arr = np.zeros(stepnumber)
x4Arr = np.zeros(stepnumber)
y4Arr = np.zeros(stepnumber)
x5Arr = np.zeros(stepnumber)
y5Arr = np.zeros(stepnumber)
x6Arr = np.zeros(stepnumber)
y6Arr = np.zeros(stepnumber)
x7Arr = np.zeros(stepnumber)
y7Arr = np.zeros(stepnumber)

def PlotSunEarth():
#SUN- EARTH PLOT 
		
	for i in xrange(0, stepnumber):
		#EulerCromer(bodies,dt)
		#VerletVelocity(bodies,dt)
		#ForwardEuler(bodies, dt)

		tstart= datetime.now()
		VerletVelocity(bodies,dt)
		tend= datetime.now()

		#print "Tid:", tend- tstart
		#ComputeForcesAndEnergy(bodies)

		xArr[i] = bodies[0].position[0]
		yArr[i] = bodies[0].position[1]

		x1Arr[i] = bodies[1].position[0]
		y1Arr[i] = bodies[1].position[1]

	plt.plot(xArr,yArr,'ro')
	hold('on')
	plt.plot(x1Arr,y1Arr, 'blue')
	xlabel(r'$x-position$ (AU)')
	ylabel(r'$y-position$ (AU)')
	title('Sun-Earth system computed by Verlet Velocity algorithm for 1500 year with beta= 3, escape velocity=8.9,  stepsize= %s' %dt )
	legend (["Sun","Earth"])
	plt.show()

#PlotSunEarth()

def PlotThreeBodySystem():
#PLOT-THREE BODY SYSTEM

	for i in range(0, stepnumber):
		VerletVelocity(bodies,dt)
		ComputeForcesAndEnergy(bodies)

		xArr[i] = bodies[0].position[0]
		yArr[i] = bodies[0].position[1]

		x1Arr[i] = bodies[1].position[0]
		y1Arr[i] = bodies[1].position[1]


		x2Arr[i] = bodies[2].position[0]
		y2Arr[i] = bodies[2].position[1]

	hold('on')
	plt.plot(xArr,yArr,'ro')
	hold('on')
	plt.plot(x1Arr,y1Arr, 'blue')
	hold('on')
	plt.plot(x2Arr,y2Arr, 'green')
	xlabel(r'$x- position$ (AU)')
	ylabel(r'$y- position$ (AU)')
	title('Three-body system computed by VerletVelocity algorithm with jupiter mass increased with 1000 for 10 year and stepsize=%s' %dt )
	legend (["Sun","Earth","Jupiter"])
	plt.show()

#PlotThreeBodySystem()

def PlotAllPlanets():
#ALL PLANETS
	for i in range(0, stepnumber):
		#EulerCromer(bodies,dt)
		VerletVelocity(bodies,dt)
		ComputeForcesAndEnergy(bodies)

		xArr[i] = bodies[0].position[0]
		yArr[i] = bodies[0].position[1]

		x1Arr[i] = bodies[1].position[0]
		y1Arr[i] = bodies[1].position[1]


		x2Arr[i] = bodies[2].position[0]
		y2Arr[i] = bodies[2].position[1]

		x3Arr[i] = bodies[3].position[0]
		y3Arr[i] = bodies[3].position[1]

		x4Arr[i] = bodies[4].position[0]
		y4Arr[i] = bodies[4].position[1]

		x5Arr[i] = bodies[5].position[0]
		y5Arr[i] = bodies[5].position[1]

		x6Arr[i] = bodies[6].position[0]
		y6Arr[i] = bodies[6].position[1]

		x7Arr[i] = bodies[7].position[0]
		y7Arr[i] = bodies[7].position[1]


	plt.plot(xArr,yArr,'ro')
	hold('on')
	plt.plot(x1Arr,y1Arr, 'blue')
	hold('on')
	plt.plot(x2Arr,y2Arr, 'green')
	hold('on')
	plt.plot(x3Arr,y3Arr, 'purple')
	hold('on')
	plt.plot(x4Arr,y4Arr, 'brown')
	hold('on')
	plt.plot(x5Arr,y5Arr, 'orange')
	hold('on')
	plt.plot(x6Arr,y6Arr, 'pink')
	hold('on')
	plt.plot(x7Arr,y7Arr, 'black')
	plt.axis('equal')
	grid('on')
	xlabel(r'$x-position$ (AU)')
	ylabel(r'$y-position$ (AU)')
	title('All planet of the Solarsystem with NASA values for position and velocities for 10 years and stepsize=%s'%dt )
	legend (["Sun","Earth","Jupiter","Mercury", "Mars", "Venus", "Saturn", "Uranus", "Neptun"])

	plt.show()

#PlotAllPlanets()


def PerihelionPrecessionOfMercury():


	SunPerihelion = Solarsystems(position=[0,0,0], velocity=[0, 0,0], mass=sun_mass)
	MercuryPerihelion =Solarsystems(position=[0.3075*np.sqrt(2.)/2.,0.3075*np.sqrt(2.)/2.,0], velocity=[-12.44*np.sqrt(2.)/2,12.44*np.sqrt(2.)/2,0], mass= mercury_mass)

	bodies= []
	bodies.append(SunPerihelion)
	bodies.append(MercuryPerihelion)

	timeFinal = 100 # 10
	stepnumber = 10000000#varier her 1000000 
	dt= (timeFinal)/float(stepnumber)


	thetaPrevious =0.0 # Perihelion angle of the previous step
	distancePreviousPrevious = 0.0 #Mercury-Sun distance two time step ago
	distancePrevious = 0.0 #the previous time step of Mercury-Sun distance 

	xPosition = np.zeros(stepnumber)
	yPosition = np.zeros(stepnumber)

	norbits=0
	arcsecs=[]
	for i in range(0, stepnumber): #0,dt
		VerletVelocityRel(bodies,dt)
		xPosition[i] = MercuryPerihelion.position[0] - SunPerihelion.position[0]
		yPosition[i] = MercuryPerihelion.position[1] - SunPerihelion.position[1]
		thetaCurrent = np.arctan(yPosition[i]/xPosition[i])
		#arcsecs = thetaCurrent*((360*3600/(2*np.pi))) #radians to arcsecs
		distanceCurrent = sqrt( (MercuryPerihelion.position[0]- SunPerihelion.position[0])**2 + (MercuryPerihelion.position[1] -SunPerihelion.position[1])**2)

		if (distanceCurrent > distancePrevious and distancePrevious < distancePreviousPrevious ): #if we are in Perihelion print angles, which is in radians
			#print "Perihelion angle is:", thetaPrevious #printer vinkelen inne i Perihelion
			#print "Perihelion position:", MercuryPerihelion.position[0], MercuryPerihelion.position[1] #printer vinkelen inne i Perihelion
			arcsecs.append(thetaPrevious*((360*3600/(2*np.pi))) -45*3600)
			norbits+=1

		distancePreviousPrevious = distancePrevious
		distancePrevious = distanceCurrent
		thetaPrevious = thetaCurrent
		
		
	time = linspace(0,timeFinal,norbits)	
	#print norbits, range(len(arcsecs)), arcsecs,time
	#plt.plot(xPosition,yPosition)
	plt.plot(time,arcsecs)
	xlabel('time (year)')
	ylabel(r'$\theta_p$ (Arcseconds)')
	title('Perihelion precession of Mercury over a century with stepsize= %s' %dt )
	#legend (["Sun","Mercury"])
	plt.show()
		


#PerihelionPrecessionOfMercury()	

