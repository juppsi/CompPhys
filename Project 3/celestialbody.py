from pylab import*
import numpy as np
#from math import pi, arctan

class Solarsystems:

	parameters = { 'position':[float,float,float], 'velocity':[float, float, float], 'force': [float, float, float],'oldAccel':[float, float, float]}
	paras= {'mass':float, 'radius':float, 'AngularMomentum':float, 'spin': float}
	
	# [x, y, z]
	position =None
	
	# [vx, vy, vz]
	velocity=None

	# [fx, fy, fz]
	force= None

	oldAccel=None

	
	def __init__(self, position,velocity,mass):
		self.mass = mass
		self.position = position
		self.velocity= velocity 
		self.force = [0,0,0]
		self.oldAccel = [0,0,0]
		self.potensialEn = 0.
		self.kineticEn = 0.
		self.AngularMomentum= 0.
		#print force
		#print oldAccel
	


def ComputeForcesAndEnergy(bodies):

	for planet in (bodies):
		planet.force[0]=0.0
		planet.force[1]=0.0
		planet.force[2]=0.0
		planet.potensialEn = 0.0
		planet.kineticEn = 0.0

		#iterere over alle bodiene og finner kraften
	for i in range (0,len(bodies)):
		body1 = bodies[i]
		for j in range(i+1,len(bodies)):
			body2 = bodies[j]

			dx= body2.position[0] - body1.position[0]
			dy = body2.position[1] - body1.position[1]
			dz= body2.position[2] - body1.position[2]

			dr_squared = dx**(2) + dy**(2) + dz**(2)
			r = sqrt(dr_squared) #avstanden
			G = 4.0*np.pi**(2)


			Force_calculate = -G*body1.mass * body2.mass/float(r**(3))
			PotEnergy = G*body1.mass * body2.mass/ float(r)
			
			fx_ny= Force_calculate*dx
			fy_ny= Force_calculate*dy
			fz_ny = Force_calculate*dz

			body1.force[0] -= fx_ny
			body1.force[1] -= fy_ny
			body1.force[2] -= fz_ny
			body1.potensialEn += PotEnergy

			body2.force[0] += fx_ny
			body2.force[1] += fy_ny
			body2.force[2] += fz_ny
			body2.potensialEn += PotEnergy



	for planet in (bodies):
		planet.kineticEn = 0.5*planet.mass*(planet.velocity[0]**(2) + planet.velocity[1]**(2) + planet.velocity[2]**(2))
		
		V= sqrt(planet.velocity[0]**(2) + planet.velocity[1]**(2) + planet.velocity[2]**(2))
		R = sqrt(planet.position[0]**(2)+ planet.position[1]**(2)+ planet.position[2]**(2))

		planet.AngularMomentum = V*R*planet.mass
	
	PotensialEn_sum,  KineticEn_sum , TotalEnergy, AngularMomentum_sum= 0.0, 0.0, 0.0, 0.0

	for planet in (bodies):
		KineticEn_sum += planet.kineticEn
		PotensialEn_sum += planet.potensialEn
		TotalEnergy += KineticEn_sum + PotensialEn_sum
		AngularMomentum_sum += planet.AngularMomentum

	#print "PotensialEn_sum", PotensialEn_sum, "KineticEn_sum", KineticEn_sum, "TotalEnergy", TotalEnergy, "AngularMomentum_sum", AngularMomentum_sum


def ComputeForcesAndEnergyRel(bodies):
	for planet in (bodies):
		planet.force[0]=0.0
		planet.force[1]=0.0
		planet.force[2]=0.0

		#iterere over alle bodiene og finner kraften
	for i in range (0,len(bodies)):
		body1 = bodies[i]
		for j in range(i+1,len(bodies)):
			body2 = bodies[j]

			dx= body2.position[0] - body1.position[0]
			dy = body2.position[1] - body1.position[1]
			dz= body2.position[2] - body1.position[2]

			dr_squared = dx**(2) + dy**(2) + dz**(2)
			r = sqrt(dr_squared) #avstanden
			G = 4.0*np.pi**(2)


			l1= (body2.position[1]*body2.velocity[2]-body2.position[2]*body2.velocity[1])**2
			l2= (body2.position[0]*body2.velocity[2]- body2.position[2]*body2.velocity[0])**2
			l3= (body2.position[0]*body2.velocity[1]- body2.position[2]*body2.velocity[0])**2
			
			l = l1 -l2 +l3
			c=63145
			Force_calculate = -G*body1.mass * body2.mass/float(r**(3))*(1.0 +(3.0*l**(2))/float(r**2*c**2)) #PERIHELION

			fx_ny= Force_calculate*dx
			fy_ny= Force_calculate*dy
			fz_ny = Force_calculate*dz

			body1.force[0] -= fx_ny
			body1.force[1] -= fy_ny
			body1.force[2] -= fz_ny

			body2.force[0] += fx_ny
			body2.force[1] += fy_ny
			body2.force[2] += fz_ny



def ForwardEuler(bodies, dt):
	ComputeForcesAndEnergy(bodies)
	for planet in (bodies):
		for i in range(0,3):
			planet.position[i] += planet.velocity[i]*dt
			planet.velocity[i] += (planet.force[i]/float(planet.mass))*dt
	


def EulerCromer(bodies, dt):
	ComputeForcesAndEnergy(bodies)
	#j=0
	for planet in (bodies):
		#print j, planet.position, planet.velocity, planet.force
		#j+=1
		for i in range(0,3):
			planet.velocity[i] += (planet.force[i]/float(planet.mass))*dt
			planet.position[i] += planet.velocity[i]*dt


def VerletVelocity(bodies,dt):
	ComputeForcesAndEnergy(bodies)

	for planet in (bodies):
		for i in range(0,3):
			planet.oldAccel[i] = planet.force[i]/float(planet.mass)
			planet.position[i] += planet.velocity[i]*dt + 0.5*(dt**(2))*planet.force[i]/float(planet.mass)

	#new acceleration and velocity

	ComputeForcesAndEnergy(bodies)
	
	for planet in (bodies):
		for i in range(0,3):
			planet.velocity[i] += (((planet.force[i]/float(planet.mass)) + planet.oldAccel[i]) * 0.5*dt)


def VerletVelocityRel(bodies,dt):
	ComputeForcesAndEnergyRel(bodies)

	for planet in (bodies[1:]):
		for i in range(0,3):
			planet.oldAccel[i] = planet.force[i]/float(planet.mass)
			planet.position[i] += planet.velocity[i]*dt + 0.5*(dt**(2))*planet.force[i]/float(planet.mass)


	#new acceleration and velocity

	ComputeForcesAndEnergyRel(bodies)
	for planet in (bodies[1:]):
		for i in range(0,3):
			planet.velocity[i] += (((planet.force[i]/float(planet.mass)) + planet.oldAccel[i]) * 0.5*dt)



