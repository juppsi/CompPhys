
import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import sys, math


def tridiag(v,nx,alpha):
	
	#tridiagonal matrix, reused from project 1

	b= 2.0 + 2.0*alpha #main diagonal fo CRANK-NICOLSEN
	# b= 1.0 + 2.0*alpha #main diagoanl for BACKWARD Euler
	a= -alpha #-alpha #low- diagonal
	c= -alpha #-alpha #upper diagonal
	
	b_hat = np.zeros(nx+2) + b


	#forward substitue
	for i in range (1,nx+1):
		b_hat[i]=b-a*c/float(b_hat[i-1])
		v[i]=v[i]- a*v[i-1]/float(b_hat[i-1])
	
	#backward substitue
	for i in reversed(range(1, nx+1)):
		v[i]=(v[i]-c*v[i+1])/float(b_hat[i])

	return v
	

def ForwardEuler(alpha,u,nt,nx):
	#Explicit forward euler SCHEME
	for it in range(1,nt+2):
		for ix in range(1,nx+1):
			u[it,ix]= alpha*u[it-1,ix-1] + (1.0-2.0*alpha)*u[it-1,ix] + alpha*u[it-1,ix+1]	
	return u

def BackwardEuler(alpha, u,nt, nx):
	#Implicit Backward euler scheme by using gauss elimination of tridiagonal matrix (reusing from project 1)
	for it in range(1,nt+2):
		u[it,:] = tridiag(u[it-1,:],nx,alpha)
	return u

def CrankNicolson(alpha,u,nt,nx):
	#Crank-Nicolson scheme
	for it in range(1,nt+2):
		for ix in range(1,nx+1):
			u[it,ix]= alpha*u[it-1,ix-1] + (2.0-2.0*alpha)*u[it-1,ix] + alpha*u[it-1,ix+1]
		u[it]=tridiag(u[it],nx,alpha)

	return u


def AnalyticalOneDim():
	#Analytical solution for the 1D diffustion equation 

	for t in [0.003, 0.03, 0.1,0.3, 0.6]: #time evolution
		x=np.linspace(0,1,400)
		u=np.zeros(len(x))
		n=1000

		for n in xrange(1,n+1):

			u += np.sin(np.pi*n*x)*np.exp(-np.pi*n**2*t)/float(np.pi*n)

		u= x + 2*u

		ts = "t =%s" % str(t)
		plt.plot(u,x, label= ts)
		plt.xlabel('$x$')
		plt.ylabel('$u(x,t)$')
		plt.axis([0, 1.01,0.0, 1.0])
		plt.savefig('analytical1Ddd.pdf')
		plt.title('Analytical soluion of one- dimenstional diffusion equation')

	plt.tight_layout()
	plt.legend(loc= 'best')
	plt.show()

#AnalyticalOneDim()



def Analytical2D():
	#Analytical solution to the 2D diffusion equation

	n= 100

	x= np.linspace(0,1,300)
	t= np.linspace(0,1,300)

	X,Y =np.meshgrid(x,t)

	u=np.zeros(X.shape)

	for n in xrange(1,2*n,2):
		#u += -(2./(pi*n))*sin(n*pi*x)*exp(-n**2*pi**2*t)
		u +=  -np.sin(np.pi*n*X)*np.sinh(np.pi*n*Y)/float(np.pi*n*np.sinh(np.pi*n))

	u= 4.0*u/float(np.pi)

	fig= plt.figure()
	plt.hold('on')
	ax = fig.add_subplot(111)
	plt.xlabel('$x$')
	plt.ylabel('$u(x,t)$')
	
	plt.title('Analytical soluion of two- dimenstional diffusion equation')
	plt.savefig('analytical2Dd.pdf')
	cont=ax.contour(X,Y,u)
	plt.show()


#Analytical2D()


def numerical1D():

	L = 1.0
	dx = 1./float(10) #h
	dt = (dx**2)/float(2.0)
	alpha = dt/float(dx**2)


	for Tfinal in [ 0.003, 0.03, 0.09, 0.3]: #time evolution

		nt = np.int(Tfinal/float(dt) -2.) #grid points
		nx = np.int(L/float(dx) - 2.) #grid points

		time = np.linspace(0,Tfinal,nt+2) #grid time

		x = np.linspace(0,L,nx+2)
		u = np.zeros((nt+2,nx+2))

		# boundary condition 
		u[:,-1] = 1
		
		v = u.copy()
		r = u.copy()
		s = u.copy()

		vv = ForwardEuler(alpha,v,nt,nx)
		rr = BackwardEuler(alpha,r,nt,nx)
		ss = CrankNicolson(alpha,u,nt,nx)
		
		FE_ = "FE" + ", t=%s" % str(Tfinal) #FE= forward Euler
		BE_ = "BE" + ", t =%s" % str(Tfinal) #BE = backward Euler
		CN_ = "CN" + ", t= %s" % str(Tfinal) #CN= Crank-Nicolson

		#CN_ = "t =%s" % str(Tfinal)
		
		plt.plot(x,v[nt+1],'-b', label= FE_)
		plt.hold('on')
		plt.plot(x,r[nt+1],'*g',label= BE_)
		plt.hold('on')
		plt.plot(x,ss[nt+1],'ro',label= CN_)
		plt.hold('on')
		plt.xlabel('$x$')
		plt.ylabel('$u(x,t)$')

	plt.title('1D diffusion equation computed by finite difference scheme with \n dx= %s, dt= %s, alpha= %s '%(dx,dt,alpha))
	plt.legend(loc= 'best')
	plt.tight_layout()
	#plt.savefig('alle1over10.pdf')
	plt.show()

numerical1D()

def JacobiImplicit(u,alpha,nx,ny,tsteps):
	#Using Implicit scheme to solve the 2D difffusion equation with Jacobi method as a iterative method.

	epsilon = 10**(-8) #tolerance
	difference = epsilon + 1 
	iterations = 0
	maxiteration = nx**3

	uPrev = np.zeros((nx+2,ny+2))
	u_temp = np.zeros((nx+2,ny+2))

	uPrev = u
	
	u_temp = uPrev

	while(difference > epsilon and iterations < maxiteration):
		difference = 0.0
		for k in xrange(tsteps):
			for i in range(1, nx+1):
				for j in range(1,ny+1):
					u[i,j]=(1.0/float(1.+4.*alpha))*(uPrev[i,j] + alpha*(u_temp[i+1,j] + u_temp[i-1,j] + u_temp[i,j+1] + u_temp[i,j-1]))

									
			iterations +=1


			difference /= nx*ny

			u_temp = uPrev
			uPrev = u

			"""
			#Plotting the Time evolution
			if iterations%100 == 0:
				L=1.0
				x= np.linspace(0,L,nx+2)
				y= np.linspace(0,L,ny+2)
				X,Y =np.meshgrid(x,y)
			
				#u=np.zeros(X.shape)

				fig= plt.figure()
				plt.hold('on')
				ax = fig.add_subplot(111)

				cont=ax.contour(X,Y,u)
				plt.tight_layout()
				#plt.savefig('timeevolution.pdf')
				plt.xlabel('$x$')
				plt.ylabel('$y$')

				plt.title('Time evolution of 2D diffusion equation')
				#plt.show()
				"""
	return u


def PlotImplicit2D():

	n=50 #already at n=30 we reach the equilibrium state
	dx= 1./float(n-1) 

	dt= 0.001 #(dx**2)/float(2.0) #from stability

	tsteps = 1000   # equilibrium: int(0.5/dt)

	L = 1.0
	Tfinal= 0.1

	alpha= dt/float(dx**2) #stability

	nt= np.int(Tfinal/float(dt)-2) #grid time

	nx = np.int(L/float(dx)-2) #grid points
	ny = np.int(L/float(dx)-2) #grid points 

	x= np.linspace(0,L,nx+2)
	y= np.linspace(0,L,ny+2)

	time= np.linspace(0,Tfinal,n)

	u = np.zeros((nx+2,ny+2))


	
	#Velg enten begge, eller en av dem.
	#boundary condition 
	for i in range(0, ny+2):
		u[i,nx+1] = 1.0
		u[i,0]=1

	#boundary condition
	for j in range(0,nx+2):
		u[ny+1,j] = 1.0
		#u[0,j]=1

	u = JacobiImplicit(u,alpha,nx,ny,tsteps)
	
	X,Y =np.meshgrid(x,y)
	


	fig= plt.figure()
	plt.hold('on')
	ax = fig.add_subplot(111)

	
	cont=ax.contour(X,Y,u)
	plt.tight_layout()
	#plt.axis([0.99, 1.0,0.0, 0.03])
	plt.xlabel('$x$')
	plt.ylabel('$y$')

	plt.title('2D diffusion equation with dx= %1.2f, dt= %1.3f, alpha= %1.1f, t= %1.1f '%(dx,dt,alpha,Tfinal))
	plt.legend(loc= 'best')
	#plt.savefig('boundaryhallo1.pdf')

	plt.show()
	
#PlotImplicit2D()



