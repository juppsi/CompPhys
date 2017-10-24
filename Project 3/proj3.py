from pylab import*
import numpy as np
from math import pi



def values():
    timeFinal = 1.0 #1 year
    stepnumber = 1000
    h = (timeFinal)/float(stepnumber) #stepsize 
    time = 0.0
    x = 1 #1AU
    y =  0.0
    vx = 0.0
    vy = 2.*pi
    r = sqrt(x**2 + y**2)

    
    return timeFinal,stepnumber,h, time, x, y,vx,vy,r

timeFinal,stepnumber,h, time, x, y,vx,vy,r = values()



def EulerForward(timeFinal,h, x, y,vx,vy,r):
	time = 0.0
	while (time <= timeFinal):

		ax= -4*pi**(2)*x/float(r**(3))
		ay = -4*pi**(2)*y/float(r**(3))

		x = x +vx*h
		y =  y + vy*h 

		vx =  vx + ax*h
		vy =  vy + ay*h

		r = sqrt(x**2 + y**2)

		time =  time + h
		
		print time, x, y, vx, vy


EulerForward(timeFinal, h, x, y,vx,vy,r)


def Verletmethod(timeFinal, h, x, y,vx,vy,r):
	time = 0.0
	while (time <= timeFinal):
		ax= -4*pi**(2)*x/float(r**(3))
		ay = -4*pi**(2)*y/float(r**(3))

		x= x + h*vx + h**(2)*ax/float(2)
		y= y + h*vy + h**(2)*ay/float(2)

		vx= vx + ax**(2)*h/float(2)
		vy= vy + ay**(2)*h/float(2)

		r = sqrt(x**2 + y**2)

		time = time + h

	plot (time, x, time, y)
	return 

Verletmethod(timeFinal, h, x, y,vx,vy,r)






