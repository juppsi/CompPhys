

from pylab import*
import sys
from datetime import datetime


n=1000000
x_n= 1
x_0= 0

x= linspace(x_0,x_n,n+2)
h= (x_n - x_0)/float(n+1)

f_tilde = zeros(size(x))
b_hat = zeros(size(x))
v = zeros(size(x))
relative_error = zeros(size(x))

b=2*ones(n+2)



f= lambda x: 100*(exp(-10*x)) #source term 
u= lambda x: 1-(1-exp(-10))*x- exp(-10*x) #exact solution


f_tilde=f(x)*(h**2)




"""
def ForwardBackwardEulerEasy():
	b_hat[0]=2
	for i in range (1,n+1):
		alpha = 1/b_hat[i-1]
		b_hat[i]=b[i]-alpha
		f_tilde[i]=f_tilde[i]- f_tilde[i-1]*alpha
	return b_hat, f_tilde
	
	v =zeros(size(x))
	
	v[n]=f_tilde[n]/float(b_hat[n])	

	for i in reversed(range(1, n-1)):
		v[i]=(f_tilde[i]-v[i+1])/float(b_hat[i])
	return v
"""
def Forward_Euler():
	b_hat = b
	for i in range (2,n+1):
		b_hat[i]=b[i]-1/float(b_hat[i-1])
		f_tilde[i]=f_tilde[i]- f_tilde[i-1]/float(b_hat[i-1])
	return b_hat, f_tilde
	
def Backward_Euler():
	v =zeros(size(x))
	v[n]=f_tilde[n]/float(b_hat[n])
	
	for i in reversed(range(1, n)):
		v[i]=(f_tilde[i]-v[i+1])/float(b_hat[i])
		#v[i+1]=(f_tilde[i]-v[i+2]*c[i])/float(b_hat[i])
	return v

tstart= datetime.now()
print tstart

b_hat, f_tilde = Forward_Euler()
v = Backward_Euler()

tend= datetime.now()
print tend

print "Endring i tid:", tend- tstart 
