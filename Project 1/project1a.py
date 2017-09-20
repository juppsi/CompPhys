from pylab import*
from matplotlib.pyplot import*
import sys
from datetime import datetime

#for n in 10, 100, 1000: 
n=10000000
x_n= 1
x_0= 0

x= linspace(x_0,x_n,n+2)
h= (x_n - x_0)/float(n+1)

f_tilde = zeros(size(x))
b_hat   = zeros(size(x))
v = zeros(size(x))
u =zeros(size(x))
relative_error = zeros(size(x))


a=-1*ones(n+2)
c=-1*ones(n+2)
b=2*ones(n+2)




f= lambda x: 100*(exp(-10*x)) #source term 
u= lambda x: 1-(1-exp(-10))*x- exp(-10*x) #exact solution

f_tilde=f(x)*(h**2)


def ForwardEuler():
	b_hat = b
	for i in range (2,n+1):
		b_hat[i]=b[i]-c[i-1]*a[i-1]/float(b_hat[i-1])
		f_tilde[i]=f_tilde[i]- f_tilde[i-1]*a[i-1]/float(b_hat[i-1])
	return b_hat, f_tilde
	

def BackwardEuler(b_hat, f_tilde):
	
	v =zeros(size(x))
	v[n]=f_tilde[n]/float(b_hat[n])
	
	for i in reversed(range(1, n)):
		v[i]=(f_tilde[i]-v[i+1]*c[i])/float(b_hat[i])
		#v[i+1]=(f_tilde[i]-v[i+2]*c[i])/float(b_hat[i])
	return v

tstart= datetime.now()
print tstart


b_hat, f_tilde = ForwardEuler()
v = BackwardEuler(b_hat,f_tilde)

tend= datetime.now()
print tend

print "Endring i tid:", tend- tstart 

"""
v_correct = zeros(size(v)+1)
v_correct[0] = 0
v_correct[1:] = v[:]
#plot(x, v_correct[:-1], "r", label = "Numerical solution")
"""

figure()

plot(x, v, "r", label = "Numerical solution")
hold('on')
plot(x, u(x), "k--", label = "Closed form term") 
title("n = %.i" % n)
legend()
show()

u_vec = u(x)


for i in range(1,n+1):
	#relative_error[i]= log(abs((u(x[i])-v[i])/abs(u(x[i]))))
	relative_error[i]= abs(u_vec[i]-v[i])/abs(u_vec[i])

relative_error = log10(relative_error[1:n+1])
relative_error = max(relative_error)
print "relative_error", relative_error

