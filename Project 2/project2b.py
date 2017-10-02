from pylab import*
from datetime import datetime
import numpy as np
from scipy.sparse import diags
from numpy import linalg as LA
import unittest 
import sys


#Fuction to find the values of cosinus and sinus
def Rotation(A,R,k,l,n):
    s,c=0,0
    
    if (A[k,l] != 0):
        tau = (A[l,l] - A[k,k])/float(2.0*A[k,l])
        if (tau > 0): 
            t = 1.0/(tau + sqrt(1.0 + tau**2))
        else:
            t = -1.0/(-tau + sqrt(1.0 + tau**2))
        
        c = 1/float(sqrt(1+t**(2)))
        s = c*t    
    else: 
        c = 1
        s = 0
    
    a_kk = float(A[k,k])
    a_ll = float(A[l,l])
    

    #Changing the matrix elements with indices k and l
    A[k,k] = c*c*a_kk -2.*c*s*A[k,l] +s*s*a_ll
    A[l,l] = s*s*a_kk + 2.*c*s*A[k,l] + c*c*a_ll
    
    

    A[k,l], A[l,k] = 0.0,0.0

    for i in range(0,n):
        if ((i != k) and (i != l)):
            a_ik =A[i,k] #definerer
            a_il= A[i,l] #definerer 

            A[i,k]= c*a_ik -s*a_il
            A[k,i]= A[i,k] 
            A[i,l]= c*a_il + s*a_ik
            A[l,i]= A[i,l]   
          
        #Finner de nye egenvektorene 

        r_ik = R[i,k]
        r_il = R[i,l]

        R[i,k] = c*r_ik -s*r_il
        R[i,l] = c*r_il +s*r_ik


    
    return A,R


 #Fuction to find maximum matrix element.
def MaximumOffDiagonal(n,A):

	maxx = 0.0
	k,l=0,0

	for i in range(0,n):
		for j in range(i+1,n):
			if (abs(A[i,j]) > maxx):
				maxx= abs(A[i,j])
				l = i
				k = j
			

	return maxx, k,l 



def JacobiMethod(A,R,n):	
	
	for i in range (0,n):
		for j in range (0,n):
			if (i==j):
				R[i,j]=1
			else:
				R[i,j]=0

	
	epsilon = 10**(-8)
	max_number_iterations = n**3
	iterations=0
	maksoff_diagonal,k,l= MaximumOffDiagonal(n, A)

	while ((abs(maksoff_diagonal) > epsilon) and (iterations < max_number_iterations) ):

		A,R = Rotation(A,R,k,l,n)
		maksoff_diagonal,k,l = MaximumOffDiagonal(n,A)
		iterations  += 1

	print "antall iterasjoner=",  iterations
	

	#print "egenverdier:", eigen_valuess #halllooooo
	return A,R,iterations

#definerer potensialene
def potensial(rho):
	return rho*rho

def potensial_twoelectron(rho, w):
	return w*w*rho*rho + 1.0/rho

def values():
	n = 150
	rho_min=0 
	rho_max = 5
	h= (rho_max - rho_min)/float(n+1)
	
	return n, rho_min, rho_max, h

n, rho_min, rho_max, h = values()
rho = linspace(rho_min,rho_max,n)
V= zeros(n)
d= zeros(n)
e= zeros(n)

e = -1/float(h**2)

A= zeros((n,n))
R= zeros((n,n))

for i in range(0,n):
	#d[i] = 2.0/float(h**2) +(i+1)*h 
	d[i]= 2.0/float(h**2) + potensial((i+1)*h)


o = np.array([e*np.ones(n-1),d*np.ones(n),e*np.ones(n-1)])
offset = [-1,0,1]

A= diags(o,offset).toarray()

B = copy(A)
#print "A=",A

s= np.array([0*np.ones(n-1),np.ones(n),0*np.ones(n-1)])
R = diags(s,offset).toarray()
#print  "R= ", R


epsilon= 10**(-8)

tstart= datetime.now()

A_jacobi,R_jacobi,iterations = JacobiMethod(A,R,n)
tend= datetime.now()

print "Tid_jacobi:", tend- tstart

#print "A_jacobi=", np.sort(np.diagonal(A_jacobi)) #egenverdi


#print "R_jacobi=", R_jacobi #riktige egenvektorer 



#Bruker numpy sin egenverdi og egenvektor
tstart1= datetime.now()
eig_vals,eig_vecs =  np.linalg.eig(B) #egenverdi, egenvektor 

eig_vals_sorted = np.sort(eig_vals)

#print "egenverdi_sortet=", eig_vals_sorted
tend1= datetime.now()
#print "Tid numpy:", tend1- tstart1

eigvec=eig_vecs
#print "Egenvektor_numpy=", eig_vecs #RIKTIG 




#UNITTEST 1

def unittest1():
	#sjekker om egenverdiene og egenvektorene til matrise d_unit er lik matlab sin egenverdi og egenvektor
	d_unit =zeros((n,n))

	#HUSK: n+1 matrise
	#d_unit= np.matrix([[2,0,3],[0,3,0], [0,0,3]])
	#d_unit = np.arange(4).reshape(2,2)
	d_unit= np.matrix([[3., 0. ,-1.], [0., 3., 0.], [-1., 0., 5.]])
	d_unit2 = copy(d_unit)
	print "d_unit", d_unit


	D_ddd,R_sss,iterations = JacobiMethod(d_unit,R,n)

	print "unit_egenverdi=",np.sort(np.diagonal(D_ddd))
	#print "unit_egenvektor=", R_sss

	eig_vals,eig_vecs =  np.linalg.eig(d_unit2)
	eig_valssorted = np.sort(eig_vals)
	print "egenverdi_sortet_numpy=", eig_valssorted
	#print "Egenvektor_numpy=", eig_vecs
	return

unittest1() #husk at n= 3, siden vi tester for 3x3 matrise

#UNITTEST NUMBER 2 


def unit_test2():
	#Normene er bevarte, derfor kan man bruke dette. Sjekker om man for identitesmatrisen
	x=np.matrix([[1,2,3],[1,2,3],[1,2,3]])
	eee= LA.norm(R_jacobi*x)
	norm_x = LA.norm(x)
	print "unittest_norm_leftside=", eee
	print "unittest_norm_rightside=", norm_x

unit_test2() #n=3, fordi x er en 3x3 matrise



#Plot non-interaction

def plotnoninteraction():

	sorted_idx = np.argsort(np.diagonal(A_jacobi))[0:4]

	f_eigvec = R_jacobi[:,sorted_idx[0]]
	s_eigvec= R_jacobi[:,sorted_idx[1]]
	t_eigvec = R_jacobi[:,sorted_idx[2]]


	figure()

	plot(rho,f_eigvec**2,'-r')
	hold ('on') 
	plot(rho,s_eigvec**2,'-b')
	hold ('on')
	plot(rho,t_eigvec**2,'-g')  
	xlabel(r'$\rho$')
	ylabel(r'$\psi^2$')
	title('Probability distributions of single non-interaction electron for four lowest-lying states')
	legend ([r"$\psi_0$","$\psi_1$","$\psi_2$","$\psi_3$"])
	show()

plotnoninteraction()


#PLOT NTERACTION

def plotinteraction():
	for w in (0.01, 0.5, 1, 5):
		for i in range(0,n):
			d[i]= 2.0/float(h**2) + potensial_twoelectron( (i+1)*h,w)

		print "w", w
		q = np.array([e*np.ones(n-1),d*np.ones(n),e*np.ones(n-1)])
		offset = [-1,0,1]

		S= diags(q,offset).toarray()

		

		r= np.array([0*np.ones(n-1),np.ones(n),0*np.ones(n-1)])
		T = diags(r,offset).toarray()


		epsilon= 10**(-8)

		

		AA_jacobi,RR_jacobi,iterations = JacobiMethod(S,T,n)
		
		sorted_idx = np.argsort(np.diagonal(AA_jacobi))[0:3]
		#print "sorted_idx=", sorted_idx
		print "sortet diagonal egenverdi= ", np.sort(np.diagonal(AA_jacobi))[0:3]

	
		first_eigvec = RR_jacobi[:,sorted_idx[0]]
		second_eigvec= RR_jacobi[:,sorted_idx[1]]
		third_eigvec = RR_jacobi[:,sorted_idx[2]]

	#plotting wavefunction

		
		figure()

		plot(rho,first_eigvec**2,'-r')
		hold ('on') 
		plot(rho,second_eigvec**2,'-b')
		hold ('on')
		plot(rho,third_eigvec**2,'-g')  
		xlabel(r'$\rho$')
		ylabel(r'$\psi^2$')
		title('Probability distributions of two interacting electrons for four lowest-lying states for $\omega$= %s' %w)
		legend ([r"$\psi_0$","$\psi_1$","$\psi_2$","$\psi_3$"])
		show()

		
	return 

plotinteraction()

