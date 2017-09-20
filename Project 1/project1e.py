

from scipy.sparse import diags
import numpy as np
import pprint
import scipy
import scipy.linalg
from datetime import datetime

"""
A= scipy.array([2, -1, 0, 0,0,..., 0 ], [-1,2, -1,0, 0,..., 0], [0, -1, 2, -1,0,..., 0], [...],[0,...,0,0,-1,2,-1], [0,0,0,...,0,-1,2])
P,L,U=scipy.linalg.lu(A)

"""


n = 10000 #10,100,1000
a= 1*np.ones(n-1)
b=-2*np.ones(n)
c=1*np.ones(n-1)




k = np.array([a,b,c])
offset = [-1,0,1]

tstart= datetime.now()
print tstart
A = diags(k,offset).toarray()
P,L,U=scipy.linalg.lu(A)

tend= datetime.now()
print tend

print "Endring i tid:", tend- tstart 


"""
print "A:"
pprint.pprint(A)

print "P:"
pprint.pprint(P)

print "L:"
pprint.pprint(L)

print "U:"
pprint.pprint(U)

print "LU:"
pprint.pprint(L*U)

"""
