import numpy
from matplotlib.pyplot import *

x = [50, 75, 100,125, 150, 175, 200]

y = [4039, 9229, 16474,25826, 37443, 50991, 66827]



coefficients = numpy.polyfit(x, y, 3)
polynomial = numpy.poly1d(coefficients)
ys = polynomial(x)
#print coefficients
print polynomial

plot(x, y, 'o' )
plot(x, ys, label= '$0.0001013n^{3} + 1.662n^{2} - 2.445n - 1.214$')
ylabel('Number of Similarity transformation')
xlabel('Grid points (n)')
legend()
title('Extracting function of the estimated number of similarity transformations.')
show()
