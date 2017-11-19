import numpy as np
import sys, math
import matplotlib.pyplot as plt

import matplotlib as mpl



listfiles20 = ["temp_T1_L20.txt", "temp_T2_L20.txt", "temp_T3_L20.txt", "temp1_T1_L20.txt", "temp1_T2_L20.txt"]  #Latticesize =20

listfiles40= ["T1_L40.txt", "T2_L40.txt", "T3_L40.txt", "temp1_T1_L40.txt", "temp1_T2_L40.txt"] #Latticesize =40


listfiles60 = ["temp_T1_L60.txt", "temp_T2_L60.txt", "temp_T3_L60.txt", "temp1_T1_L60.txt", "temp1_T2_L60.txt"]  #Latticesize =60

listfiles80 = ["temp_T1_L80.txt", "temp_T2_L80.txt", "temp_T3_L80.txt", "temp1_T1_L80.txt", "temp1_T2_L80.txt"]  #Latticesize =80

listfiles100 = ["temp_T1_L100.txt", "temp_T2_L100.txt", "temp_T3_L100.txt", "temp1_T4_L100.txt", "temp1_T6_L100.txt"]  #Latticesize =100


temperature = []
EAverage = []
MagAbsAverage= []
SpecificHeat =[]
Susceptibility = []


for file in listfiles40: #change listfiles to chosen latticesize
	lines = open(file, 'r').readlines()


	for line in lines:
		line = line.split()

		temperature.append(float(line[0]))
		EAverage.append(float(line[1]))
		MagAbsAverage.append(float(line[2]))
		SpecificHeat.append(float(line[3]))
		Susceptibility.append(float(line[4]))




#PLOTTING EXPECTATION VALUES MEAN ENERGY, MEAN ABOSLUTE MAGNETIZATION, SPECIFIC HEAT, SUSEPTIBILITY

fig, ax = plt.subplots(nrows=2, ncols=2)
plt.suptitle('40 x 40 Ising model with $1000000$ MCS and $\Delta t = 0.02$', fontsize=15)


plt.subplot(2, 2, 1)
plt.plot(temperature,EAverage) 
plt.xlabel('Temperature')
plt.ylabel("$<E/N>$")
plt.grid('on')
#plt.savefig('100x100En.pdf')
#plt.title('100 x 100 Ising model with $100000$ MCS and $\Delta t = 0.05$')
plt.tight_layout()


plt.subplot(2, 2, 2)
plt.plot(temperature,MagAbsAverage)
plt.xlabel('Temperature')
plt.ylabel("$<|M|/N>$")
#plt.title('100 x 100 Ising model with $100000$ MCS and $\Delta t = 0.05$')
plt.grid('on')
#plt.savefig('80x80Mag.pdf')
plt.tight_layout()


plt.subplot(2, 2, 3)
plt.plot(temperature,SpecificHeat)
plt.xlabel('Temperature')
plt.ylabel("Specific heat $C_V$")
#plt.title('100 x 100 Ising model with $100000$ MCS and $\Delta t = 0.05$')
plt.tight_layout()
plt.grid('on')
#plt.savefig('100x100CV.pdf')


plt.subplot(2, 2, 4)
plt.plot(temperature,Susceptibility)
plt.xlabel('Temperature')
plt.ylabel("Susceptibility $\chi$")

plt.tight_layout()
plt.grid('on')

plt.subplots_adjust(top=0.9)
plt.savefig('40x40.pdf')
plt.show()







