from Maps4 import *
import numpy as np
import matplotlib.pyplot as plt

x11=[]
x22=[]
e=7
K=e/2*np.pi
a=0.5
b=0.1
x1=a
x2=b
for i in range(10000):
    x1n,x2n=chirikov(x1,x2,K)
    x1=x1n 
    x2=x2n 
    #print(x1)

    x11.append(x1)
    x22.append(x2)  

plt.plot(x11,x22,'o', markersize=1, color='red')
plt.xlabel('x1')
plt.ylabel('x2')
plt.title("Poincare section of 2D Standard Map with K =" + str(e) + " and initial conditions x1= "+ str(a)+", x2= "+str(b)+"")
plt.show()