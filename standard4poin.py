from Maps4 import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

x11=[]
x22=[]
x33=[]
x44=[]


K=0.5
B= 0.05
a=0.55
b=0.10
c=0.054 #0.54 #0.005
d=0.01
x1=a
x2=b
x3=c
x4=d

for i in range(1000):
    x1n,x2n, x3n, x4n=standard4d(K, B, x1,x2,x3,x4)
    print(x2n)
    x1=x1n %1
    x2=x2n %1
    print(x2)
    x3=x3n %1
    x4=x4n %1
    #print(x1)

    x11.append(x1)
    x22.append(x2) 
    x33.append(x3)
    x44.append(x4)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

sol =ax.scatter(x11,x22, x33, 'o', c= x44, cmap=cm.jet)
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('x3')
fig.colorbar(sol, shrink=0.5, aspect=5, label='x4')
plt.title("Phase space of 4D Standard Map with K =" + str(K) + " and initial conditions x1= "+ str(a)+", x2= "+str(b)+", x3= "+str(c)+", x4= "+str(d)+"")
plt.show()
"""
plt.scatter(x11, x22)
plt.show()"""