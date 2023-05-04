from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

alpha=[]
t=50000
x11=[]
x22=[]
e=0.8
K=e/2*np.pi
a=0.11
b=0.8
x1=a
x2=b

dp=0.01
dq=0.01

for i in range(t):
    x1n,x2n=chirikov(x1,x2,K)
    x1=x1n 
    x2=x2n 
    #print(x1)

    x11.append(x1)
    x22.append(x2)  

    Mat2=np.array([[1, K*np.cos(x2)],[1, 1+K*np.cos(x2)]]) 

    dp,dq=np.dot(Mat2, np.array([dp,dq]))

    #print(dp, dq)

    a1=np.sqrt(dp**2+dq**2)
    alpha.append(a1)
    dp=dp/a1
    dq=dq/a1 #Renormalization """

    #print(a)
    #plt.plot(x1,x2,'o', markersize=1, color='red')

#plt.show()

#print(alpha)
alpha=alpha[1:]
lalpha=np.log(alpha)
#print(lalpha)

lyap=[]
for i in range (1,len(lalpha)):
    l=sum (lalpha[0:i])/(i+1)
    lyap.append(l)

print(lyap[-1])

ll=np.arange(len(lyap))
plt.plot(ll,lyap, 'r')
plt.xlabel('t')
plt.ylabel('Lyapunov exponent')
plt.title("Lyapunov exponent for the 2D standard map with K= "+str(e)+", x1= "+str(a)+" and x2= "+str(b)+"")
plt.show()