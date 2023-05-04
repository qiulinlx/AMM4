from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

alpha=[]
t=3000
K=0.5
B= 0.05
a1=0.55
b1=0.10
c1=0.54 #0.54 #0.005
d1=0.01
x1=a1
x2=b1
x3=c1
x4=d1

dx1=0.01
dx2=0.01
dx3=0.01
dx4=0.01


for i in range(t):
    x1n,x2n, x3n, x4n=standard4d(K, B, x1,x2,x3,x4)
    #print(x2n)
    x1=x1n %1
    x2=x2n %1
    #print(x2)
    x3=x3n %1
    x4=x4n %1

    Mat4=np.array([[1+K*np.cos(2*np.pi*x1)+B*np.cos(2*np.pi*(x3-x1)), 1, -B*np.cos(2*np.pi*(x3-x1)), 0], 
                  [K*np.cos(2*np.pi*x1)+ B*np.cos(2*np.pi*(x1-x3)) , 1, -B*np.cos(2*np.pi*(x3-x1)), 0],
                  [-B*np.cos(2*np.pi*(x1-x3)), 0, 1+K*np.cos(2*np.pi*x3)+B*np.cos(2*np.pi*(x1-x3)), 1],
                  [-B*np.cos(2*np.pi*(x1-x3)), 0, K*np.cos(2*np.pi*x3)+B*np.cos(2*np.pi*(x1-x3)), 1]])

    dx1, dx2, dx3, dx4 =np.dot(Mat4, np.array([dx1,dx2, dx3, dx4]))

    #print(dp, dq)

    a=np.sqrt(dx1**2+dx2**2+dx3**2+dx4**2)
    alpha.append(a)
    dx1=dx1/a
    dx2=dx2/a #Renormalization 
    dx3=dx3/a
    dx4=dx4/a

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
plt.plot(np.log10(ll), np.log10(lyap), 'g')
plt.xlabel('Log(t)')
plt.ylabel('Log(Lyapunov exponent)')
plt.title("mLE of the 4D standard map with K= "+str(K)+", x1= "+str(a1)+" and x2= "+str(b1)+", x3= "+str(c1)+" and x4= "+str(d1)+"")
plt.show()