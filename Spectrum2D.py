from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

t=2000
x11=[]
x22=[]
e=3
K=e/2*np.pi
a=0.5
b=0.1
x1=a
x2=b

dq=np.random.uniform(0,0.05)
dp=np.random.uniform(0,0.05)
w0=np.array([dq,dp])
dq1=np.random.uniform(0,0.05)
dp1=np.random.uniform(0,0.05)
w1=np.array([dq1,dp1])
dict={}

alpha0=[]
alpha1=[]


for i in range(t):
    #Initial state
    x1n,x2n=chirikov(x1,x2,K)
    x1=x1n 
    x2=x2n 

    x11.append(x1)
    x22.append(x2)  

    Mat2=np.array([[1, K*np.cos(x2)],[1, 1+K*np.cos(x2)]]) 

    dp,dq=np.dot(Mat2,w0)

    dp1, dq1= np.dot(Mat2, w1)

    w0=np.array([dp,dq])
    w1=np.array([dp1,dq1])

    #Gram-Schmidt
    
    u0=w0
    a0= np.sqrt(dp**2+dq**2)
    w0= u0/a0
    alpha0.append(a0)

    u1=w1-np.inner(w1,w0)*w0
    a1= np.sqrt(dq1**2+dp1**2)
    w1= u1/a1
    alpha1.append(a1)
    
#print(alpha1, alpha0)

lalpha0=np.log(alpha0)
lalpha1=np.log(alpha1)

lyap0=[]
for i in range (1,len(lalpha0)):
    l0=sum (lalpha0[0:i])/((i))
    lyap0.append(l0)


lyap1=[]
for i in range (1,len(lalpha1)):
    l1=sum (lalpha1[0:i])/((i))
    lyap1.append(l1)

lyap0=lyap0[1:]
lyap1=lyap1[1:]


plt.plot((lyap0), label="Lyap1")
plt.plot((lyap1), label='Lyap2')

#plt.ylim(-1,9)
plt.legend()
plt.ylabel('Lyapunov exponent')
plt.title("Lyapunov Spectrum of the 2D Standard Map with K= "+str(e)+", and initial conditions x1= "+str(a)+", x2= "+str(b)+"")
plt.xlabel('steps')
plt.show()

ll=np.arange(len(lyap0))
S=[]

for i in range (len(lyap0)):
    s=lyap0[i]+lyap1[i]
    S.append(s)
#print(S)
sl= np.arange(len(S))

plt.plot(sl, S,  "r", label='Sum')
plt.title("Symmetry check of Lyapunov Spectrum of the 2D Standard Map with K= "+str(e)+", and initial conditions x1= "+str(a)+", x2= "+str(b)+"")
plt.ylabel("Difference between LEs")
plt.xlabel("Iteration")
plt.show()
