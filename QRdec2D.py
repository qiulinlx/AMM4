from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

t=10
x11=[]
x22=[]
e=4
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

    dp,dq=np.dot(Mat2, w0)

    dp1, dq1= np.dot(Mat2, w1)

    w0=np.array([dq,dp])
    w1=np.array([dq1,dp1])
    #print(dq)
    W=np.array([[dp, dp1],[dq, dq1]])
    #print(W)
    #QR-decomposition
    q,r =np.linalg.qr(W)
    a0=r[0,0]
    a1=r[1,1]

    alpha0.append(a0)
    alpha1.append(a1)

    w0=

    
print(alpha1, alpha0)

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


plt.plot((lyap0), label="Lyap0")
plt.plot((lyap1), label='Lyap1')

#plt.ylim(-1,9)
plt.legend()
plt.show()
