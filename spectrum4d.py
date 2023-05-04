from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

alpha0=[]
alpha1=[]
alpha2=[]
alpha3=[]

t=3000
K=0.5
B= 0.05
a=0.55
b=0.10
c=0.005 #0.54 #0.005
d=0.01
x1=a
x2=b
x3=c
x4=d

dx1=np.random.uniform(0,0.02)
dx2=np.random.uniform(0,0.02)
dx3=np.random.uniform(0,0.02)
dx4=np.random.uniform(0,0.02)
w0=np.array([dx1,dx2, dx3, dx4])

dx11=np.random.uniform(0,0.02)
dx12=np.random.uniform(0,0.02)
dx13=np.random.uniform(0,0.02)
dx14=np.random.uniform(0,0.02)
w1=np.array([dx1,dx2, dx3, dx4])

dx21=np.random.uniform(0,0.02)
dx22=np.random.uniform(0,0.02)
dx23=np.random.uniform(0,0.02)
dx24=np.random.uniform(0,0.02)
w2=np.array([dx1,dx2, dx3, dx4])

dx31=np.random.uniform(0,0.02)
dx32=np.random.uniform(0,0.02)
dx33=np.random.uniform(0,0.02)
dx34=np.random.uniform(0,0.02)
w3=np.array([dx1,dx2, dx3, dx4])

dict={}

for i in range(t):
    x1n,x2n, x3n, x4n=standard4d(K, B, x1,x2,x3,x4)
    #print(x2n)
    x1=x1n %1
    x2=x2n %1
    #print(x2)
    x3=x3n %1
    x4=x4n %1

    Mat4=np.array([[1+K*np.cos(2*np.pi*x1)+B*np.cos(2*np.pi*(x3-x1)), 1, -B*np.cos(2*np.pi*(x3-x1)), 0],
                  [K*np.cos(2* np.pi*x1)+B*np.cos(2*np.pi*(x1-x3)) , 1, -B*np.cos(2*np.pi*(x3-x1)), 0],
                  [-B*np.cos(2*np.pi*(x1-x3)), 0, 1+K*np.cos(2*np.pi*x3)+B*np.cos(2*np.pi*(x1-x3)), 1],
                  [-B*np.cos(2*np.pi*(x1-x3)), 0, K*np.cos(2*np.pi*x3)+B*np.cos(2*np.pi*(x1-x3)), 1]])

    dx1, dx2, dx3, dx4 =np.dot(Mat4, w0)
    dx11, dx12, dx13, dx14 =np.dot(Mat4, w1)
    dx21, dx22, dx23, dx24 =np.dot(Mat4, w2)
    dx31, dx32, dx33, dx34 =np.dot(Mat4, w3)

    w0=np.array([dx1,dx2, dx3, dx4])
    w1=np.array([dx11,dx12, dx13, dx14])
    w2=np.array([dx21,dx22, dx23, dx24])
    w3=np.array([dx31,dx32, dx33, dx34])

    #Gram-Schmidt Orthogonalization
    u0=w0
    a0= np.sqrt(dx1**2+dx2**2+dx3**2+dx4**2)
    w0= u0/a0
    alpha0.append(a0)

    u1=w1-np.inner(w1,w0)*w0
    a1= np.sqrt(dx11**2+dx12**2+dx13**2+dx14**2)
    w1= u1/a1
    alpha1.append(a1)

    u2=w2-np.inner(w2, w1)*w1-np.inner(w2,w0)*w0
    a2= np.sqrt(dx21**2+dx22**2+dx23**2+dx24**2)
    w2= u2/a2
    alpha2.append(a2)

    u3=w3- np.inner(w3,w2)*w2-np.inner(w3,w1)*w1-np.inner(w3,w0)*w0
    a3=np.sqrt(dx31**2+dx32**2+dx33**2+dx34**2)
    w3= u3/a3
    alpha3.append(a3)


lalpha0=np.log(alpha0)
lalpha1=np.log(alpha1)
lalpha2=np.log(alpha2)
lalpha3=np.log(alpha3)


#print(lalpha)

lyap0=[]
for i in range (1,len(lalpha0)):
    l0=sum (lalpha0[0:i])/((i))
    lyap0.append(l0)


lyap1=[]
for i in range (1,len(lalpha1)):
    l1=sum (lalpha1[0:i])/((i))
    lyap1.append(l1)

lyap2=[]
for i in range (1,len(lalpha2)):
    l2=sum (lalpha2[0:i])/((i))
    lyap2.append(l2)

lyap3=[]
for i in range (1,len(lalpha3)):
    l3=sum (lalpha3[0:i])/((i))
    lyap3.append(l3)

lyap0=lyap0[1:]
lyap1=lyap1[1:]
lyap2=lyap2[1:]
lyap3=lyap3[1:]
ll=np.arange(len(lyap0))

S=[]

for i in range (len(lyap0)):
    s=lyap0[i]+lyap1[i]+lyap2[i]+lyap3[i]
    S.append(s)
print(S)
sl= np.arange(len(S))

plt.plot((lyap0), label="Lyap1")
plt.plot((lyap1), label='Lyap2')
plt.plot((lyap2), label="Lyap3")
plt.plot((lyap3), label="Lyap4")
plt.ylim(-1,0.8)
plt.legend()
plt.ylabel('Lyapunov exponent')
plt.title("Lyapunov Spectrum of the 4D Standard Map with K= "+str(K)+", B="+str(B)+" x1= "+str(a)+", x2= "+str(b)+", x3= "+str(c)+", x4= "+str(d))
plt.show()

plt.plot(sl, S, "g", label='Sum')
plt.title("Symmetry check of Lyapunov Spectrum for the 4D Standard Map with K= "+str(K)+", B="+str(B)+" x1= "+str(a)+", x2= "+str(b)+", x3= "+str(c)+", x4= "+str(d))
plt.ylabel("Difference between LEs")
plt.xlabel("Iteration")
plt.show()

