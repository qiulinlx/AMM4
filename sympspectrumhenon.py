from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

E=1/8
x=0
a=0.5
b=-0.1
y= a #np.random.uniform(-0.05, 0.05)
py= b #np.random.uniform(-0.05, 0.05)
px= np.sqrt(2*E-y**2+2/3*y**3-py**2)

tau=0.1
c1=0.5*(1-1/np.sqrt(3))
c2= np.sqrt(3)/3
d1=0.5
cc=(2-np.sqrt(3))/24

k=2000
dict={}
q=np.array([x,y,px,py])

#Create 2N deviation vectors
for p in range(8):
    v1=np.random.uniform(0,0.05)
    v2=np.random.uniform(0,0.05)
    v3=np.random.uniform(0,0.05)
    v4=np.random.uniform(0,0.05)
    var="dict"+str(p)

    dict[var]=[v1,v2,v3,v4]

w0=dict["dict0"]
w1=dict["dict1"]
w2=dict["dict2"]
w3=dict["dict3"]

alpha0=[]
alpha1=[]
alpha2=[]
alpha3=[]

def eA(c, x, y, px, py, dx, dy, dpx, dpy):
    x=x+tau*c*px
    y=y+tau*c*py
    px=px
    py=py
    dx=dx+dpx*tau*c
    dy=dy+dpy*tau*c
    dpx=dpx
    dpy=dpy
    return (x, y, px, py, dx, dy, dpx, dpy)

def eB(d, x, y, px, py, dx, dy, dpx, dpy):
    x=x
    y=y
    px=px-tau*d*x*(1+2*y)
    py=py+tau*d*(y**2-x**2-y)
    dx=dx 
    dy=dy
    dpx=dpx-tau*d*((1+2*y)*dx+2*x*dy)
    dpy= dpy+d*tau*(-2*x*dx+(-1+2*y)*dy)
    return (x,y,px,py, dx, dy, dpx, dpy)

def eC(cc, x,y, px, py, dx, dy, dpx, dpy):
    x=x
    y=y
    px=px-2*x*(1+2*x**2+6*y+2*y**2)*cc*tau
    py=py-2*(y-3*y**2+2*y**3+3*x**2+2*x**2*y)*cc*tau
    dx=dx
    dy=dy
    dpx=dpx-2*((1+6*x**2+2*y**2+6*y)*dx+2*x*(3+2*y)*dy)*cc*tau
    dpy=dpy-2*(2*x*(3+2*y)*dx+(1+2*x**2+6*y**2-6*y)*dy)*cc*tau 
    return (x,y,px,py, dx, dy, dpx, dpy)


for i in range(k):
#Initial state
    z=eC(cc, x, y, px, py, w0[0], w0[1], w0[2], w0[3])
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]



    z=eA(c1,x,y, px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eB(d1, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eA(c2, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z= eB(d1, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eA(c1,x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    
    z=eC(cc, x, y, px, py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    v01=z[4]
    v02=z[5]
    v03=z[6]
    v04=z[7]


    w0 = np.array([v01,v02,v03,v04])

    z=eC(cc, x, y, px, py, w1[0], w1[1], w1[2], w1[3])
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]

    z=eA(c1,x,y, px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eB(d1, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eA(c2, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z= eB(d1, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eA(c1,x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    
    z=eC(cc, x, y, px, py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    v11=z[4]
    v12=z[5]
    v13=z[6]
    v14=z[7]

    w1 = np.array([v11,v12,v13,v14])

    z=eC(cc, x, y, px, py, w2[0], w2[1], w2[2], w2[3])
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]

    z=eA(c1,x,y, px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eB(d1, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eA(c2, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z= eB(d1, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    

    z=eA(c1,x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    
    z=eC(cc, x, y, px, py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    v21=z[4]
    v22=z[5]
    v23=z[6]
    v24=z[7]
    w2 = np.array([v21,v22,v23,v24])

    z=eC(cc, x, y, px, py, w3[0], w3[1], w3[2], w3[3])
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]

    z=eA(c1,x,y, px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]

    z=eB(d1, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]

    z=eA(c2, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]

    z= eB(d1, x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]

    z=eA(c1,x,y,px,py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    
    z=eC(cc, x, y, px, py, dx, dy, dpx, dpy)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]
    v31=z[4]
    v32=z[5]
    v33=z[6]
    v34=z[7]
    w3 = np.array([v31,v32,v33,v34])


#Gram-Schmidt Orthogonalization
    u0=w0
    a0= np.sqrt(v01**2+v02**2+v03**2+v04**2)
    w0= u0/a0
    alpha0.append(a0)

    u1=w1-np.inner(w1,w0)*w0
    a1= np.sqrt(v11**2+v12**2+v13**2+v14**2)
    w1= u1/a1
    alpha1.append(a1)

    u2=w2-np.inner(w2, w1)*w1-np.inner(w2,w0)*w0
    a2= np.sqrt(v21**2+v22**2+v23**2+v24**2)
    w2= u2/a2
    alpha2.append(a2)

    u3=w3- np.inner(w3,w2)*w2-np.inner(w3,w1)*w1-np.inner(w3,w0)*w0
    a3=np.sqrt(v31**2+v32**2+v33**2+v34**2)
    w3= u3/a3
    alpha3.append(a3)


""" print(alpha1, 'Break')
print(alpha3) """

""" plt.plot(alpha0, label='alpha0')
plt.plot(alpha1, label='alpha1')
plt.plot(alpha2, label='alpha2')
plt.plot(alpha3, label='alpha3')
plt.legend()
plt.show()
 """

alpha0=alpha0[1:]
alpha1=alpha1[1:]
alpha2=alpha2[1:]
alpha3=alpha3[1:]


lalpha0=np.log(alpha0)
lalpha1=np.log(alpha1)
lalpha2=np.log(alpha2)
lalpha3=np.log(alpha3)


#print(lalpha)

lyap0=[]
for i in range (1,len(lalpha0)):
    l0=sum (lalpha0[0:i])/(tau*(i))
    lyap0.append(l0)


lyap1=[]
for i in range (1,len(lalpha1)):
    l1=sum (lalpha1[0:i])/(tau*(i))
    lyap1.append(l1)

lyap2=[]
for i in range (1,len(lalpha2)):
    l2=sum (lalpha2[0:i])/(tau*(i))
    lyap2.append(l2)

lyap3=[]
for i in range (1,len(lalpha3)):
    l3=sum (lalpha3[0:i])/(tau*(i))
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
plt.title("Lyapunov Spectrum for the Henon Heiles system with E= "+str(E)+", x1= "+str(a)+" and x2= "+str(b)+"")
plt.show()

plt.plot(sl, S, label='Sum')
plt.title("Symmetry check of Lyapunov Spectrum for the Henon Heiles system with E= "+str(E)+", x1= "+str(a)+" and x2= "+str(b)+"")
plt.ylabel("Difference between LEs")
plt.xlabel("Iteration")
plt.show()
