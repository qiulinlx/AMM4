from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

E=0.167
x=0
y= 0.3 #np.random.uniform(-0.05, 0.05)
py= 0.0 #np.random.uniform(-0.05, 0.05)
px= np.sqrt(2*E-y**2+2/3*y**3-py**2)

n=4000
tau=0.02

alpha=[]

v1=0.01
v2=0.01
v3=0.01
v4=0.01

dx=v1
dy=v2
dpx=v3
dpy=v4

tau=0.1
c1=0.5*(1-1/np.sqrt(3))
c2= np.sqrt(3)/3
d1=0.5
cc=(2-np.sqrt(3))/24

xl=[]
xl.append(x)
yl=[]
yl.append(y)
pyl=[]
pyl.append(py)
pxl=[]
pxl.append(px)
dxl=[]
dxl.append(v1)
dyl=[]
dyl.append(v2)
dpxl=[]
dpxl.append(v3)
dpyl=[]
dpyl.append(v4) 


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


for i in range(n):
    z=eC(cc, x, y, px, py, dx, dy, dpx, dpy)
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
    dx=z[4]
    dy=z[5]
    dpx=z[6]
    dpy=z[7]
    
    a=np.sqrt(dx**2+dy**2+dpx**2+dpy**2)
    alpha.append(a)

    dx=dx/a
    dy=dy/a
    dpx=dpx/a
    dpy=dpy/a

    xl.append(x)
    yl.append(y)
    pxl.append(px)
    pyl.append(py)
    dxl.append(dx)
    dyl.append(dy)
    dpxl.append(dpx)
    dpyl.append(dpy)

#print(alpha)

alpha=alpha[1:]
lalpha=np.log(alpha)
#print(lalpha)

lyap=[]
for i in range (1,len(lalpha)):
    l=sum (lalpha[0:i])/(i+1)
    lyap.append(l)

#print(lyap)

lyap=lyap[1:]
ll=np.arange(len(lyap))
ll=np.arange(len(lyap))
plt.plot(ll,np.log10(lyap))
plt.xlabel('t')
plt.ylabel('Lyapunov exponent')
plt.title("Lyapunov exponent for the Henon-Heiles system with E= "+str(E)+" and tau= "+str(tau)+"")
plt.show()

#plt.plot(yl, pyl)
#plt.show() 