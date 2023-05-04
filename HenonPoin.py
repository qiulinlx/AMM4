import scipy as sp
import numpy as np
import matplotlib.pyplot as plt


def hamiltonian(x, y, px, py): #Henon
    h=0.5*(px**2+py**2)+0.5*(x**2+y**2)+x**2*y-(y**3)/3    
    return (h)

tau=0.02
c1=0.5*(1-1/np.sqrt(3))
c2= np.sqrt(3)/3
d1=0.5
cc=(2-np.sqrt(3))/24

E=0.167
x=0
y= 0.3#np.random.uniform(-0.05, 0.05)
py= 0.0 #np.random.uniform(-0.05, 0.05)
px= np.sqrt(2*E-y**2+2/3*y**3-py**2)


def f(x,y):
    return -x-2*x*y

def g(x,y):
    return -y-x**2+y**2 


def eA(c, x, y, px, py):
    x=x+tau*c*px
    y=y+tau*c*py
    px=px
    py=py
    return (x,y,px,py)

def eB(d, x, y, px, py):
    x=x
    y=y
    px=px+tau*d*f(x,y)
    py=py+tau*d*g(x,y)
    return (x,y,px,py)

def eC(cc, x,y,px,py):
    x=x
    y=y
    px=px+(2*x*(1+2*x**2+6*y+2*y**2))*(tau**3*cc/2)
    py=py+(2*(y-3*y**2+2*y**3+3*x**2+2*x**2*y))*(tau**3*cc//2)
    return (x,y,px,py) 
    

xl=[]
xl.append(x)
yl=[]
yl.append(y)
pyl=[]
pyl.append(py)
pxl=[]
pxl.append(px)

error=[]

for i in range(4000):

    z=eC(cc, x, y, px, py)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]

    z=eA(c1,x,y, px,py)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]

    z=eB(d1, x,y,px,py)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]

    z=eA(c2, x,y,px,py)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]

    z= eB(d1, x,y,px,py)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]

    z=eA(c1,x,y,px,py)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]

    z=eC(cc, x, y, px, py)
    x=z[0]
    y=z[1]
    px=z[2]
    py=z[3]

    xl.append(x)
    yl.append(y)
    pxl.append(px)
    pyl.append(py)

    e=hamiltonian(x, y, px, py)
    ee= np.abs(E-e)/E
    error.append(ee)

plt.plot(yl, pyl, 'o', markersize=1)
plt.xlabel('y')
plt.ylabel('py')
plt.title('Poincare section of Henon Heiles system with E=0.167')
plt.show()