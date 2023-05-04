from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

E=1/6
x=0
y= 0.8 #np.random.uniform(-0.05, 0.05)
py= 0.12 #np.random.uniform(-0.05, 0.05)
px= np.sqrt(2*E-y**2+2/3*y**3-py**2)

tau=0.02

v1=0.01
v2=0.01
v3=0.01
v4=0.01

k=5000
""" t=tau*np.arange(0,k)

z0 = [v1,v2,v3,v4, x, y, px, py]
sol= solve_ivp(Henonvar, [0, k*tau],  z0 ,  t_eval=t, method="DOP853" )

v1=sol.y[0]
v2=sol.y[1]
v3=sol.y[2]
v4=sol.y[3]
x=sol.y[4]
y=sol.y[5]
px=sol.y[6]
py=sol.y[7] """

alpha=[]

#Initial position
for i in range(1,k):
    t= [tau*(i-1), i*tau]
   # print(state)
    q=np.array([x,y,px,py])
    state=  np.array([v1,v2,v3,v4, x,y])

    #Initial state
    Di= solve_ivp( Henonvar ,t, state, method="DOP853" ) #Evolving the orbit and deviation vectori
    Traj=solve_ivp(henonmotion, t, q, method="DOP853") #Evolving the orbit

    Di0= Di.y[0]
    Di1=Di.y[1]
    Di2=Di.y[2]
    Di3=Di.y[3]

    v1=Di0[-1]
    v2=Di1[-1]
    v3=Di2[-1]
    v4=Di3[-1]

    Trajx=Traj.y[0]
    Trajy=Traj.y[1]
    Trajpx=Traj.y[2]
    Trajpy=Traj.y[3]

    x=Trajx[-1]
    y=Trajy[-1]
    px=Trajpx[-1]
    py=Trajpy[-1]

    a=np.sqrt(v1**2+v2**2+v3**2+v4**2)
    #print(a)

    alpha.append(a)

    #Renormalization
    v1=v1/a
    v2=v2/a
    v3=v3/a
    v4=v4/a






#print(alpha)

#print(alpha)
alpha=alpha[1:]
lalpha=np.log(alpha)
#print(lalpha)

lyap=[]
for i in range (len(lalpha)):
    l=sum (lalpha[1:i])/(tau*(i+1))
    lyap.append(l)

print(lyap[-1])
lyap=lyap[1:]
ll=np.arange(len(lyap))


""" plt.plot(ll,lyap, 'b')
plt.xlabel('t')
plt.ylabel('Lyapunov exponent')
plt.title("Lyapunov exponent for the Henon-Heiles system with E= "+str(E)+" and tau= "+str(tau)+"")
plt.show() """