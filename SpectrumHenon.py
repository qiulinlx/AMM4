from Maps4 import *
from Variationalequations import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

E=1/12
x=0
a=-0.15
b=0.3
y= a #np.random.uniform(-0.05, 0.05)
py= b #np.random.uniform(-0.05, 0.05)
px= np.sqrt(2*E-y**2+2/3*y**3-py**2)

tau=0.05


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

for i in range(k):
#Initial state
    state0= np.concatenate((w0, [x,y]))
    state1= np.concatenate((w1, [x,y]))
    state2= np.concatenate((w2, [x,y]))
    state3= np.concatenate((w3, [x,y]))


    t= [tau*(i-1), i*tau]
    Di0 = solve_ivp( Henonvar ,t, state0, method="DOP853" ) #Evolving the orbit and deviation vector0
    Di1 =solve_ivp( Henonvar ,t, state1, method="DOP853" ) #Evolving the orbit and deviation vector1
    Di2 =solve_ivp( Henonvar ,t, state2, method="DOP853" ) #Evolving the orbit and deviation vector2
    Di3 =solve_ivp( Henonvar ,t, state3, method="DOP853" ) #Evolving the orbit and deviation vector3

    Traj=solve_ivp(henonmotion, t, q, method="DOP853") #Evolving the orbit

    Trajx=Traj.y[0]
    Trajy=Traj.y[1]
    Trajpx=Traj.y[2]
    Trajpy=Traj.y[3]

    x=Trajx[-1]
    y=Trajy[-1]
    px=Trajpx[-1]
    py=Trajpy[-1]

    q=np.array([x,y,px,py])


    w01= Di0.y[0]
    w02=Di0.y[1]
    w03=Di0.y[2]
    w04=Di0.y[3]

    v01=w01[-1]
    v02=w02[-1]
    v03=w03[-1]
    v04=w03[-1]

    w11= Di1.y[0]
    w12=Di1.y[1]
    w13=Di1.y[2]
    w14=Di1.y[3]

    v11=w11[-1]
    v12=w12[-1]
    v13=w13[-1]
    v14=w14[-1]

    w21= Di2.y[0]
    w22=Di2.y[1]
    w23=Di2.y[2]
    w24=Di2.y[3]

    v21=w21[-1]
    v22=w22[-1]
    v23=w23[-1]
    v24=w24[-1]

    w31= Di3.y[0]
    w32=Di3.y[1]
    w33=Di3.y[2]
    w34=Di3.y[3]

    v31=w31[-1]
    v32=w32[-1]
    v33=w33[-1]
    v34=w34[-1]


    w0 = np.array([v01,v02,v03,v04])
    w1 = np.array([v11,v12,v13,v14])
    w2 = np.array([v21,v22,v23,v24])
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


plt.plot((lyap0), label="Lyap1")
plt.plot((lyap1), label='Lyap2')
plt.plot((lyap2), label="Lyap3")
plt.plot((lyap3), label="Lyap4")
plt.ylim(-1,0.8)
plt.legend()
plt.ylabel('Lyapunov exponent')
plt.title("Lyapunov Spectrum for the Henon Heiles system with E= "+str(E)+", x1= "+str(a)+" and x2= "+str(b)+"")
plt.show()
plt.show()
