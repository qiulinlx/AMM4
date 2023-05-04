import Maps4 as mp
import numpy as np

#Henon-Heiles

def Henonvar(t,z):
    [v1, v2, v3, v4, x, y]=z
    dv1=v3
    dv2=v4
    dv3= -v1-2*x*v2-2*y*v1
    dv4= -v2-2*x*v1+2*y*v2 
    return [dv1, dv2, dv3, dv4, x,y]


def Standard2Var(K, p, q, dp, dq):
    Mat2=np.array([[1, K*np.cos(q)],[1, 1+K*np.cos(q)]]) %(2*np.pi)
    return np.dot(Mat2, np.array([dp,dq]))



def Standard4Var(K, B, dx1,dx2,dx3,dx4, x1, x2, x3, x4):
    Mat4=np.array([[1+K*np.cos(2*np.pi*x1)+B*np.cos(2*np.pi*(x3-x1)), 1, -B*np.cos(2*np.pi(x3-x1)), 0],
                  [-B*np.sin(2*np.pi*(x1-x3)) , 1, -B*np.sin(2*np.pi*(x3-x1)), 0],
                  [-B*np.sin(2*np.pi*(x1-x3)), 0, 1+K*np.cos(2*np.pi*x3)+B*np.cos(2*np.pi*(x1-x3)), 1],
                  [-B*np.sin(2*np.pi(x1-x3)), 0, K*np.cos(2*np.pi*x3)+B*np.cos(2*np.pi*(x1-x3)), 1]])
    return np.dot(Mat4, np.array([dx1,dx2,dx3,dx4]))
