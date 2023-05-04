import numpy as np

def henon(x, y, px, py): #Henon
    h=0.5*(px**2+py**2)+0.5*(x**2+y**2)+x**2*y-(y**3)/3    
    return (h)

def henonmotion(t,z):
    [x,y,px,py]=z
    dx = px
    dy = py
    dpx = -x - 2*x*y
    dpy = -y - x**2 + y**2
    return [dx, dy, dpx, dpy]

def chirikov(p, q, K):
    """Applies Chirikov's standard map once with parameter e"""
    
    p_prime = (p + (K)*np.sin(q)) % (2*np.pi)
    q_prime = (p + q + (K)*np.sin(q)) % (2*np.pi)
    
    return p_prime, q_prime

def standard4d(K, B, x1,x2,x3,x4):
    x4n=x4+(K/(2*np.pi))*(np.sin(2*np.pi*x3))-(B/(2*np.pi))*np.sin((2*np.pi*(x1-x3))) % (1)
    x3n= x3+x4+(K/(2*np.pi))*(np.sin(2*np.pi*x3))-(B/(2*np.pi))*np.sin(2*np.pi*(x1-x3))  % (1)
    x2n=x2+(K/(2*np.pi))*(np.sin(2*np.pi*x1))-(B/(2*np.pi)) *np.sin(2*np.pi*(x3-x1)) % (1)
    x1n=x1+x2+(K/(2*np.pi))*(np.sin(2*np.pi*x1))-(B/(2*np.pi))*np.sin(2*np.pi*(x3-x1)) % (1)
    return x1n,x2n,x3n, x4n

