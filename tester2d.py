import numpy as np
import matplotlib.pyplot as plt
from Maps4 import chirikov

plt.rcParams['figure.dpi'] = 200


def frame(K, N, r):
    '''generates data for N iterations of the Chirikov map with parameter e and 
    resolution r
    
    constant colour represents the same original point but after multiple iterations'''

    p0 = np.arange(0, 2*np.pi, r)
    q0 = np.arange(0, 2*np.pi, r)

    P, Q = np.meshgrid(p0, q0)

    P = P.reshape(len(p0)*len(q0))
    Q = Q.reshape(len(p0)*len(q0))
    
    # The original colour of a point is uniquely determined by an RGB value determined by
    # the starting point where Red is given by P0, Blue is given by Q0, and
    # Green is given by P0+Q0.
    colours = np.array([P/max(P), (P+Q)/max(P+Q), Q/max(Q)]).T

    Pnew = [P]
    Qnew = [Q]
    for i in range(N):
        P, Q = chirikov(P,Q, K=K)
        Qnew.append(Q)
        Pnew.append(P)
    return np.hstack(Pnew), np.hstack(Qnew), np.vstack([colours]*(N+1))

def plot(K, N, r):
    '''plots N iterations of the Chirikov map with parameter e and resolution r
    
    constant colour represents the same original point but after multiple iterations'''
    P, Q, colours = frame(K, N, r)

    plt.scatter(Q,P, s=.1, c=colours)
    
    plt.xlabel('Q')
    plt.ylabel('P')
    plt.title("2D Standard map with K = " + str(e) + ", N = " + str(N) + "")
    plt.show()

e=0.71

plot(K=e/2*np.pi, N=1000, r=.5)

