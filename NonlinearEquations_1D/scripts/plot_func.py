import matplotlib.pyplot as plt
import numpy as np

def func(x):
    if (x==0):
        return 0;
    return x + x*x*np.sin(2/x)

def derivative(func, x, delta=1e-10):
    return (func(x+delta)-func(x-delta))/(2*delta)

if __name__ == "__main__":
        
    pts = np.linspace(-10, 10, 50000)
    vals = np.zeros(len(pts))
    der = np.zeros(len(pts))
    for i in range(len(vals)):
        vals[i] = func(pts[i])
        der[i] = derivative(func, pts[i])
    
    plt.figure()
    plt.grid()
#    plt.plot(pts, vals, label="func")
    plt.plot(pts, der, label="deriv")
    plt.legend()
    plt.show()