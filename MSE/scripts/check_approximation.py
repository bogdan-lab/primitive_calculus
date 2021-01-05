import numpy as np
import matplotlib.pyplot as plt
import os


def polynomial(x, weights):
    res = 0
    for i in range(len(weights)):
        res += x**i*weights[i]
    return res

def legandre_series(x, weights):
    res = 0
    lg_val = []
    for i in range(len(weights)):
        if i==0:
            res += weights[0]
            lg_val.append(1)
            continue
        if i==1:
            res += weights[1]*x
            lg_val.append(x)
            continue
        tmp = 2*(i-1)/i*x*lg_val[i-1] - (i-1)/i*lg_val[i-2]
        lg_val.append(tmp)
        res += weights[i]*tmp
    return res


def plot_approximation(x, y, weights, func):
    y_pts = np.zeros(len(x))
    for i in range(len(y_pts)):
        y_pts[i] = func(x[i], weights)
    plt.plot(x, y_pts, label="SE = %.3e, %s" % (sum((y-y_pts)**2), func.__name__))
    return 0


if __name__ == "__main__":
    
    os.chdir("../program_runs/")
    
    x = np.loadtxt("x.txt")
    y = np.loadtxt("y.txt")
    w_p = np.loadtxt("weight_pol.txt")
    w_l = np.loadtxt("weight_leg.txt")
    
    
    plt.figure()
    plt.grid()
    plt.plot(x, y, '.')
    plot_approximation(x, y, w_p, polynomial)
    plot_approximation(x, y, w_l, legandre_series)
    plt.legend()