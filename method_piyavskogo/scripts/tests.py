import numpy as np
import matplotlib.pyplot as plt


def func(x):
    return abs(x)*np.sin(x)

def minorant(x, y, L):
    return func(y) - L*abs(x - y) 


def get_G_func(func_res):
    data = []
    for j in range(len(func_res[0])):
        max_val = -np.inf
        for i in range(len(func_res)):
            max_val = max(max_val, func_res[i][j][1])
        data.append([func_res[0][j][0], max_val])
    return np.array(data)



if __name__ == "__main__":
    

    
    L = 5
    check_pts = [-6, -4, -2, 0, 2, 4, 6]
    pts = np.linspace(-2*np.pi, 2*np.pi, 100)
    
    
    
    func_res = []
    for y in check_pts:
        func_res.append(np.column_stack((pts, minorant(pts, y, L))))
    
    G = get_G_func(func_res)
    
    
    
    
    plt.figure()
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.plot(pts, func(pts))
    #for el in func_res:
    #    plt.plot(el[:,0], el[:,1], "tab:orange")
    plt.plot(G[:,0], G[:,1], "tab:green")
    plt.show()


