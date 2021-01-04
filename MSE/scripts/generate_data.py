import numpy as np



def save_columns(x, x_name, y, y_name):
    np.savetxt(x_name, x)
    np.savetxt(y_name, y)
    return 0





if __name__ == "__main__":
    test = np.arange(0, 100, 1)
    np.savetxt("test.txt", test)