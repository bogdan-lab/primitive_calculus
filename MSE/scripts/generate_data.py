import numpy as np
import matplotlib.pyplot as plt


def save_columns(x, x_name, y, y_name):
    np.savetxt(x_name, x)
    np.savetxt(y_name, y)
    return 0





if __name__ == "__main__":
    x = np.linspace(-1, 1, 200)
    y = 1 + x + x**2 + np.sin(x)
    save_columns(x, "x.txt", y, "y.txt")
    plt.plot(x, y)
    