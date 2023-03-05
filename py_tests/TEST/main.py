import matplotlib.pyplot as plt
import numpy as np


def plot(x, y, t_itle, name, x_label, y_label):
    plt.figure(figsize=(12, 7))
    plt.minorticks_on()
    plt.grid(
        which='major'
    )
    plt.grid(
        which='minor',
        linestyle='--'
    )
    plt.title(t_itle)
    plt.plot(x, y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(name)
    plt.show()

def number3():
    with open("/home/ilya/slae_lab/slae/tests/TEST_04_03/test3_data.txt", 'r') as f:
        y = f.readline()
        data = y.split()
        data = np.array(data, dtype=float)
        tau = 0.001
        x = []
        for i in range(len(data)):
            x.append(tau)
            tau+=0.001
        plot(x, data)
def MPI():
    with open("/home/ilya/slae_lab/slae/tests/TEST_04_03/test4_data_MPI.txt", 'r') as f:
        y = f.readline()
        data = y.split()
        data = np.array(data, dtype=float)
        #ind = np.where(data == 0)[0][0]
        x = np.arange(0, len(data))
        plot(x, np.log(data), "MPI, tau = 0.0001", "MPI_plot.png", "iteration", "ln(discrepancy)")

def Jacobi():
    with open("/home/ilya/slae_lab/slae/tests/TEST_04_03/test4_data_jacobi.txt", 'r') as f:
        y = f.readline()
        data = y.split()
        data = np.array(data, dtype=float)
        ind = np.where(data == 0)[0][0]
        x = np.arange(0, len(data[:ind]))
        plot(x, np.log(data[: ind]), "Jacobi", "jacobi_plot.png","iteration", "ln(discrepancy)")

def Gauss():
    with open("/home/ilya/slae_lab/slae/tests/TEST_04_03/test4_data_Gauss.txt", 'r') as f:
        y = f.readline()
        data = y.split()
        data = np.array(data, dtype=float)
        ind = np.where(data == 0)[0][0]
        x = np.arange(0, len(data[:ind]))
        plot(x, np.log(data[: ind]), "Gauss_Siedel", "Gauss_Siedel_plot.png","iteration", "ln(discrepancy)")

Gauss()