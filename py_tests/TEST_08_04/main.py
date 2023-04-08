import numpy as np
import scipy as sp
import  matplotlib.pyplot as plt
def create_test():
    N = 289
    n = 1
    A = 12
    B = 30
    C = 1
    fA = open("test_data.txt", "a")
    fi = open("test_i.txt", "a")
    fj = open("test_j.txt", "a")
    fb = open("test_b.txt", "a")
    fx = open("test_x.txt", "a")
    fL_max = open("test_L_max.txt", "a")
    fL_min = open("test_L_min.txt", "a")
    for i in range(n):
        matrix = []
        b = []
        for p in range(N):
            line = []
            for k in range(N):
                t = 0
                if k == p:
                    t = B * 2
                    line.append(t)
                    fA.write(str(t))
                    fi.write(str(p))
                    fj.write(str(k))
                    fA.write(" ")
                    fi.write(" ")
                    fj.write(" ")

                if p == k + 1 or  p  == k - 1:
                    t = A
                    line.append(t)
                    fA.write(str(t))
                    fi.write(str(p))
                    fj.write(str(k))
                    fA.write(" ")
                    fi.write(" ")
                    fj.write(" ")
                if p == k + 17 or p == k - 17:
                    t = A
                    line.append(t)
                    fA.write(str(t))
                    fi.write(str(p))
                    fj.write(str(k))
                    fA.write(" ")
                    fi.write(" ")
                    fj.write(" ")
                if t==0:
                    line.append(0)

            matrix.append(line)
        fA.write("\n")
        fi.write("\n")
        fj.write("\n")
        rank = np.linalg.matrix_rank(np.array(matrix))

        b = np.ones(N)
        b*=C

        lmbda, va = np.linalg.eigh(np.array(matrix))
        mmax = np.max(np.abs(lmbda))
        mmin = np.min(lmbda)
        x = np.linalg.solve(matrix, b)

        for p in b:
            fb.write(str(p))
            fb.write(" ")
        fb.write("\n")

        for p in x:
            fx.write(str(p))
            fx.write(" ")
        fx.write("\n")
        fL_max.write(str(mmax))
        fL_max.write("\n")
        fL_min.write(str(mmin))
        fL_min.write("\n")
        matrix.clear()
        b.clear()


#create_test()
def plot():
    plt.figure(figsize=(12, 7))
    plt.minorticks_on()
    plt.grid(
        which='major'
    )
    plt.grid(
        which='minor',
        linestyle='--'
    )

    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_MPI.txt", 'r') as f:
        y = f.readline()
        data = y.split()
        data_MPI = np.array(data, dtype=float)
        it_MPI = np.arange(len(data_MPI))
    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_MPI_opt.txt", 'r') as f:
        y = f.readline()
        data = y.split()
        data_MPI_opt = np.array(data, dtype=float)
        it_MPI_opt = np.arange(len(data_MPI_opt))
    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_CHEB_MPI.txt", 'r') as f:
        y = f.readline()
        data = y.split()
        data_CHEB_MPI = np.array(data, dtype=float)
        it_CHEB_MPI = np.arange(len(data_CHEB_MPI))
    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_SOR.txt", 'r') as f:
        y = f.readline()
        data = y.split()
        data_SOR = np.array(data, dtype=float)
        it_SOR = np.arange(len(data_SOR))
    plt.plot(it_MPI, np.log(data_MPI), label="MPI")
    plt.plot(it_MPI_opt, np.log(data_MPI_opt), label="MPI_opt")
    plt.plot(it_CHEB_MPI, np.log(data_CHEB_MPI), label="CHEB_MPI")
    plt.plot(it_SOR, np.log(data_SOR), label="SOR")
    plt.legend()
    plt.savefig("first_ex.png")
    plt.show()

#plot()
def plot_second():
    plt.figure(figsize=(12, 7))
    plt.minorticks_on()
    plt.grid(
        which='major'
    )
    plt.grid(
        which='minor',
        linestyle='--'
    )

    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_CHEB_MPI_second.txt", 'r') as f:
        data_MPI = []
        y = f.readlines()
        for line in y:
            data = line.split()
            data_MPI.append(np.array(data, dtype=float))

    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_count_CHEB_MPI_second.txt", 'r') as f:
        it_MPI = []
        y = f.readlines()
        for line in y:
            data = line.split()
            it_MPI.append(np.array(data, dtype=float))
    plt.plot(it_MPI[0], data_MPI[0], label = "delta = 96")
    plt.plot(it_MPI[1], data_MPI[1], label = "delta = 98")
    plt.plot(it_MPI[2], data_MPI[2], label = "delta = 100")
    plt.plot(it_MPI[3], data_MPI[3], label = "delta = 102")
    plt.plot(it_MPI[4], data_MPI[4], label = "delta = 104")
    plt.xlabel("iteration")
    plt.ylabel(r"$r$")
    plt.legend()
    plt.savefig("delta_lambda_1.png")
    plt.show()
#plot_second()
def second_part_plot():
    plt.figure(figsize=(12, 7))
    plt.minorticks_on()
    plt.grid(
        which='major'
    )
    plt.grid(
        which='minor',
        linestyle='--'
    )
    data_MPI_matrix = []
    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_MPI_ex2.txt", 'r') as f:
        y = f.readlines()
        for line in y:
            data = line.split()
            data_MPI_matrix.append(np.array(data, dtype=float))
        it_MPI = np.arange(len(data_MPI_matrix))

    data_MPI_opt_matrix = []
    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_MPI_opt_ex2.txt", 'r') as f:
        y = f.readlines()
        for line in y:
            data = line.split()
            data_MPI_opt_matrix.append(np.array(data, dtype=float))
        it_MPI_opt = np.arange(len(data_MPI_opt_matrix))
    data_CHEB_MPI_matrix = []
    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_CHEB_MPI_ex2.txt", 'r') as f:
        y = f.readlines()
        for line in y:
            data = line.split()
            data_CHEB_MPI_matrix.append(np.array(data, dtype=float))
        it_CHEB_MPI = np.arange(len(data_CHEB_MPI_matrix))
    data_Steepest_Descent_matrix = []
    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_Steepest_Descent_ex2.txt", 'r') as f:
        y = f.readlines()
        for line in y:
            data = line.split()
            data_Steepest_Descent_matrix.append(np.array(data, dtype=float))
        it_Steepest_Descent = np.arange(len(data_Steepest_Descent_matrix))
    data_Conjugate_Gradient_matrix = []
    with open("/home/ilya/slae_lab/slae/tests/TEST_08_04/test_r_Conjugate_Gradient_ex2.txt", 'r') as f:
        y = f.readlines()
        for line in y:
            data = line.split()
            data_Conjugate_Gradient_matrix.append(np.array(data, dtype=float))
        it_Conjugate_Gradient= np.arange(len(data_Conjugate_Gradient_matrix))

    lambda_max_v = {0, 0, 0, 1}
    lambda_min_v = {1, 0, 0, 0}
    x_MPI = []
    y_MPI = []
    x_MPI_opt = []
    y_MPI_opt = []
    x_CHEB_MPI = []
    y_CHEB_MPI = []
    x_Steepest_descent = []
    y_Steepest_descent = []
    x_Conjugate_Gradient = []
    y_Conjugate_Gradient = []
    for it in range(len(it_MPI)):
        x_MPI.append(data_MPI_matrix[it][0])
        y_MPI.append(data_MPI_matrix[it][3])

    for it in range(len(it_MPI_opt)):
        x_MPI_opt.append(data_MPI_opt_matrix[it][0])
        y_MPI_opt.append(data_MPI_opt_matrix[it][3])

    for it in range(len(it_CHEB_MPI)):
        x_CHEB_MPI.append(data_CHEB_MPI_matrix[it][0])
        y_CHEB_MPI.append(data_CHEB_MPI_matrix[it][3])

    for it in range(len(it_Steepest_Descent)):
        x_Steepest_descent.append(data_Steepest_Descent_matrix[it][0])
        y_Steepest_descent.append(data_Steepest_Descent_matrix[it][3])
    for it in range(len(it_Conjugate_Gradient)):
        x_Conjugate_Gradient.append(data_Conjugate_Gradient_matrix[it][0])
        y_Conjugate_Gradient.append(data_Conjugate_Gradient_matrix[it][3])


    plt.plot(x_MPI[1:], y_MPI[1:], "--o", label = "MPI")
    plt.plot(x_MPI_opt[1:], y_MPI_opt[1:],"--o", label = "MPI_opt")
    plt.plot(x_CHEB_MPI[1:], y_CHEB_MPI[1:],"--o", label = "CHEB_MPI")
    plt.plot(x_Steepest_descent[1:], y_Steepest_descent[1:],"--o", label = "Steepest_descent")
    plt.plot(x_Conjugate_Gradient[1:], y_Conjugate_Gradient[1:],"--o", label = "Conjugate_Gradient")
    plt.xlabel(r'$\lambda_{min}$')
    plt.ylabel(r'$\lambda_{max}$')
    plt.legend()
    plt.savefig("second_ex_with_mpi.png")
    plt.show()

second_part_plot()