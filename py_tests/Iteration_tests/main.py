import numpy as np
import scipy as sp

def create_test():
    N = 100
    n = 10
    fA = open("test_data.txt", "a")
    fi = open("test_i.txt", "a")
    fj = open("test_j.txt", "a")
    fb = open("test_b.txt", "a")
    fx = open("test_x.txt", "a")
    fL = open("test_L.txt", "a")
    fL_min = open("test_l_min.txt", "a")
    for i in range(n):
        matrix = []
        b = []
        for p in range(N):
            line = []
            for k in range(N):
                if k == p:
                    t = np.round(N*N*np.random.uniform(0.5, 1), 4)
                    line.append(t)
                    fA.write(str(t))
                    fi.write(str(p))
                    fj.write(str(k))
                    fA.write(" ")
                    fi.write(" ")
                    fj.write(" ")
                else:
                    if np.random.random() < 0.1:
                        t = np.round(np.random.uniform(-N, N),4)
                        line.append(t)
                        fA.write(str(t))
                        fi.write(str(p))
                        fj.write(str(k))
                        fA.write(" ")
                        fi.write(" ")
                        fj.write(" ")
                    else:
                        line.append(0)
            matrix.append(line)
        fA.write("\n")
        fi.write("\n")
        fj.write("\n")
        rank = np.linalg.matrix_rank(np.array(matrix))

        for p in range(rank):
            b.append(np.random.uniform(-N, N))
        for p in range(n - rank):
            b.append(0)

        lmbda, va = np.linalg.eigh(matrix)
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
        fL.write(str(mmax))
        fL.write("\n")
        fL_min.write(str(mmin))
        fL_min.write("\n")
        matrix.clear()
        b.clear()

create_test()