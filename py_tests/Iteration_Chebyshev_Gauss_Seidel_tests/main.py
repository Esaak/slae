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

    for i in range(n):
        matrix = [0] * N
        D = [0] * N
        U = [0] * N
        L = [0] * N
        for j in range(N):
            matrix[j] = [0] * N
            D[j] = [0] * N
            U[j] = [0] * N
            L[j] = [0] * N
        b = []
        for p in range(N):
            for k in range(N):
                if k == p:
                    t = np.round(N * N * np.random.uniform(0.5, 1), 4)
                    matrix[k][p] = t
                    D[k][p] = t
                else:
                    if np.random.random() < 0.1:
                        t = np.round(np.random.uniform(-N, N), 4)
                        matrix[p][k] = t
                        if k > p :
                            U[p][k] = t
                        else:
                            L[p][k] = t

        for p in range(N):
            for k in range(N):
                if matrix[p][k] != 0:
                    fA.write(str(matrix[p][k]))
                    fi.write(str(p))
                    fj.write(str(k))
                    fA.write(" ")
                    fi.write(" ")
                    fj.write(" ")
        fA.write("\n")
        fi.write("\n")
        fj.write("\n")
        rank = np.linalg.matrix_rank(np.array(matrix))

        for p in range(rank):
            b.append(np.random.uniform(-N, N))
        for p in range(n - rank):
            b.append(0)

        x = np.linalg.solve(matrix, b)

        for p in b:
            fb.write(str(p))
            fb.write(" ")
        fb.write("\n")

        for p in x:
            fx.write(str(p))
            fx.write(" ")
        fx.write("\n")
        P = np.linalg.inv(np.array(D) + np.array(U)) @ np.array(L) @ np.linalg.inv(np.array(D) + np.array(L)) @ np.array(U)
        eig_vals= np.linalg.eigvals(P)
        fL.write(str(np.max(np.abs(eig_vals))))
        fL.write("\n")

        b.clear()


create_test()
