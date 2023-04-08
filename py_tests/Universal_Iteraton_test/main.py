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
    fL_SSOR = open("test_L_SSOR.txt", "a")
    fW = open("test_w.txt", "a")
    fL_max = open("test_L_max.txt", "a")
    fL_min = open("test_L_min.txt", "a")
    fL_SYM_GS = open("test_L_SYM_GS.txt", "a")
    fL_CHEB_GS = open("test_L_CHEB_GS.txt", "a")
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
            for k in range(p, N):
                if k == p:
                    t = np.round(N * N * np.random.uniform(0.5, 1), 4)
                    D[k][p] = t
                else:
                    if np.random.random() < 0.1:
                        t = np.round(np.random.uniform(-np.sqrt(N), np.sqrt(N)), 4)
                        matrix[p][k] = t
                        if k > p :
                            U[p][k] = t
        matrix =np.array(matrix).T @ np.array(matrix)
        L = np.array(U).T
        matrix = matrix + np.array(D)
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
        mu_max = np.max(np.abs(np.linalg.eigvals(np.eye(N) - np.linalg.inv(D) @ matrix)))
        omega = 1 + (mu_max/(1 + np.sqrt(1 - mu_max ** 2))) ** 2
        P = np.linalg.inv(np.array(D) + omega * np.array(U)) @ (omega * np.array(L) + (omega - 1) * np.array(D)) @ np.linalg.inv(np.array(D) + omega * np.array(L)) @ (omega * np.array(U) + (omega - 1) * np.array(D))
        eig_vals= np.linalg.eigvals(P)
        fL_SSOR.write(str(np.max(np.abs(eig_vals))))
        fW.write(str(omega))
        fW.write("\n")
        fL_SSOR.write("\n")

        eig_vals = np.linalg.eigvals(matrix)
        fL_max.write(str(np.max(np.abs(eig_vals))))
        fL_min.write(str(np.min(np.abs(eig_vals))))
        fL_max.write("\n")
        fL_min.write("\n")
        P = np.linalg.inv(np.array(D) + np.array(U)) @ np.array(L) @ np.linalg.inv(
            np.array(D) + np.array(L)) @ np.array(U)
        eig_vals = np.linalg.eigvals(P)
        fL_CHEB_GS.write(str(np.max(np.abs(eig_vals))))
        fL_CHEB_GS.write("\n")

        mu_max = np.max(np.abs(np.linalg.eigvals(np.eye(N) - np.linalg.inv(D) @ matrix)))
        omega = 1 + (mu_max / (1 + np.sqrt(1 - mu_max ** 2))) ** 2
        P = - np.linalg.inv(np.array(D) + omega * np.array(L)) @ (omega * np.array(U) + (omega - 1) * np.array(D))
        eig_vals = np.linalg.eigvals(P)
        fL_SYM_GS.write(str(np.max(np.abs(eig_vals))))
        fL_SYM_GS.write("\n")

        b.clear()


create_test()
