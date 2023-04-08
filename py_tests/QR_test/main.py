import scipy as sp
import numpy as np



def create_test():

    minn = 1
    maxx = 10
    N = 10
    dec = 3
    fM = open("matrix.txt", "w")
    fQ = open("Q.txt", "w")
    fR = open("R.txt", "w")

    #matrixs= np.random.rand(N,N)*(maxx-minn) + minn
    #matrixs = np.round(matrixs, dec)

    for p in range(N):
        matrix = np.random.rand(N, N) * (maxx - minn) + minn
        matrix = np.round(matrix, dec)

        for i in matrix:
            for j in i:
                fM.write(str(j))
                fM.write(" ")
            fM.write("\n")

        for i in np.linalg.qr(matrix)[0]:
            for j in i:
                fQ.write(str(j))
                fQ.write(" ")
            fQ.write("\n")

        for i in np.linalg.qr(matrix)[1]:
            for j in i:
                fR.write(str(j))
                fR.write(" ")
            fR.write("\n")
    fQ.close()
    fM.close()
    fR.close()

def simple_test():
    o = np.eye(3, dtype=float)
    q = np.arange(1, 10).reshape(3,3)
    print(np.array(np.linalg.qr(q)[0]))
    print(np.array(np.linalg.qr(q)[1]))
#create_test()
simple_test()