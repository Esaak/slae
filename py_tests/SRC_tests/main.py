import scipy as sp
import numpy as np
def create_test():
    N = 10
    minn = 1
    maxx = 10
    dim = 3
    fA = open("test_data.txt", "a")
    fd = open("test_indices.txt", "a")
    fx = open("test_indptr.txt", "a")
    fi = open("test_i.txt", "a")
    fj = open("test_j.txt", "a")
    fD = open("test_D.txt", "a")
    fX = open("test_X.txt", "a")
    for p in range(N):
        data = []
        D = []
        for i in range(N):
            line = []
            for j in range(N):
                temp = -1
                if np.random.rand() < 0.2:
                    temp = np.random.uniform(minn, maxx)
                    fi.write(str(i))
                    fj.write(str(j))
                    fi.write(' ')
                    fj.write(' ')
                else:
                    temp = 0
                line.append(np.round(temp, dim))
            D.append(np.random.uniform(minn, maxx))
            data.append(line)
        for it in sp.sparse.csr_matrix(data).data:
            fA.write(str(it))
            fA.write(' ')
        fA.write('\n')
        for it in sp.sparse.csr_matrix(data).indices:
            fd.write(str(it))
            fd.write(' ')
        fd.write('\n')
        for it in sp.sparse.csr_matrix(data).indptr:
            fx.write(str(it))
            fx.write(' ')
        for it in D:
            fD.write(str(it))
            fD.write(str(' '))

        for it in np.dot(data,D):
            fX.write(str(it))
            fX.write(str(' '))
        fX.write('\n')
        fD.write('\n')
        fx.write('\n')
        fi.write('\n')
        fj.write('\n')
    fD.close()
    fX.close()
    fA.close()
    fd.close()
    fx.close()

create_test()
# data = np.array( [8.943, 4.212, 6.965, 3.956, 6.448, 8.476, 1.279, 3.074, 2.49, 2.31, 3.842, 1.086, 3.214, 1.574, 3.556, 3.948])
# i =np.array( [0, 1, 2, 2, 2, 3, 4, 5, 5, 5, 6, 6, 7, 7, 9, 9])
# j = np.array([3, 6, 0, 8, 9, 6, 0, 0, 2, 8, 2, 3, 0, 6, 4, 9])
# arr = sp.sparse.csr_matrix((data, (i,j)), shape=(10,10)).toarray()
# o = np.ones(10)
# print(arr@o)