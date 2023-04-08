import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

N = 10
low = 1
high = 10
dec = 3
def test(N, low, high, dec):
    a = [np.random.uniform(low, high) for _ in range(N-1)]
    c = [np.random.uniform(low, high) for _ in range(N-1)]
    a = np.round(a, dec)
    c = np.round(c, dec)

    b = np.round(np.abs(a) + np.abs(c) + np.abs(np.round(np.random.uniform(low, high), dec)) * 10, dec)
    b = np.round(np.append(b, np.abs(np.round(np.random.uniform(low, high), dec)) * 10), dec)
    d = np.round([np.random.uniform(low, high) for _ in range(N)], dec)

    A = np.diag(a, -1) + np.diag(b, 0) + np.diag(c, 1)
    x = np.round(np.linalg.linalg.solve(A, d),dec)
    return a, b, c, d, x


fA = open("testA.txt", "a")
fd = open("testd.txt", "a")
fx = open("testx.txt", "a")
for i in range(10):
    a, b, c, d, x = test(N, low, high, dec)
    fA.write(str(0))
    fA.write(' ')
    for i in range(N - 1):
        fA.write(str(a[i]))
        fA.write(' ')
    fA.write('\n')
    for i in range(N):
        fA.write(str(b[i]))
        fA.write(' ')
    fA.write('\n')
    for i in range(N - 1):
        fA.write(str(c[i]))
        fA.write(' ')
    fA.write(str(0))
    fA.write('\n')

    for i in range(N):
        fd.write(str(d[i]))
        fd.write(' ')
    fd.write('\n')

    for i in range(N):
        fx.write(str(x[i]))
        fx.write(' ')
    fx.write('\n')
fA.close()
fd.close()
fx.close()


