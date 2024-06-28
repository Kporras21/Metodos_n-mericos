import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import legendre

def function(x):
    return x**4 + np.sin(x)**2

def I_Riemman(a, b, N):
    
    h = (b - a) / N
    integral = 0.0
    
    for i in range(N):
        x_i = a + i * h
        integral += function(x_i) * h
    return integral

#I1 = I_Riemman(-1, 1, 10)
#print(I1)

#I2 = I_Riemman(-1, 1, 100)
#print(I2)

#I3 = I_Riemman(-1, 1, 200)
#print(I3)

def I_Trapezoides(a, b, N):

    Sum = 0
    h = (b - a) / N

    for k in range(N):
        Sum += function(a+(k*h))
    return h * (1/2 * function(a) + 1/2 * function(b) + Sum)

#I4 = I_Trapezoides(-1, 1, 10)
#print(I4)

#I5 = I_Trapezoides(-1, 1, 100)
#print(I5)

#I6 = I_Trapezoides(-1, 1, 200)
#print(I6)

def I_Simpson(a, b, N):

    Sum_impar = 0
    Sum_par = 0
    h = (b - a) / N

    for k in range(1, N, 2):
        Sum_impar += function(a + (k * h))

    for p in range(2, N-1, 2):
        Sum_par += function(a + (p * h))

    return (h/3) * (function(a) + function(b) + 4 * Sum_impar + 2 * Sum_par)

#I7 = I_Simpson(-1, 1, 10)
#print(I7)

#I8 = I_Simpson(-1, 1, 100)
#print(I8)

#I9 = I_Simpson(-1, 1, 200)
#print(I9)

def gaussxw(N):

    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))
    epsilon = 1e-15
    delta = 1.0

    while delta > epsilon:
        p0 = np.ones(N, dtype = float)
        p1 = np.copy(x)

        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = np.max(np.abs(dx))

    w = 2 * (N + 1) * (N + 1)/(N * N * (1 - x * x) * dp * dp)

    return x,w

def gaussxwab(a, b, x, w):

    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

x2, w2 = gaussxw(2)

x3, w3 = gaussxw(3)

x4, w4 = gaussxw(4)

x2r , w2r = gaussxwab(-1, 1, x2, w2)

x3r , w3r = gaussxwab(-1, 1, x3, w3)

x4r , w4r = gaussxwab(-1, 1, x4, w4)

#I10 = np.sum(w2r * function(x2r))
#print(I10)

#I11 = np.sum(w3r * function(x3r))
#print(I11)

#I12 = np.sum(w4r * function(x4r))
#print(I12)

def error(I_N, I_an = 0.945351):
    return (I_N - I_an) / I_an

err_Riemman = np.zeros(100)

for N in range(1, 101):
    I_N = I_Riemman(-1, 1, N)
    err = error(I_N)
    err_Riemman[N-1] = err

#print(err_Riemman)
log_err_Riemman = np.log(err_Riemman)
#print(log_err_Riemman)

err_Trapezoides = np.zeros(100)

for N in range(1, 101):
    I_N = I_Trapezoides(-1, 1, N)
    err = error(I_N)
    err_Trapezoides[N-1] = err

log_err_Trapezoides = np.log(err_Trapezoides)

err_Simpson = np.zeros(100)

for N in range(1, 101):
    I_N = I_Simpson(-1, 1, N)
    err = error(I_N)
    err_Simpson[N-1] = err

log_err_Simpson = np.log(err_Simpson)

plt.plot(range(1, 101), log_err_Riemman, linestyle='-', color='b', label='Riemann')
plt.plot(range(1, 101), log_err_Trapezoides, linestyle='-', color='r', label='Trapezoides')
#plt.plot(range(1, 101), log_err_Simpson, linestyle='-', color='y', label='Simpson')
plt.xlabel('N')
plt.ylabel('Error')
plt.show()
