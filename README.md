# Calculo-Numerico
Codigo das materia de Calculo Numerico

import numpy as np
import matplotlib.pyplot as plt


def thomas(A, B, C, D):
    B, D = B.copy(), D.copy()
    n = len(D)
    y = np.zeros(n)
    for i in range(1, n):
        if abs(B[i-1]) < 1e-14:
            raise ValueError('Pivo nulo ou muito pequeno.')
        m = A[i]/B[i-1]
        B[i] -= m*C[i-1]
        D[i] -= m*D[i-1]
    if abs(B[-1]) < 1e-14:
        raise ValueError('Sistema singular ou pivo muito pequeno.')
    y[-1] = D[-1]/B[-1]
    for i in range(n-2, -1, -1):
        y[i] = (D[i]-C[i]*y[i+1])/B[i]
    return y


def PVC(p, q, r, a, b, N, esquerda, direita):
    if N < 1 or b <= a:
        raise ValueError('Intervalo ou numero de nos invalidos.')
    te, alfa = esquerda
    td, beta = direita
    if te not in ('D', 'N') or td not in ('D', 'N'):
        raise ValueError('Use D ou N nas condicoes de contorno.')
    h = (b-a)/(N+1)
    x = np.linspace(a, b, N+2)
    A = np.zeros(N+2)
    B = np.zeros(N+2)
    C = np.zeros(N+2)
    D = np.zeros(N+2)
    for i in range(N+2):
        A[i] = 1+h*p(x[i])/2
        B[i] = -(2+h**2*q(x[i]))
        C[i] = 1-h*p(x[i])/2
        D[i] = h**2*r(x[i])
    if te == 'D':
        B[0], C[0], D[0] = 1, 0, alfa
    else:
        D[0] += 2*h*alfa*A[0]
        C[0] = 2
    A[0] = 0
    if td == 'D':
        A[-1], B[-1], D[-1] = 0, 1, beta
    else:
        D[-1] -= 2*h*beta*C[-1]
        A[-1] = 2
    C[-1] = 0
    y = thomas(A, B, C, D)
    return x, y


def p(x):
    return 0


def q(x):
    return 1


def r(x):
    return 0


if __name__ == '__main__':
    a, b = 0, 1
    N = 3
    esquerda = ('N', 1)
    direita = ('D', np.exp(1))
    exata = np.exp

    x, y = PVC(p, q, r, a, b, N, esquerda, direita)
    print('x          Numerico       Exato          Erro')
    for xi, yi in zip(x, y):
        print(f'{xi:.6f}   {yi:.8f}     {exata(xi):.8f}     {abs(yi-exata(xi)):.8f}')
    for n in (N, 2*N+1, 4*N+3):
        xx, yy = PVC(p, q, r, a, b, n, esquerda, direita)
        print(f'h = {(b-a)/(n+1):.6f}  Erro maximo = {max(abs(yy-exata(xx))):.8e}')
    plt.plot(x, y, 'o-', label='Diferencas finitas')
    xx = np.linspace(a, b, 200)
    plt.plot(xx, exata(xx), '--', label='Solucao exata')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()