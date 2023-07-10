# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Para la simulación de una reaccion con un sustrato A que pasa a ser el producto B con una probabilidad de P
#siendo la reacción reversible

#importamos las librerias que vamos a usar

import matplotlib.pyplot as plt
import numpy as np
import random

#creamos la matriz de eventos A B en un tiempo determinado y con un numero de simulaciones

def mat_simul_ab(p_ab, p_ba, n_eventos=100, t=100):
    mat = np.full((n_eventos, t), 'A')
    for i in range(n_eventos):
        for j in range(1,t):
            num = random.random()
            if mat[i,j-1]=='A':
                if num <= p_ab:
                    mat[i,j] = 'B'
                else:
                    mat[i,j] = 'A'
            elif mat[i,j-1]=='B':
                if num <= p_ba:
                    mat[i,j] = 'A'
                else:
                    mat[i,j] = 'B'
    return mat

#Definimos una función que cuente el número de apariciones de A y B en cada instante de tiempo
#Devueve una matriz cuya primera fila es el número de A y la segunda el numero de B

def contar_ab(p_ab, p_ba, n_eventos=100, t=100):
  mat = mat_simul_ab(p_ab, p_ba, n_eventos, t)
  count = np.zeros((2, t))
  for j in range(len(mat[0])):
    count[0, j] = np.count_nonzero(mat[:,j] == 'A')
    count[1, j] = np.count_nonzero(mat[:,j] == 'B')
  return(count)

#representamos esta proporción de reactivos y productos en un plot

matriz = contar_ab(0.02 , 0.01, 10000, 200)

plt.title('Evolucion de especies en el tiempo')
plt.plot(matriz[0], 'r', label='A')
plt.plot(matriz[1], 'b', label='B')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.legend(loc='best')
plt.show()