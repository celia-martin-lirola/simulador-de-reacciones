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

def contar_ab(mat):
  count = np.zeros((2, len(mat[0])))
  for j in range(len(mat[0])):
    count[0, j] = np.count_nonzero(mat[:,j] == 'A')
    count[1, j] = np.count_nonzero(mat[:,j] == 'B')
  return(count)

#representamos esta proporción de reactivos y productos en un plot

p_ab = 0.02
p_ba = 0.01
n_eventos = 10000
t = 200

matriz_ab = mat_simul_ab(p_ab, p_ba, n_eventos, t)
print(matriz_ab)

matriz_conteo = contar_ab(matriz_ab)
print(matriz_conteo)

plt.title('Evolucion de especies en el tiempo')
plt.plot(matriz_conteo[0], 'r', label='A')
plt.plot(matriz_conteo[1], 'b', label='B')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.legend(loc='best')

plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/Plots/evolucion_ab.jpg')
plt.show()

#para comprobar que la tendencia de la simulacion es correcta, representamos la relacion Nb/Na
#esta relacion debe quedar en torno a la division de las probabilidades (p_ab/p_ba)

y = matriz_conteo[1]/matriz_conteo[0]
k = p_ab/p_ba

plt.plot(y, 'b', label = 'Nb/Na')
plt.axhline(y = k, color = 'r', label = 'p_ab/p_ba')
plt.title('Tendencia de la simulacion')
plt.legend(loc = 'best')
plt.xlabel('Tiempo')

plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/Plots/tendencia_ab.jpg')
plt.show()