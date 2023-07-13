# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:29:35 2023

@author: celia
"""

#Simulando un sistema de transicion entre 3 estados: A, B, C
#El sistema puede transicionar de cualquier a cualquier estado con una probabilidad que se define en una lista

#Importamos las librerias necesarias
import numpy as np
from scipy.ndimage import gaussian_filter
import random


#Creamos una funcion que recoja en una matriz los estados de la simulacion en cada instante de tiempo y repetida un numero n de veces (n_eventos)
#Definimos la lista de probabilidades de transicion:
    #Lista de probabilidades: p_lista = [p_ab, p_ac, p_ba, p_bc, p_ca, p_cb]
#Asignamos un peso probabilistico a cada uno de los estados para modificar su porporcion inicial:
    #weights=(A, B, C)
#Añadimos la opcion de elegir el numero semilla para trabajar con la misma matriz
def mat_simul(p_lista, weights=(1/3, 1/3, 1/3), n_eventos=100, t=100, seed=8462836):
  mat = np.full((n_eventos, t), 'A')
  estados = ['A', 'B', 'C']

  if seed != 8462836:
    random.seed(seed)
  mat[:, 0] = random.choices(estados, weights=weights, k=n_eventos)

  for i in range(n_eventos):
    for j in range(1,t):
      x = random.random()

      if mat[i,j-1]=='A':
        if x > (p_lista[0]+p_lista[1]):
          mat[i,j] = 'A'
        elif x <= p_lista[0]:
          mat[i,j] = 'B'
        else:
          mat[i,j] = 'C'

      elif mat[i,j-1]=='B':
        if x > (p_lista[2]+p_lista[3]):
          mat[i,j] = 'B'
        elif x <= p_lista[2]:
          mat[i,j] = 'A'
        else:
          mat[i,j] = 'C'

      elif mat[i,j-1]=='C':
        if x > (p_lista[4]+p_lista[5]):
          mat[i,j] = 'C'
        elif x <= p_lista[4]:
          mat[i,j] = 'A'
        else:
          mat[i,j] = 'B'

  return mat


#Definimos otra funcion que cuente la cantidad de eventos que se encuentra en cada uno de los estados A, B, C en cada instante de tiempo
#Devuelve una matriz de tres filas correspondientes a los estados A, B y C respectivamente, siendo las columnas los instantes de tiempo
def count_abc(mat):
  t = len(mat[0])
  count = np.zeros((3, t))
  for j in range(len(mat[0])):
    count[0, j] = np.count_nonzero(mat[:,j] == 'A')
    count[1, j] = np.count_nonzero(mat[:,j] == 'B')
    count[2, j] = np.count_nonzero(mat[:,j] == 'C')
  return(count)


#Calculo de la entropia de Shannon
#Definimos una funcion que calcule la entropia en cada instante de tiempo, siendo la suma de las entropias de los tres estados
#Añadimos el parametro sigma que realiza un suavizado gaussiano sobre esta entropia
def entrop_simul(mat, n_eventos, sigma=3):
  t = len(mat[0])
  St = []
  for i in range(t):
    S = 0
    for j in range(len(mat[:,0])):
      if mat[j,i] != 0:
        Pi = mat[j,i]/n_eventos
        S = S - Pi*np.log(Pi)
    St.append(S)

  St = gaussian_filter(St,sigma=sigma)
  return St