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


#Calculo de la variacion de entropia de Shannon
#Devuelve una lista con los valores de variacion segun el tiempo
def var_entrop(mat):
  delta_St = []
  for x in range(1, len(mat)):
    delta_S = mat[x]-mat[x-1]
    delta_St.append(delta_S)
  return(delta_St)


#Para la divergencia de Kullback-Leibler
#Calculo de la entropia relativa o divergencia de informacion, propocional a la energia que debe invertir el sistema para lograr cada estado
#Realizamos un conteo de transiciones de cualquier a cualquier estado, distinguiendo entre cada una y recogiendolo en una matriz
def contar_trans(mat):
  t = len(mat[0])
  mat_trans = np.zeros((6, t-1))
  for j in range(1, len(mat[0])):
    for i in range(len(mat[:, 0])):
      if mat[i, j-1]=='A' and mat[i,j-1]!=mat[i,j]:
        if mat[i,j]=='B':
          mat_trans[0,j-1] += 1
        elif mat[i,j]=='C':
          mat_trans[2,j-1] += 1
      elif mat[i, j-1]=='B' and mat[i,j-1]!=mat[i,j]:
        if mat[i,j]=='A':
          mat_trans[1,j-1] += 1
        elif mat[i,j]=='C':
          mat_trans[4,j-1] += 1
      elif mat[i, j-1]=='C' and mat[i,j-1]!=mat[i,j]:
        if mat[i,j]=='A':
          mat_trans[3,j-1] += 1
        elif mat[i,j]=='B':
          mat_trans[5,j-1] += 1
  return (mat_trans)


#Presentamos varias maneras de normalizar estas transiciones:
    
#Normalizando el numero de transiciones dividiendo entre el total
def trans_norm(mat):
  mat = contar_trans(mat)
  for i in range(len(mat[0])):
    tot = np.sum(mat[:,i])
    for j in range(len(mat[:,0])):
      if mat[j,i] != 0:
        mat[j,i] = mat[j,i]/tot
  return(mat)

#Normalizando por tandas de intervalos de tiempo
#Se debe especificar el intervalo deseado
#Se divide entre el total de los intervalos sin solapar
def trans_norm_interv(mat, intervalo=5):
  mat = contar_trans(mat)
  t = len(mat[0])
  n = t//intervalo
  mat_norm = np.zeros((6,n))
  for x in range(n):    #normalizar transiciones
    tot = np.sum(mat[:,(x*intervalo):((x+1)*intervalo)])
    if tot != 0:
      for j in range(len(mat[:,0])):
        trans_j = sum(mat[j,(x*intervalo):((x+1)*intervalo)])
        mat_norm[j,x] = trans_j/tot
  return(mat_norm)

#Normalizando por tandas de intervalos de tiempo solapando
def trans_norm_interv_sol(mat, intervalo=5):
  t = len(mat[0])
  mat = contar_trans(mat)
  mat_norm = np.zeros((6,t-intervalo+1))
  for i in range(t-intervalo+1):
    tot = np.sum(mat[:,i:(i+intervalo)])
    #print(tot)
    if tot != 0:
      for j in range(len(mat[:,0])):
        trans_j = sum(mat[j,i:(i+intervalo)])
        #print(trans_j)
        mat_norm[j,i] = trans_j/tot
  return(mat_norm)



#Calculamos la divergencia de Kullback-Leibler
#Especificando la matriz sobre la que actua la funcion segun el normalizado que se desee
def entrop_prod(mat, sigma=1):
  D_t = np.zeros(len(mat[0]))
  for i in range(len(mat[0])):
    D = 0
    for j in range(0, len(mat[:,0]), 2):
      if mat[j,i]!=0 and mat[j+1,i]!=0:
        D = D - (mat[j,i] * np.log(mat[j+1,i]/float(mat[j,i])))
        D = D - (mat[j+1,i] * np.log(mat[j,i]/float(mat[j+1,i])))
    D_t[i] = D
  D_t = gaussian_filter(D_t,sigma=sigma)
  return(D_t)