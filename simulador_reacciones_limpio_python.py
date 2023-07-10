# -*- coding: utf-8 -*-
"""Simulador_reacciones_limpio.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1wtZH34lYDnfXfzFf4pW_vKx9Nyllsjfr
"""

import math
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import numpy as np
import random

"""# **Simulador de reaccion A,B**

Para la simulación de una reaccion con un sustrato A que pasa a ser el producto B con una probabilidad de P, siendo la reacción reversible:

1. Primero escribimos una funcion que cree una matriz de apariciones de A y B según la probabilidad especificada tanto para la reaccion directa (p_ab), como para la reversible (p_ba).
Dicha matriz estará compuesta por tantas filas como número de simulaciones que se especifiquen y tantas columnas como tiempo de simulacion.
Esta matriz será reflejo de la cantidad de reactivos o productos hay en cada instante de la reaccion.
"""

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

mat_simul_ab(0.01 , 0.02, 10000, 200)

"""2. Definimos una función que cuente el número de apariciones de A y B en cada instante de tiempo, devolviendo una matriz cuya primera fila es el número de A y la segunda el numero de B."""

def contar_ab(p_ab, p_ba, n_eventos=100, t=100):
  mat = mat_simul_ab(p_ab, p_ba, n_eventos, t)
  count = np.zeros((2, t))
  for j in range(len(mat[0])):
    count[0, j] = np.count_nonzero(mat[:,j] == 'A')
    count[1, j] = np.count_nonzero(mat[:,j] == 'B')
  return(count)

contar_ab(0.01 , 0.02, 100, 20)

"""3. Por último, representamos esta proporción de reactivos y productos en un plot."""

matriz = contar_ab(0.02 , 0.01, 10000, 200)

plt.title('Evolucion de especies en el tiempo')
plt.plot(matriz[0], 'r', label='A')
plt.plot(matriz[1], 'b', label='B')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.legend(loc='best')
plt.show()

"""Para comprobar que la simulación se ha realizado correctamente, representamos la relación Nb/Na, que debería seguir una tendencia que se estabilice en torno al producto de la división de las probabilidades (p_ab/p_ba)."""

def comprobar_simulacion(probabilidad_d, probabilidad_r, n_simulaciones = 100, t = 100):

  na_nb_array = evolucion_reaccion(probabilidad_d, probabilidad_r, n_simulaciones, t)

  y = na_nb_array[1]/na_nb_array[0]
  x = list(range(0, len(na_nb_array[0])))
  k = (probabilidad_d/probabilidad_r)
  y_k = []
  for i in range(0,len(na_nb_array[0])):
    y_k.append(k)

  plt.plot(x,y)
  plt.plot(x,y_k)
  plt.show()

p_ab = 0.02
p_ba = 0.01

mat = contar_ab(p_ab , p_ba, 200, 10000)

y = mat[1]/mat[0]
k = p_ab/p_ba

plt.plot(y, 'b', label = 'Nb/Na')
plt.axhline(y = k, color = 'r', label = 'p_ab/p_ba')
plt.title('Tendencia de la simulacion')
plt.legend(loc = 'best')
plt.xlabel('Tiempo')
plt.show()

"""# **Simulador de reaccion A,B,C**

Realizamos una simulación igual a la anterior pero añadiendo un tercer estado o especie (C).
En este caso, se producirán 3 reacciones reversibles al mismo tiempo, por lo que se tendrán que definir 6 probabilidades.
Dichas probabilidades se definirán en una lista:

p_lista = [p_ab, p_ac, p_ba, p_bc, p_ca, p_cb]

1. Primero definimos una función que cree una matriz que recoja las apariciones de las distintas especies reactivas a lo largo del tiempo de simulacion, igual que en la simulacion de 2 especies.
"""

#p_lista = [p_ab, p_ac, p_ba, p_bc, p_ca, p_cb]

def mat_simul(p_lista, n_eventos=100, t=100):
  mat = np.full((n_eventos, t), 'A')
  for i in range(n_eventos):
    for j in range(1,t):
      x = random.random()

      if mat[i,j-1]=='A':
        if x > (p_lista[0]+p_lista[1]):
          mat[i,j] = 'A'
        elif x <= p_lista[0]:
          mat[i,j] = 'B'
        else: #(p_lista[0] < x >= (p_lista[0]+p_lista[1])):
          mat[i,j] = 'C'

      elif mat[i,j-1]=='B':
        if x > (p_lista[2]+p_lista[3]):
          mat[i,j] = 'B'
        elif x <= p_lista[2]:
          mat[i,j] = 'A'
        else: #p_lista[2] < x or x >= (p_lista[2]+p_lista[3]):
          mat[i,j] = 'C'

      elif mat[i,j-1]=='C':
        if x > (p_lista[4]+p_lista[5]):
          mat[i,j] = 'C'
        elif x <= p_lista[4]:
          mat[i,j] = 'A'
        else: #p_lista[4] < x >= (p_lista[4]+p_lista[5]):
          mat[i,j] = 'B'

  return mat

p_lista = [0.1, 0.2, 0.2, 0.1, 0.1, 0.05]

mat_simul(p_lista, 200, 1000)

"""2. Definimos una función que cuente las apariciones de las distintas especies (N), igual que con 2. Las filas 1, 2 y 3 se corresponderán con A, B y C respectivamente."""

def count_abc(p_lista, n_eventos=100, t=100):
  mat = mat_simul(p_lista, n_eventos, t)
  count = np.zeros((3, t))
  for j in range(len(mat[0])):
    count[0, j] = np.count_nonzero(mat[:,j] == 'A')
    count[1, j] = np.count_nonzero(mat[:,j] == 'B')
    count[2, j] = np.count_nonzero(mat[:,j] == 'C')
  return(count)

p_lista = [0.1, 0.2, 0.2, 0.1, 0.1, 0.05]
count_abc(p_lista, 20, 10)

"""3. Representamos la evolución de las especies en el tiempo en un plot."""

p_lista = [0.0045, 0.0015, 0.006, 0.0015, 0.004, 0.005]
mat = count_abc(p_lista, 500, 1500)

plt.title('Evolucion de especies en el tiempo')
plt.plot(mat[0], 'r', label='A')
plt.plot(mat[1], 'b', label='B')
plt.plot(mat[2], 'g', label='C')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.legend(loc='best')
plt.show()
plt.show()

"""4. Para visualizarlo mejor representamos las apariciones de unas especies frente a otras:"""

p_lista = [0.0045, 0.0015, 0.006, 0.0015, 0.004, 0.005]
mat = count_abc(p_lista, 500, 1500)

fig, ax = plt.subplots(2,2, sharex = True, sharey = True)
fig.suptitle('Proporcion de las especies')

ax[0,0].plot(mat[0], mat[1]) #A frente B
ax[0,0].set_xlabel('Na')
ax[0,0].set_ylabel('Nb')

ax[1,0].plot(mat[2], mat[1]) #C frente B
ax[1,0].set_xlabel('Nc')
ax[1,0].set_ylabel('Nb')

ax[0,1].plot(mat[0], mat[2]) #A frente C
ax[0,1].set_xlabel('Na')
ax[0,1].set_ylabel('Nc')

plt.show()

"""## **Modificando las condiciones iniciales**

Con el método anterior, se simula la reacción con A como único reactivo inicial. Para comenzarcon distintas concentraciones de las especies en el estado inicial añadimos una función que decida de manera aleatoria el reactivo inicial (A, B o C) con distintos pesos o probabilidades de salir.
"""

#p_lista = [p_ab, p_ac, p_ba, p_bc, p_ca, p_cb]
#weights=(A, B, C) -> probabilidad de que salga cada una de las especies (=proporcion inicial de los reactivos)

def mat_simul(p_lista, weights=(1/3, 1/3, 1/3), n_eventos=100, t=100):
  mat = np.full((n_eventos, t), 'A')
  estados = ['A', 'B', 'C']

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

def count_abc(p_lista, weights, n_eventos=100, t=100):
  mat = mat_simul(p_lista, weights, n_eventos, t)
  count = np.zeros((3, t))
  for j in range(len(mat[0])):
    count[0, j] = np.count_nonzero(mat[:,j] == 'A')
    count[1, j] = np.count_nonzero(mat[:,j] == 'B')
    count[2, j] = np.count_nonzero(mat[:,j] == 'C')
  return(count)

p_lista = [0.0045, 0.0015, 0.006, 0.0015, 0.004, 0.005]
weights = (10, 20, 70)
mat = count_abc(p_lista, weights, 500, 1500)

plt.title('Evolucion de especies en el tiempo')
plt.plot(mat[0], 'r', label='A')
plt.plot(mat[1], 'b', label='B')
plt.plot(mat[2], 'g', label='C')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.legend(loc='best')

fig, ax = plt.subplots(2,2, sharex = True, sharey = True)
fig.suptitle('Proporcion de las especies')

ax[0,0].plot(mat[0], mat[1]) #A frente B
ax[0,0].set_xlabel('Na')
ax[0,0].set_ylabel('Nb')

ax[1,0].plot(mat[2], mat[1]) #C frente B
ax[1,0].set_xlabel('Nc')
ax[1,0].set_ylabel('Nb')

ax[0,1].plot(mat[0], mat[2]) #A frente C
ax[0,1].set_xlabel('Na')
ax[0,1].set_ylabel('Nc')

plt.show()

"""## **Mismo número semilla**

Para poder reproducir las simulaciones con unas mismas condiciones iniciales y evolución en el tiempo añadimos un parámetro nuevo, el número semilla, para que los números aleatorios generados sean iguales si se repite la simulación.
"""

#p_lista = [p_ab, p_ac, p_ba, p_bc, p_ca, p_cb]
#weights=(A, B, C) -> probabilidad de que salga cada una de las especies (=proporcion inicial de los reactivos)

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

#para comprobar si todos salen igual con la misma semilla: si que salen todos igual
p_lista = [0.01, 0.02, 0.02, 0.01, 0.01, 0.05]
weights = (60, 10, 30)

print('Iguales:')
for x in range(3):
  print(mat_simul(p_lista, weights, n_eventos=100, t=100, seed=8462836), '\n')

print('Diferentes')
for x in range(3):
  print(mat_simul(p_lista, weights, n_eventos=100, t=100, seed=25), '\n')

"""Para poder reproducir los experimentos, también extraemos la matriz de las funciones de conteo y posteriores para poder usar la misma matriz para calcular distintos parámetros."""

def count_abc(mat):
  t = len(mat[0])
  count = np.zeros((3, t))
  for j in range(len(mat[0])):
    count[0, j] = np.count_nonzero(mat[:,j] == 'A')
    count[1, j] = np.count_nonzero(mat[:,j] == 'B')
    count[2, j] = np.count_nonzero(mat[:,j] == 'C')
  return(count)

"""## **Cálculo de entropía Shannon**

La entropía presente en un sistema en el que se da una reacción se puede calcular a partir de la probabilidad

Añadimos además el parámetro sigma, con el que vamos a ser capaces de realizar un suavizado más o menos brusco aplicando un filtro gaussiano de n puntos.
"""

def entrop_simul(mat, n_eventos, sigma=3):
  mat = count_abc(mat)
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

"""Representamos la evolución de la entropía de la reacción en el tiempo junto con la evolución de las especies en el tiempo."""

p_lista = [0.001, 0.02, 0.02, 0.01, 0.01, 0.05]
weights = (90, 10, 30)

mat = mat_simul(p_lista, weights, n_eventos=200, t=1000, seed=8462836)
mat_N = count_abc(mat)
mat_S = entrop_simul(mat, n_eventos=200, sigma=11)

plt.title('Evolucion de especies en el tiempo')   #evolucion de poblaciones
plt.plot(mat_N[0], 'r', label='Na')
plt.plot(mat_N[1], 'b', label='Nb')
plt.plot(mat_N[2], 'g', label='Nc')
plt.legend(loc='best')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.show()

plt.plot(mat_S)  #evolucion de entropia
plt.title('Evolucion de entropia')
plt.xlabel('Tiempo')
plt.ylabel('S')
plt.show()

"""## **Variación de entropía**

Representamos también la derivada de la entropía de la reacción para observar la variación de la entropía del entorno durante el transcurso de la reacción.

La derivada se calculará a partir de la matriz con los valores de entropía del apartado anterior.
"""

def var_entrop(mat):
  delta_St = []
  for x in range(1, len(mat)):
    delta_S = mat[x]-mat[x-1]
    delta_St.append(delta_S)
  return(delta_St)

p_lista = [0.001, 0.02, 0.02, 0.01, 0.01, 0.05]
weights = (90, 10, 30)

mat = mat_simul(p_lista, weights, n_eventos=200, t=1000, seed=8462836)
mat_S = entrop_simul(mat, n_eventos=200, sigma=5)
mat_St = var_entrop(mat_S)

fig, ax = plt.subplots(2,1, sharex = True)

ax[0].plot(mat_S, label = 'Evolucion de entropia')
ax[0].legend()
ax[0].set_ylabel('S')

ax[1].plot(mat_St, 'r', label = 'Variacion de entropia')
ax[1].legend()
ax[1].set_ylabel('dS/dt')
ax[1].set_xlabel('Tiempo')

plt.show()

"""## **Mostrar la media de un número n de simulaciones**"""

def graf_media_especies(p_lista, weights=(1/3, 1/3, 1/3), n_simulaciones=5, n_eventos=100, t=100, sigma=3):
  count_media = np.zeros((3, t))
  count_media_2 = np.zeros((3, t))
  fig, ax = plt.subplots()
  for x in range(n_simulaciones):
    mat = mat_simul(p_lista, weights, n_eventos, t, seed=8462836)
    count = count_abc(mat)
    for i in range(len(count[:,0])):
      count[i] = gaussian_filter(count[i],sigma=sigma)
    plt.plot(count[0],'lightsalmon', linewidth=1)
    plt.plot(count[1], 'lightblue', linewidth=1)
    plt.plot(count[2], 'lightgreen', linewidth=1)

    count_media = count_media + count
    count_media_2 = count_media_2 + count**2
  count_media = count_media/n_simulaciones
  count_media_2 = count_media_2/n_simulaciones

  plt.plot(count_media[0], 'r', label='media_Na')
  plt.plot(count_media[1], 'b', label='media_Nb')
  plt.plot(count_media[2], 'g', label='media_Nc')

  stdr_dev = np.sqrt(count_media_2 - count_media**2)
  plt.plot(count_media[0] + stdr_dev[0], '--r', linewidth=1)
  plt.plot(count_media[1] + stdr_dev[1], '--b', linewidth=1)
  plt.plot(count_media[2] + stdr_dev[2], '--g', linewidth=1)
  plt.plot(count_media[0] - stdr_dev[0], '--r', linewidth=1)
  plt.plot(count_media[1] - stdr_dev[1], '--b', linewidth=1)
  plt.plot(count_media[2] - stdr_dev[2], '--g', linewidth=1)

  ax.legend(loc='best')
  plt.title('Evolucion de especies reactivas')
  plt.xlabel('Tiempo')
  plt.ylabel('Cantidad')
  plt.show()

p_lista = [0.001, 0.002, 0.00, 0.001, 0.01, 0.005]
weights = (60, 10, 30)

graf_media_especies(p_lista, weights, n_simulaciones=20, n_eventos=100, t=100, sigma=1)

def graf_entrop_media(p_lista, weights=(1/3, 1/3, 1/3), n_simulaciones=5, n_eventos=100, t=100, sigma=3):
  fig, ax = plt.subplots()
  S_media = np.zeros((t))
  S_media_2 = np.zeros((t))

  for x in range(n_simulaciones):
    mat = mat_simul(p_lista, weights, n_eventos, t)
    mat = entrop_simul(mat, n_eventos, sigma)
    plt.plot(mat, 'lightblue', linewidth=1)
    S_media = S_media + mat
    S_media_2 = S_media_2 + mat**2

  S_media = S_media/n_simulaciones
  plt.plot(S_media, 'b')

  S_media_2 = S_media_2/n_simulaciones
  stdr_dev = np.sqrt(S_media_2 - S_media**2)
  plt.plot(S_media + stdr_dev, '--b', linewidth=1)
  plt.plot(S_media - stdr_dev, '--b', linewidth=1)

  plt.title('Evolucion de entropia')
  plt.xlabel('Tiempo')
  plt.ylabel('S')
  plt.show()

p_lista = [0.001, 0.002, 0.02, 0.001, 0.01, 0.005]
weights = (60, 10, 30)

graf_entrop_media(p_lista, weights, n_simulaciones=20, n_eventos=100, t=100, sigma=1)

p_lista = [0.001, 0.002, 0.002, 0.001, 0.001, 0.005]
weights = (60, 10, 30)

graf_media_especies(p_lista, weights, n_simulaciones=20, n_eventos=100, t=500, sigma=1)
graf_entrop_media(p_lista, weights, n_simulaciones=20, n_eventos=100, t=500, sigma=1)

"""## **Conteo de transiciones**"""

#mat_trans[:, 0] = (Nab, Nba, Nac, Nca, Nbc, Ncb)
#posiciones      =   0    1    2    3    4    5

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

"""Normalizamos el número de transiciones"""

def trans_norm(mat):
  mat = contar_trans(mat)
  #print(mat)

  for i in range(len(mat[0])):
    tot = np.sum(mat[:,i])
    for j in range(len(mat[:,0])):
      if mat[j,i] != 0:
        mat[j,i] = mat[j,i]/tot

  return(mat)

"""## **Divergencia de Kullback-Leibler**

Calculamos la entropía realtiva o divergencia de la información, que se traduce en la información como la energía que debe invertir el sistema para lograr el estado que observa en el transcurso del tiempo.
"""

#mat_trans[:, 0] = (Nab, Nba, Nac, Nca, Nbc, Ncb)
#posiciones      =   0    1    2    3    4    5

def entrop_prod(mat, sigma=1):
  mat = trans_norm(mat)
  #print(np.round(mat, decimals=2))
  D_t = np.zeros(len(mat[0]))

  for i in range(len(mat[0])):
    D = 0
    for j in range(0, len(mat[:,0]), 2):
      if mat[j,i]!=0 and mat[j+1,i]!=0:
        D = D - (mat[j,i] * np.log(mat[j+1,i]/mat[j,i]))
        D = D - (mat[j+1,i] * np.log(mat[j,i]/mat[j+1,i]))
    D_t[i] = D

  D_t = gaussian_filter(D_t,sigma=sigma)

  return(D_t)

"""Representamos la divergencia de Kullback-Leibler"""

#p_lista = [0.1, 0.2, 0.2, 0.1, 0.1, 0.05]
#weights = (20, 20, 30)

p_lista = [0.001, 0.002, 0.002, 0.001, 0.001, 0.005]
weights = (90, 10, 30)

mat = mat_simul(p_lista, weights, n_eventos=1000, t=1000, seed=8462836)
mat = entrop_prod(mat, sigma=10)

plt.plot(mat, 'b')
plt.title('Divergencia de Kullback-Leibler')
plt.xlabel('Tiempo')
plt.ylabel('dS/dt')
plt.show()

"""## **Entropia normalizando transiciones por tandas de intervalos n**

Normalizando en cada instante de tiempo da lugar a una matriz con muchos ceros. Para evitar esto normalizamos en tandas de n intervalos de tiempo. En vez de usar sigma se hace la media con el número de transiciones en n instantes de tiempo.

### Sin solapar
"""

#normalizando por tandas sin solapar

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

#mat_trans[:, 0] = (Nab, Nba, Nac, Nca, Nbc, Ncb)
#posiciones      =   0    1    2    3    4    5

def entrop_prod_inter(mat, sigma=1, intervalo=5):
  mat = trans_norm_interv(mat, intervalo)
  #print(np.round(mat, decimals=2))
  D_t = np.zeros(len(mat[0]))

  for i in range(len(mat[0])):
    D = 0
    for j in range(0, len(mat[:,0]), 2):
      if mat[j,i]!=0 and mat[j+1,i]!=0:
        D = D - (mat[j,i] * np.log(mat[j+1,i]/mat[j,i]))
        D = D - (mat[j+1,i] * np.log(mat[j,i]/mat[j+1,i]))
    D_t[i] = D

  D_t = gaussian_filter(D_t,sigma=sigma)

  return(D_t)

p_lista = [0.01, 0.02, 0.02, 0.01, 0.01, 0.05]
weights = (60, 10, 30)

matriz = mat_simul(p_lista, weights, n_eventos=1000, t=500)
mat = entrop_prod_inter(matriz, sigma=1, intervalo=2)
plt.plot(mat, 'b')
plt.title('Divergencia de Kullback-Leibler')
plt.xlabel('Tiempo')
plt.ylabel('dS/dt')
plt.show()

"""### Solapando"""

#normalizando por intervalos solapantes

def trans_norm_interv_sol(mat, intervalo=5):
  t = len(mat[0])
  mat = contar_trans(mat)
  #print(mat)
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

#mat_trans[:, 0] = (Nab, Nba, Nac, Nca, Nbc, Ncb)
#posiciones      =   0    1    2    3    4    5

def entrop_prod_inter_sol(mat, intervalo=5, sigma=1):
  mat = trans_norm_interv_sol(mat, intervalo)
  #print(np.round(mat, decimals=2))
  D_t = np.zeros(len(mat[0]))

  for i in range(len(mat[0])):
    D = 0
    for j in range(0, len(mat[:,0]), 2):
      if mat[j,i]!=0 and mat[j+1,i]!=0:
        D = D - (mat[j,i] * np.log(mat[j+1,i]/mat[j,i]))
        D = D - (mat[j+1,i] * np.log(mat[j,i]/mat[j+1,i]))
    D_t[i] = D

  D_t = gaussian_filter(D_t,sigma=sigma)

  return(D_t)

p_lista = [0.01, 0.02, 0.02, 0.01, 0.01, 0.05]
weights = (60, 10, 30)

matriz = mat_simul(p_lista, weights, n_eventos=1000, t=500, seed=8462836)
mat = entrop_prod_inter_sol(matriz, intervalo=2, sigma=1)

plt.plot(mat, 'b')
plt.title('Divergencia de Kullback-Leibler')
plt.xlabel('Tiempo')
plt.ylabel('dS/dt')
plt.show()

"""### Representación de la media"""

def graf_entrop_prod(p_lista, weights=(1/3, 1/3, 1/3), n_eventos=100, t=100, intervalo=5 ,sigma=3, n_simul=5):
  fig, ax = plt.subplots()
  matriz = mat_simul(p_lista, weights, n_eventos, t, seed=8462836)
  mat = entrop_prod_inter(matriz, intervalo, sigma)
  media = np.zeros((len(mat)))
  media_2 = np.zeros((len(mat)))
  #eje_x = range(0, t, t//intervalo-1)

  for x in range(n_simul):
    matriz = mat_simul(p_lista, weights, n_eventos, t, seed=8462836)
    mat = entrop_prod_inter(matriz, intervalo, sigma)
    plt.plot(mat, 'lightblue', linewidth=1)
    media = media + mat
    media_2 = media_2 + mat**2

  media = media/n_simul
  plt.plot(media, 'b')

  media_2 = media_2/n_simul
  stdr_dev = np.sqrt(media_2 - media**2)
  plt.plot(media + stdr_dev, '--b', linewidth=1)
  plt.plot(media - stdr_dev, '--b', linewidth=1)

  fig.suptitle('Divergencia de Kullback-Leibler')
  plt.title('Sin solapamiento', fontsize=10)
  plt.xlabel('Tiempo')
  plt.ylabel('dS/dt')
  plt.show()

def graf_entrop_prod_sol(p_lista, weights=(1/3, 1/3, 1/3), n_eventos=100, t=100, intervalo=5 ,sigma=3, n_simul=5):
  fig, ax = plt.subplots()
  media = np.zeros((t-intervalo+1))
  media_2 = np.zeros((t-intervalo+1))

  for x in range(n_simul):
    matriz = mat_simul(p_lista, weights, n_eventos, t, seed=8462836)
    mat = entrop_prod_inter_sol(matriz, intervalo, sigma)
    plt.plot(mat, 'lightblue', linewidth=1)
    media = media + mat
    media_2 = media_2 + mat**2

  media = media/n_simul
  plt.plot(media, 'b')

  media_2 = media_2/n_simul
  stdr_dev = np.sqrt(media_2 - media**2)
  plt.plot(media + stdr_dev, '--b', linewidth=1)
  plt.plot(media - stdr_dev, '--b', linewidth=1)

  fig.suptitle('Divergencia de Kullback-Leibler')
  plt.title('Con solapamiento', fontsize=10)
  plt.xlabel('Tiempo')
  plt.ylabel('dS/dt')
  plt.show()

p_lista = [0.01, 0.02, 0.02, 0.01, 0.01, 0.05]
weights = (60, 10, 30)

graf_entrop_prod(p_lista, weights, n_eventos=1000, t=200, intervalo=5 ,sigma=3, n_simul=20)
graf_entrop_prod_sol(p_lista, weights, n_eventos=1000, t=200, intervalo=5 ,sigma=3, n_simul=20)

"""# **Entropia por tiempos de residencia en un estado**"""

import pandas as pd

"""## **Para una reacción con dos estados y reversible**

Utilizamos como matriz inicial la matriz de A y B a lo largo del tiempo y tantas veces como n eventos.
"""

def mat_simul_ab(p_ab, p_ba, weights, n_eventos=100, t=100):
  mat = np.full((n_eventos, t), 'A')
  mat[:, 0] = random.choices(['A','B'], weights=weights, k=n_eventos)
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

weights=(20,10)
mat_simul_ab(p_ab=0.1, p_ba=0.05, weights=weights, n_eventos=15, t=10)

"""Diseñamos una función que calcule los tiempos de residencia en cada estado y cada evento.

Juntamos cada fila de la matriz en un solo elemento string y sustituimos las transiciones AB o BA por A#B o B#A. Luego utilizamos la función split para separar en los # e incluimos los tiempos de residencia en un diccionario.

Una vez recogido todo en el diccionario representamos un histograma y calculamos la varianza y media de los tiempos de residencia.


"""

def tiempos_residencia(mat):
  t_res = pd.DataFrame(columns=['estado', 'tiempo_residencia', 'frecuencia'])  #tiempos de residencia

  for j in range(len(mat[:,0])):
    event = []

    for i in range(1,len(mat[j])):
      if mat[j,i-1]==mat[j,i]:
        event.append(mat[j,i-1])
      elif mat[j,i-1]!=mat[j,i]:
        event += mat[j, i-1] + '#'
    event.append(mat[j, len(mat[0])-1])

    event = ''.join(event)
    event = event.split(sep='#')

    for x in range(len(event)):
      est = event[x][0]
      rep = len(event[x])
      prev = t_res.loc[(t_res.estado == est)&(t_res.tiempo_residencia == rep)]
      if len(prev) != 0:
        f = int(prev.frecuencia) + 1
        t_res.loc[(t_res.estado == est)&(t_res.tiempo_residencia == rep)] = [est, rep, f]
      else:
        f = 1
        t_res_2 = pd.DataFrame({'estado':est, 'tiempo_residencia':[rep], 'frecuencia':[f]})
        t_res = pd.concat([t_res, t_res_2], ignore_index=True)

  t_res = t_res.sort_values(['estado', 'tiempo_residencia'])
  return(t_res)

weights=(20,10)
mat = mat_simul_ab(p_ab=0.1, p_ba=0.05, weights=weights, n_eventos=100, t=200)
tiempos_residencia(mat)

"""### **Gráfico de frecuencias**

Representamos las frecuencias en un histograma, que debería seguir una distribución normal. Hacemos dos graficos separados para A y B.

Primero creamos una función que extrae los datos del diccionario para cada una de las especies de reacción en forma de matriz para su representación en las gráficas.
"""

def mat_frec(dic, estado):   #El estado debe aparecer entre comillas y mayuscula
  dic_n = dic.loc[dic.estado == estado]
  f_n = np.array([dic_n.frecuencia, dic_n.tiempo_residencia])
  return f_n

fig, ax = plt.subplots(figsize=(25,10))
weights=(20,10)
mat = mat_simul_ab(p_ab=0.2, p_ba=0.1, weights=weights, n_eventos=200, t=200)
dic = tiempos_residencia(mat)

f_A = mat_frec(dic, 'A')
f_B = mat_frec(dic, 'B')
print(f_A)

plt.bar(f_A[1],f_A[0], color='b', alpha=0.5)
plt.bar(f_B[1],f_B[0], color='r', alpha=0.5)
#(counts, bins, patches) = plt.hist(f_B)
plt.title('Grafico de frecuencias', fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('tiempo de residencia', fontsize=15)
plt.ylabel('N', fontsize=15)
plt.show()

"""Calculamos las entropías para cada una de los estados A y B, siendo esta la varianza de las distribución entre la media al cuadrado."""

tAmean = np.sum(f_A[0]*f_A[1])/np.sum(f_A[0])
t2Amean = np.sum(f_A[0]*f_A[1]**2)/np.sum(f_A[0])
vartA = t2Amean-tAmean**2

print(vartA/(tAmean**2))

tBmean = np.sum(f_B[0]*f_B[1])/np.sum(f_B[0])
t2Bmean = np.sum(f_B[0]*f_B[1]**2)/np.sum(f_B[0])
vartB = t2Bmean-tBmean**2

print(vartB/(tBmean**2))

"""## **Cálculo de S**

Calcular S en tramos de tiempo.

Añadimos un bucle for con k tramos de tiempo y distinguimos en el diccionario con una columna más en qué tramo de tiempo está ese tiempo de residencia. Para hacer esto la matriz se corta en los n tramos de tiempo de forma que si había un tiempo de residencia muy largo que abarca dos tramos de tiempo es separado.
"""

def tiempos_residencia_tramos(mat, k):
  t_res = pd.DataFrame(columns=['estado', 'tramo_tiempo', 'tiempo_residencia', 'frecuencia'])  #tiempos de residencia
  n = len(mat[0])//k
  #print(k)
  #print(len(mat[0])//n)

  for j in range(len(mat[:,0])):

    for k in range(len(mat[0])//n):
      event = []
      tramo = 0

      for i in range(1+k*n, (k+1)*n):
        if mat[j,i-1]==mat[j,i]:
          event.append(mat[j,i-1])
        elif mat[j,i-1]!=mat[j,i]:
          event += mat[j, i-1] + '#'
      event.append(mat[j, len(mat[0])-1])

      event = ''.join(event)
      event = event.split(sep='#')

      for x in range(len(event)):
        est = event[x][0]
        rep = len(event[x])
        tramo = k
        prev = t_res.loc[(t_res.estado == est)&(t_res.tiempo_residencia == rep)&(t_res.tramo_tiempo == k)]
        if len(prev) != 0:
          f = int(prev.frecuencia) + 1
          t_res.loc[(t_res.estado == est)&(t_res.tiempo_residencia == rep)&(t_res.tramo_tiempo == k)] = [est, tramo, rep, f]
        else:
          f = 1
          t_res_2 = pd.DataFrame({'estado':est, 'tramo_tiempo':tramo, 'tiempo_residencia':[rep], 'frecuencia':[f]})
          t_res = pd.concat([t_res, t_res_2], ignore_index=True)

  t_res = t_res.sort_values(['estado', 'tramo_tiempo', 'tiempo_residencia'])
  return(t_res)

"""Modificamos la función anterior para que también distinga en el tramo de tiempo."""

def mat_frec_tramos(dic, estado, tramo_tiempo):   #El estado debe aparecer entre comillas y mayuscula
  dic_n = dic.loc[(dic.estado == estado)&(dic.tramo_tiempo == tramo_tiempo)]
  f_n = np.array([dic_n.frecuencia, dic_n.tiempo_residencia])
  return f_n

"""Creamos una función que calcule S con el array de frecuencias y tiempos de residencia.

A partir de f_n se calcula la varianza y la media y se calcula S como var/media^2
"""

def calculo_S(f_n):
  ti_mean = np.sum(f_n[0]*f_n[1])/np.sum(f_n[0])
  t2i_mean = np.sum(f_n[0]*f_n[1]**2)/np.sum(f_n[0])
  var_ti = t2i_mean-ti_mean**2

  return(var_ti/(ti_mean**2))

"""Creamos un bucle que vaya calculando las entropias de cada tramo de residencia y de cada especie y luego los represente."""

#volvemos a ejecutar la función de conteo de cada estado en la matriz inicial.

def count_abc(mat):
  t = len(mat[0])
  count = np.zeros((3, t))
  for j in range(len(mat[0])):
    count[0, j] = np.count_nonzero(mat[:,j] == 'A')
    count[1, j] = np.count_nonzero(mat[:,j] == 'B')
    count[2, j] = np.count_nonzero(mat[:,j] == 'C')
  return(count)

#p_lista = [p_ab, p_ac, p_ba, p_bc, p_ca, p_cb]
p_lista = [0.001, 0.002, 0.002, 0.001, 0.001, 0.005]

weights = (60, 10, 10)
n = 5

mat = mat_simul(p_lista, weights, n_eventos=50, t=1000, seed=8462836)
dic = tiempos_residencia_tramos(mat, n)

mat_N = count_abc(mat)
#mat_S = entrop_simul(mat, n_eventos=200, sigma=11)

plt.title('Evolucion de especies en el tiempo')   #evolucion de poblaciones
plt.plot(mat_N[0], 'r', label='Na')
plt.plot(mat_N[1], 'b', label='Nb')
plt.plot(mat_N[2], 'g', label='Nc')
plt.legend(loc='best')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.show()

#Entropia calculada en cada tramo y representada en conjunto
fig, ax = plt.subplots()
estados = ['A', 'B', 'C']
colores = ['r', 'b', 'g']
for i in range(len(estados)):
  S_list = []
  for j in range(n):
    f_n = mat_frec_tramos(dic, estados[i], j)
    S = calculo_S(f_n)
    S_list.append(S)
  plt.plot(S_list, color=colores[i], label=estados[i])
plt.title('Entropia por tramos')
plt.xlabel('Tramos tiempo')
plt.ylabel('S')
ax.legend(loc='best')
plt.show()

#Representación de las distribuciones de los tiempos de residencia en cada tramo
fig, ax = plt.subplots(n, 1, sharex=True)
for j in range(n):
  for i in range(len(estados)):
    f_n = mat_frec_tramos(dic, estados[i], j)
    ax[j].bar(f_n[1],f_n[0], color=colores[i], label=estados[i], alpha=0.5)
    ax[j].set_ylabel('N'+str(j))

fig.suptitle('Grafico de frecuencias')
ax[n-1].set_xlabel('tiempo de residencia')
plt.show()



#p_lista = [p_ab, p_ac, p_ba, p_bc, p_ca, p_cb]
p_lista = [0.001, 0.002, 0.002, 0.001, 0.001, 0.005]

weights = (60, 10, 10)
n = 5

mat = mat_simul(p_lista, weights, n_eventos=50, t=1000, seed=8462836)
mat = mat[:, 0:200]
print(mat)
dic = tiempos_residencia(mat)

mat_N = count_abc(mat)
#mat_S = entrop_simul(mat, n_eventos=200, sigma=11)

plt.title('Evolucion de especies en el tiempo')   #evolucion de poblaciones
plt.plot(mat_N[0], 'r', label='Na')
plt.plot(mat_N[1], 'b', label='Nb')
plt.plot(mat_N[2], 'g', label='Nc')
plt.legend(loc='best')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.show()

#Entropia calculada en cada tramo y representada en conjunto
fig, ax = plt.subplots()
estados = ['A', 'B', 'C']
colores = ['r', 'b', 'g']
for i in range(len(estados)):
  f_n = mat_frec(dic, estados[i])
  S = calculo_S(f_n)
  print('S de', estados[i], ': ', S)

#Representación de las distribuciones de los tiempos de residencia en cada tramo
for i in range(len(estados)):
  f_n = mat_frec(dic, estados[i])
  plt.bar(f_n[1],f_n[0], color=colores[i], label=estados[i], alpha=0.5)

plt.ylabel('N')
plt.title('Grafico de frecuencias')
plt.xlabel('tiempo de residencia')
plt.show()

for i in mat:
  #print(i)
  print(np.unique(i))
  #print(len(i))
  break

from collections import Counter

dictionary_eventos = Counter(i)
print(dictionary_eventos)
print(len(dictionary_eventos.keys()))

