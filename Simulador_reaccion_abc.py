 # -*- coding: utf-8 -*-

"""
Created on Tue Jul 11 11:48:03 2023

@author: celia
"""

#Simulando un sistema de transicion entre 3 estados: A, B, C
#El sistema puede transicionar de cualquier a cualquier estado con una probabilidad que se define en una lista

#Importamos las librerias necesarias
import matplotlib.pyplot as plt
import numpy as np
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


#Definimos los parametros de simulacion
p_lista = [0.0045, 0.0015, 0.006, 0.0015, 0.004, 0.005]
weights = (60, 10, 30)
n_eventos = 500
t = 1500
seed=8462836

mat = mat_simul(p_lista, weights, n_eventos, t, seed)
mat = count_abc(mat)
print(mat)

#Representamos la evolucion de la cantidad de cada uno de los estados a lo largo del tiempo de simulacion
plt.title('Evolucion de especies en el tiempo')
plt.plot(mat[0], 'r', label='A')
plt.plot(mat[1], 'b', label='B')
plt.plot(mat[2], 'g', label='C')
plt.xlabel('Tiempo')
plt.ylabel('N')
plt.legend(loc='best')
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/evolucion_abc.jpg', 
            dpi=300)
plt.show()

#Representamos las cantidades de unos estados frente a otros
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

plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/proporcion_abc.jpg',
            dpi=300)
plt.show()


