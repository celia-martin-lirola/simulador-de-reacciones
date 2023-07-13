 # -*- coding: utf-8 -*-

"""
Created on Tue Jul 11 11:48:03 2023

@author: celia
"""
import f_simul
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

#Definimos los parametros de simulacion
p_lista = [0.0045, 0.0015, 0.006, 0.0015, 0.004, 0.005]
weights = (60, 10, 30)
n_eventos = 500
t = 1500
seed=8462836
sigma=3
n_simulaciones = 20

mat = f_simul.mat_simul(p_lista, weights, n_eventos, t, seed)
print(mat)
mat = f_simul.count_abc(mat)
mat_S = f_simul.entrop_simul(mat, n_eventos, sigma)
mat_St = f_simul.var_entrop(mat_S)

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

#Representacion de la entropia de Shannon
plt.plot(mat_S)  #entropia a lo largo del tiempo
plt.title('Evolucion de entropia')
plt.xlabel('Tiempo')
plt.ylabel('S')
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/entropia_shannon.jpg',
            dpi=300)
plt.show()

#Representacion de la variacion de la entropia de Shannon
fig, ax = plt.subplots(2,1, sharex = True)
ax[0].plot(mat_S, label = 'Evolucion de entropia')
ax[0].legend()
ax[0].set_ylabel('S')
ax[1].plot(mat_St, 'r', label = 'Variacion de entropia')
ax[1].legend()
ax[1].set_ylabel('dS/dt')
ax[1].set_xlabel('Tiempo')
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/entropia_shannon_var.jpg',
            dpi=300)
plt.show()


#Representacion de varias simulaciones a la vez junto con su media y desviacion estandar
count_media = np.zeros((3, t))
count_media_2 = np.zeros((3, t))
fig, ax = plt.subplots()
for x in range(n_simulaciones):
  mat = f_simul.mat_simul(p_lista, weights, n_eventos, t, seed=8462836)
  count = f_simul.count_abc(mat)
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
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/media_evolucion_abc.jpg',
            dpi=300)
plt.show()
