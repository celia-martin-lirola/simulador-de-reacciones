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
#p_lista = [p_ab, p_ac, p_ba, p_bc, p_ca, p_cb]
p_lista = [0.0045, 0.0015, 0.006, 0.0015, 0.004, 0.005]
#weights=(A, B, C)
weights = (60, 10, 30)
n_eventos = 500
t = 1500
seed=8462836
sigma=3

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


#Repetir la simulacion n veces para mostrar promedio
n_simulaciones = 20

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


#Representacion de la media de entopias de Shannon calculadas
#Junto con su media y desviacion estandar
fig, ax = plt.subplots()
S_media = np.zeros((t))
S_media_2 = np.zeros((t))

for x in range(n_simulaciones):
  mat = f_simul.mat_simul(p_lista, weights, n_eventos, t)
  mat = f_simul.count_abc(mat)
  mat = f_simul.entrop_simul(mat, n_eventos, sigma)
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
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/media_entropia_shannon.jpg',
            dpi=300)
plt.show()


#Representacion de la divergencia de Kullback-Leibler
intervalo = 5
sigma = 15
#Se representa el conjunto y la media de n_simulaciones
#Sin solapamiento
fig, ax = plt.subplots()
mat_norm = f_simul.trans_norm_interv(mat)
mat_Skl = f_simul.entrop_prod(mat_norm, sigma)
media = np.zeros((len(mat_Skl)))
media_2 = np.zeros((len(mat_Skl)))

for x in range(n_simulaciones):
  mat_n = f_simul.mat_simul(p_lista, weights, n_eventos, t, seed=8462836)
  mat_norm = f_simul.trans_norm_interv(mat_n)
  mat_Skl = f_simul.entrop_prod(mat_norm, sigma)
  plt.plot(mat_Skl, 'lightblue', linewidth=1)
  media = media + mat_Skl
  media_2 = media_2 + mat_Skl**2

media = media/n_simulaciones
plt.plot(media, 'b')

media_2 = media_2/n_simulaciones
stdr_dev = np.sqrt(media_2 - media**2)
plt.plot(media + stdr_dev, '--b', linewidth=1)
plt.plot(media - stdr_dev, '--b', linewidth=1)

fig.suptitle('Divergencia de Kullback-Leibler')
plt.title('Sin solapamiento', fontsize=10)
plt.xlabel('Tiempo')
plt.ylabel('dS/dt')
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/divergencia_KL.jpg',
            dpi=300)
plt.show()

#Se representa el conjunto y la media de n_simulaciones
#Con solapamiento
fig, ax = plt.subplots()
mat_norm = f_simul.trans_norm_interv_sol(mat)
mat_Skl = f_simul.entrop_prod(mat_norm, sigma)
media = np.zeros((len(mat_Skl)))
media_2 = np.zeros((len(mat_Skl)))

for x in range(n_simulaciones):
  mat_n = f_simul.mat_simul(p_lista, weights, n_eventos, t, seed=8462836)
  mat_norm = f_simul.trans_norm_interv_sol(mat_n)
  mat_Skl = f_simul.entrop_prod(mat_norm, sigma)
  plt.plot(mat_Skl, 'lightblue', linewidth=1)
  media = media + mat_Skl
  media_2 = media_2 + mat_Skl**2

media = media/n_simulaciones
plt.plot(media, 'b')

media_2 = media_2/n_simulaciones
stdr_dev = np.sqrt(media_2 - media**2)
plt.plot(media + stdr_dev, '--b', linewidth=1)
plt.plot(media - stdr_dev, '--b', linewidth=1)

fig.suptitle('Divergencia de Kullback-Leibler')
plt.title('Sin solapamiento', fontsize=10)
plt.xlabel('Tiempo')
plt.ylabel('dS/dt')
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/divergencia_KL_solapamiento.jpg',
            dpi=300)
plt.show()


#Representamos la distribucion de los tiempos de residencia de cada una de las especies
#Definimos n como el numero de tramos en los que dividimos la simulacion
n=10
dic = f_simul.tiempos_residencia_tramos(mat, n)
fig, ax = plt.subplots()
estados = ['A', 'B', 'C']
colores = ['r', 'b', 'g']
for i in range(len(estados)):
  S_list = []
  for j in range(n):
    f_n = f_simul.mat_frec_tramos(dic, estados[i], j)
    S = f_simul.calculo_S(f_n)
    S_list.append(S)
  plt.plot(S_list, color=colores[i], label=estados[i])
plt.title('Entropia por tramos')
plt.xlabel('Tramos tiempo')
plt.ylabel('S')
ax.legend(loc='best')
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/entropia_wt_tramos.jpg',
            dpi=300)
plt.show()

#Representación de las distribuciones de los tiempos de residencia en cada tramo
fig, ax = plt.subplots(n, 1, sharex=True)
for j in range(n):
  for i in range(len(estados)):
    f_n = f_simul.mat_frec_tramos(dic, estados[i], j)
    ax[j].bar(f_n[1],f_n[0], color=colores[i], label=estados[i], alpha=0.5)
    ax[j].set_ylabel('N'+str(j))
fig.suptitle('Grafico de frecuencias')
ax[n-1].set_xlabel('tiempo de residencia')
plt.savefig('/Users/celia/OneDrive/Escritorio/BioCompLab/simulador/simulador-de-reacciones/plots_abc/distribucion_wt_tramos.jpg',
            dpi=300)
plt.show()