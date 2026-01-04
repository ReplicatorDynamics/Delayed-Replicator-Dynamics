import math

import numpy as np
import matplotlib.pyplot as plt
import symengine

from jitcdde import jitcdde, y, t

res =50
duration =15000    #30000
delay = 50
   #183
F =0.85                  #Fertility reward
m = 1                 #mortality cost
K = 50000
alpha = 0 # amplitude of seasonal mortality
theta = 50 #duration of seasonal mortality
Phi = 0
Psi =0.1253      #0.125*m+0.001-alpha     #.158
Omega = 0.000
bp = 0     # 0.3
dp = 0     #  0.8
pd = 0    #   0.6



frac = 0.6    #  base of the birth mortality
OOO = 50000   #   scale parameter, evolving prey population size at which juvenile survival equals frac

d = symengine.symbols("d")

bm = frac ** (y(1, t - d) / OOO)
logistic = (1 - (y(1, t) / K))
dellogistic = (1 - (y(1, t - d) / K))
suppr = bm

#initial history

initial = 0  #0-constant history, 1 function history

qinit = 0.7   #initial conditions, dove frequency
ninit = 45000    #prey population size
xinit = 0   #pedator population size




if initial == 1:
 def initial_history_func_q(t):
  return [qinit + 0.0001 * t, ninit - 0.0001 * t, xinit + 0 * t]


#def seasons(t) :
 #    return[alpha * math.sin(((2 * math.pi) / theta) * t)]
 #   (y(1, t - d) / y(1, t)) * (0.5 * F * (y(0, t - d) ** 2 - y(0, t)) + (y(0, t - d) - y(0, t)) * Phi) * suppr + 0.5 * m * y(0, t) * ((1 - y(0, t)) ** 2)
 #     y(1, t - d) * (0.5 * F + Phi) * suppr - y(1, t) * (0.5 * m * (1-y(0, t)) ** 2 + Psi + alpha * symengine.sin(((2 * math.pi) / theta) * t) + Omega * y(1, t) + pd * y(2, t))
 #      y(2, t - d) * bp * y(1, t - d) - y(2, t) * dp
equations=[
    (y(1, t - d) / y(1, t)) * (0.5 * F * (y(0, t - d) ** 2 - y(0, t)) + (y(0, t - d) - y(0, t)) * Phi) * suppr + 0.5 * m * y(0, t) * ((1 - y(0, t)) ** 2)
    ,
    y(1, t - d) * (0.5 * F + Phi) * suppr - y(1, t) * (0.5 * m * (1-y(0, t)) ** 2 + Psi + alpha - alpha * symengine.sin(((2 * math.pi) / theta) * t) + Omega * y(1, t) + pd * y(2, t)),
     y(2, t) * bp * y(1, t) - y(2, t) * dp]

ddesys = jitcdde(equations, control_pars=[d], max_delay=800.)

plt.rcParams['font.size'] = 8
fig, axs = plt.subplots(2, 1)
fig.tight_layout(rect=[0, 0, 1, 0.95], pad=3.0)
#fig.suptitle(" solved by jitcdde")

ts = np.linspace(0, duration, res * duration)

if initial == 1:
 ddesys.past_from_function(initial_history_func_q)

if initial == 0:
 ddesys.constant_past([qinit , ninit, xinit])

params = [delay]
ddesys.set_parameters(*params)
ys = []
for t in ts:
    ys.append(ddesys.integrate(t))
ys=np.array(ys)

ndif = ys[:, 1].copy()

count = delay
while (count < res * duration):
    ndif[count] = ys[count - delay, 1] / ys[count, 1]
    count = count + 1

qdif = ys[:, 0].copy()

count = delay
while (count < res * duration):
    qdif[count] = ys[count - delay, 0] - ys[count, 0]
    count = count + 1



#nullclines and subnullclines

if xinit == 0 and suppr == logistic and Omega == 0 :

 qs = np.linspace(0, 1, 1000)
 ns = np.linspace(0, K, 1000)

   #nullclines

 qnull = 1 - ((1 - (ns / K)) * F) / m

 nnull = K * (1 - (Psi + alpha + 0.5 * m * (1 - qs) ** 2) / (0.5 * F + Phi))

     #subnullclines

 Anull = ((0.5 * F + Phi) ** 2)/K
 Bnull = -((0.5 * F + Phi) ** 2 - (0.5 * F + Phi) * (((1 - qs) ** 2) * 0.5 * m + Psi + alpha) - 0.5 * qs * ((1 - qs) ** 2) * F * m)
 Cnull = 0.5 * qs * ((1 - qs) ** 2) * m * K * ((1 - qs) * m - F)

 DELTAnull = Bnull ** 2 - 4 * Anull * Cnull

 SQDELTA = (DELTAnull)  ** 0.5

 NnullB = (-Bnull / (2 * Anull)) - ((SQDELTA) / (2 * Anull))

     #q subnullcline

 Aqnull = ((0.5 * F + Phi) * F) / (m * (K ** 2))
 Bqnull = (0.5 * qs * (1 - qs) - ((0.5 * F + Phi) - (1 - 2 * qs + qs ** 2) * 0.5 * m - Psi - alpha) / m) * (F / K)
 Cqnull = 0.5 * ((qs - 2 * (qs ** 2) + qs ** 3) * m - qs * (1 - qs) * F)

 DELTAQnull = (Bqnull) ** 2 - 4 * (Aqnull) * (Cqnull)

 Nqnullslope = -(Bqnull) / (2 * Aqnull) - ((DELTAQnull) ** 0.5) / (2 * Aqnull)

fig1 = plt.figure(1)
axs[0].plot(ts, ys[:, 0], color='black', linewidth=0.4)      #ys[:, 0]
axs[0].plot(ts, ys[:, 2], color='grey', linewidth=0.4)
if alpha != 0:
   axs[0].plot(ts,  0.01 + 0.01 * np.sin(((2 * math.pi) / theta) * ts), color='red', linewidth=0.4)
#axs[0].plot(ts, ys[:, 2] -K + 1, color='red', linewidth=1)
#axs[0].plot(ts, ts * 0 + 1, color='green', linewidth=1)
#axs[0].plot(ts, ts * 0, color='green', linewidth=1)
axs[0].set_title('Dove strategy frequency')


axs[1].plot(ts, ys[:, 1], color='black', linewidth=0.4)
axs[1].plot(ts, ys[:, 2], color='grey', linewidth=0.4)

 #  axs[1].plot(ts,  alpha + alpha * np.sin(((2 * math.pi) / theta) * ts), color='red', linewidth=0.4)
#axs[1].plot(ts, ts * 0, color='green', linewidth=1)
axs[1].set_title('population size')

#axs[2].plot(ts, frac ** (ys[:,1] / OOO), color='blue', linewidth=1)
#axs[2].plot(ts, (1 - (ys[:,1] / K)), color='red', linewidth=1)
#axs[2].plot(ts, ts * 0 + 1, color='green', linewidth=1)
#axs[2].plot(ts, ts * 0, color='green', linewidth=1)
#axs[2].set_title('mortalities')





fig2 = plt.figure(2)
plt.plot(ys[:, 1], ys[:, 0], color='black', linewidth=0.4)
if xinit == 0 and suppr == logistic and Omega == 0 :
 plt.plot(ns[:], qnull[:], color='red', linewidth=0.4)
 plt.plot(nnull[:], qs[:], color='blue', linewidth=0.4)
 plt.plot(NnullB[:], qs[:], color='violet', linewidth=0.4)
 plt.plot(Nqnullslope[:], qs[:], color='orange', linewidth=0.4)
plt.title('Dove frequency/population size phase portrait')
plt.xlabel('population size')
plt.ylabel('Dove frequency')
plt.ylim([0, 1])

fig3 = plt.figure(4)
plt.plot(ys[:, 1], ys[:, 2], color='black', linewidth=0.4)
plt.title('prey vs. predator population size phase portrait')
plt.xlabel('prey population size')
plt.ylabel('predator population size')

fig5 = plt.figure(5)

# syntax for 3-D projection
ax = plt.axes(projection='3d')

# defining all 3 axis
z = ys[:, 0]
x = ys[:, 1]
y = ys[:, 2]

# plotting
ax.plot3D(x, y, z, 'black', linewidth=0.4)
#ax.set_title('3D line plot geeks for geeks')
ax.set_xlabel('prey population size', fontsize=12)
ax.set_ylabel('predator population size', fontsize=12)
ax.set_zlabel('Dove frequency', fontsize=12)

fig6 = plt.figure(6)
plt.plot(ts, qdif, color='black', linewidth=0.4)
#axs[0].plot(ts, ys[:, 2] -K + 1, color='red', linewidth=1)
#axs[0].plot(ts, ts * 0 + 1, color='green', linewidth=1)
#axs[0].plot(ts, ts * 0, color='green', linewidth=1)
#set_title('Dove strategy frequency')

axs[1].plot(ts, ys[:, 1], color='black', linewidth=0.4)
axs[1].plot(ts, ys[:, 2], color='grey', linewidth=0.4)
#axs[1].plot(ts, ts * 0, color='green', linewidth=1)
axs[1].set_title('population size')

fig7 = plt.figure(7)
plt.plot(ts, ndif, color='black', linewidth=0.4)

plt.show()

print(ndif)