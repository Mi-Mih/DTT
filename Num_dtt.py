import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math

#начальные условия
start_solution=np.array([Начальные условия, в том числе и для тензора поворота])

#Матрица правых частей диффур системы
def Right_part(solution,time):
    правая часть 1 уравнения
    правая часть 2 уравнения
    правая часть 3 уравнения
    .........................
    правая часть N уравнения
    return np.array([N выражений из правой части диффур])

#Массив времени, где t_0=0, a t_n=любое
t=np.arange(0,t_n,0.1)

#сама функция odeint, по аналогии с ode45
solution=odeint(F,s0,t)

#Когда найдём \vec(V_c(t)), проинтегриурем покоординатно и получим x(t), y(t), z(t)
........................................

#Строим графики

#omega_a(t)
fig1,ax1=plt.subplots()
ax1.grid()
ax1.legend()

#V_c(t)
fig2,ax2=plt.subplots()
ax2.grid()
ax2.legend()

#x(t)
fig2,ax2=plt.subplots()
ax2.grid()
ax2.legend()

#y(t)
fig2,ax2=plt.subplots()
ax2.grid()
ax2.legend()

#z(t)
fig2,ax2=plt.subplots()
ax2.grid()
ax2.legend()
