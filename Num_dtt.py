import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
import itertools
#входные параметры
m=100
lambd=0.1
R=100
M=17500
B_1=350
B_2=300
d_alpha_dt=45
teta=np.array([[((2/5)*m*R**2)/5+B_2,0,0],[0,((2/5)*m*R**2)/5,0],[0,0,B_1+((2/5)*m*R**2)/5]])
teta_obr=np.linalg.inv(teta)
#начальные условия
start_solution=np.array([0,0,0,0,0,0,1,0,0,0,1,0,0,0,1]*27)
nks=itertools.product('012', repeat=3)

def levi(k,n,s):
    if n==k or k==s or s==n:
        return 0
    elif (k,n,s)==(0,1,2) or (k,n,s)==(1,2,0) or (k,n,s)==(2,0,1):
         return 1
    elif (k,n,s)==(2,1,0) or (k,n,s)==(1,0,2) or (k,n,s)==(0,2,1):
         return -1 
def kron(n,k):
    if n==k:
       return 1
    elif n!=k:
       return 0
def first(n,k,s,sol,time):
    dpdt=sol[6+n+k]*sol[k]*levi(n,k,s)
    return dpdt
def second(n,k,s,sol,time):
    global m
    global lambd
    dvdt=(1/m)*(-lambd*sol[3+s]-sol[6+n+k]*sol[k]*kron(n,k)*levi(n,k,s))
    return dvdt
def third(n,k,s,sol,time):
    global lambd
    global R
    global M
    global d_alpha_dt
    global B_1
    global teta
    global teta_obr
    d_omega_dt=teta_obr[k][n]*kron(n,s)*(-lambd*(-R*sol[6+k+s]*sol[3+s]
                                                 *levi(2,s,n)*kron(k,n)+sol[s]*(1-R**2))-d_alpha_dt*B_1*sol[s]+sol[6+k+s]
                                         *M*kron(k,2)-(sol[n])**2*teta[k][n]*levi(n,k,s))
    return d_omega_dt
#Матрица правых частей диффур системы
def Right_part(solution,time):
    global nks
    F=np.array([])
    for i in nks:
        y=list(''.join(i))
        help = np.array([int(i) for i in y])
        print(help)
        dpdt=first(help[0],help[1],help[2],solution,time)
        dvdt=second(help[0],help[1],help[2],solution,time)
        d_omega_dt=third(help[0],help[1],help[2],solution,time)
        F=np.append(F,dpdt)
        F=np.append(F,dvdt)
        F=np.append(F,d_omega_dt)
        #print(F)
    return F
    '''
    правая часть 1 уравнения
    правая часть 2 уравнения
    правая часть 3 уравнения
    .........................
    правая часть N уравнения
    return np.array([N выражений из правой части диффур])
'''
#Массив времени, где t_0=0, a t_n=любое
t_n=5.1
time=np.arange(0,t_n,0.1)
#сама функция odeint, по аналогии с ode45
#solution=odeint(Right_part,start_solution,time)

#Когда найдём \vec(V_c(t)), проинтегриурем покоординатно и получим x(t), y(t), z(t)
#........................................

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
