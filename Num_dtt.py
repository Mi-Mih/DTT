import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
import itertools
from scipy import integrate

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
start_solution=np.array([0,2,3,0,2,3,0,2,3,0,2,3,0,2,3,0,2,3,1,2,3,1,2,3,1,2,3])

#набор для символа леви-чевита
help=np.array([[0,1,2],[1,2,0],[2,0,1],[2,1,0],[1,0,2],[0,2,1],[0,0,0],[1,1,1],[2,2,2]])
#nks=itertools.product('012', repeat=3)

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

def choice_p(sol,i,j):
    return sol[3*(i+2)+j]

def first(n,k,s,sol,time):
    dpdt=choice_p(sol,k,n)*sol[k]*levi(n,k,s)
    return dpdt

def second(n,k,s,sol,time):
    global m
    global lambd
    dvdt=(1/m)*(-lambd*sol[s+3]-choice_p(sol,k,n)*sol[k]*kron(n,k)*levi(n,k,s))
    return dvdt

def third(n,k,s,sol,time):
    global lambd
    global R
    global M
    global d_alpha_dt
    global B_1
    global teta
    global teta_obr
    d_omega_dt=teta_obr[k][n]*kron(n,s)*(-lambd*(-R*choice_p(sol,k,s)*sol[s+3]*levi(2,s,n)*kron(k,n)+sol[s]*(1-R**2))-d_alpha_dt*B_1*sol[s]+choice_p(sol,k,s)*M*kron(k,2)-((sol[n])**2)*teta[k][n]*levi(n,k,s))
    return d_omega_dt

#Матрица правых частей диффур системы
def Right_part(solution,time):
    global nks
    F=np.array([])
    for i in range(9):
        dpdt=first(help[i][0],help[i][1],help[i][2],solution,time)
        dvdt=second(help[i][0],help[i][1],help[i][2],solution,time)
        d_omega_dt=third(help[i][0],help[i][1],help[i][2],solution,time)
        F=np.append(F,dpdt)
        F=np.append(F,dvdt)
        F=np.append(F,d_omega_dt)
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
t_n=5
time=np.arange(0,t_n,0.1)

#функция odeint, по аналогии с ode45
solution=odeint(Right_part,start_solution,time)

#Вытащим из solution нужные массивы

#Линейная скорость шара
V_1=solution[:,4]
V_2=solution[:,7]
V_3=solution[:,1]
V_mod=(V_1**2+V_2**2+V_3**2)**(1/2)

#Правая угловая скорость шара
spin_1=solution[:,2]
spin_2=solution[:,5]
spin_3=solution[:,8]
spin_mod=(spin_1**2+spin_2**2+spin_3**2)**(1/2)

#Тензор поворота
P_11=solution[:,18]
P_12=solution[:,6]
P_13=solution[:,12]
P_21=solution[:,9]
P_22=solution[:,21]
P_23=solution[:,24]
P_31=solution[:,3]
P_32=solution[:,15]
P_33=solution[:,2]
P=np.zeros((3,3))
spin=np.zeros((3))
#Угловая скорость шара, вычисляем по формуле omega=P*spin 
omega_1=np.array([])
omega_2=np.array([])
omega_3=np.array([])
for i in range(len(time)):
    P[0][0]=P_11[i]
    P[0][1]=P_12[i]
    P[0][2]=P_13[i]
    P[1][0]=P_21[i]
    P[1][1]=P_22[i]
    P[1][2]=P_23[i]
    P[2][0]=P_31[i]
    P[2][1]=P_32[i]
    P[2][2]=P_33[i]
    spin[0]=spin_1[i]
    spin[1]=spin_2[i]
    spin[2]=spin_3[i]
    om_help=P.dot(spin)
    omega_1=np.append(omega_1,om_help[0])
    omega_2=np.append(omega_2,om_help[1])
    omega_3=np.append(omega_3,om_help[2])
omega_mod=(omega_1**2+omega_2**2+omega_3**2)**(1/2)


omega_1_rot=np.array([])
omega_2_rot=np.array([])
omega_3_rot=np.array([])
omega_rot=np.zeros((3))
#Угловая скорость ротора omega_b
for i in range(len(time)):
    P[0][0]=P_11[i]
    P[0][1]=P_12[i]
    P[0][2]=P_13[i]
    P[1][0]=P_21[i]
    P[1][1]=P_22[i]
    P[1][2]=P_23[i]
    P[2][0]=P_31[i]
    P[2][1]=P_32[i]
    P[2][2]=P_33[i]
    omega_rot[0]=omega_1[i]
    omega_rot[1]=omega_2[i]
    omega_rot[2]=omega_3[i]
    om_help= np.array([omega_1[i],omega_2[i],omega_3[i]]) +  d_alpha_dt * P.dot(np.array([0,0,1]))
    omega_1_rot=np.append(omega_1_rot,om_help[0])
    omega_2_rot=np.append(omega_2_rot,om_help[1])
    omega_3_rot=np.append(omega_3_rot,om_help[2])
omega_rot_mod=(omega_1_rot**2+omega_2_rot**2+omega_3_rot**2)**(1/2)
#Когда найдём \vec(V_c(t)), проинтегриурем покоординатно и получим x(t), y(t), z(t)
#решаем символьно, экстраполируя V_1(t) V_2(t) V_3(t)........................................



#Строим графики

#omega_a(t)
fig1,ax1=plt.subplots()
ax1.plot(time,omega_mod,label='|omega_ball|(t)')
ax1.set_xlabel('time')
ax1.set_ylabel('|omega_ball|')
ax1.set_title('Зависимость угловой скорости шара от времени')
ax1.grid()
ax1.legend()

#V_c(t)
fig2,ax2=plt.subplots()
ax2.plot(time,V_mod,label='|V_c|(t)')
ax2.set_xlabel('time')
ax2.set_ylabel('|V_c|')
ax2.set_title('Зависимость скорости цента масс шара от времени')
ax2.grid()
ax2.legend()

#omega_rot(t)
fig3,ax3=plt.subplots()
ax3.plot(time,omega_rot_mod,label='|omega_rotor|(t)')
ax3.set_xlabel('time')
ax3.set_ylabel('|omega_rotor|')
ax3.set_title('Зависимость угловой скорости ротора от времени')
ax3.grid()
ax3.legend()

#x(t)
fig4,ax4=plt.subplots()
ax4.set_xlabel('time')
ax4.set_ylabel('X')
ax4.set_title('Зависимость X коорд центра масс шара от времени')
ax4.grid()
ax4.legend()

#y(t)
fig5,ax5=plt.subplots()
ax5.set_xlabel('time')
ax5.set_ylabel('Y')
ax5.set_title('Зависимость Y коорд центра масс шара от времени')
ax5.grid()
ax5.legend()

#z(t)
fig6,ax6=plt.subplots()
ax6.set_xlabel('time')
ax6.set_ylabel('Z')
ax6.set_title('Зависимость Z коорд центра масс шара от времени')
ax6.grid()
ax6.legend()

plt.show()
