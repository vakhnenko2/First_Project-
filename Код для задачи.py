import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from matplotlib.animation import ArtistAnimation
from sympy import sin, tan, cos, pi, symbols, diff, Heaviside

x1,y1,D1,D2,Y1 = symbols('x1 y1 D1 D2 Y1')

# МОЖНО ИССЛЕДОВАТЬ ФОРМУ ГОРКИ
h = 5 # Высота подъема горки
D0 = - 0.01 * (x1-30)**2 + 6 # Уравнение параболы
L1 = 20
L2 = 40

#Функция возвращает первые и вторые производные (f1,f2) уравнений
#горки на разных участках и высоту горки
def mountain(x, L1, L2):
    if x < L1:
        f1,f2 = 0,0
        Y = h
    if x > L2:
        f1, f2 = 0, 0
        Y = h
    if L1 <= x <= L2:
        D1 = diff(D0,x1)
        D2 = diff(diff(D0, x1), x1)
        f1 = D1.subs(x1, x)
        f2 = D2.subs(x1, x)
        Y = D0.subs(x1, x)
    return f1, f2, Y

# print(mountain(25,L1,L2))

# Определяем функцию для системы диф. уравнений
def move_func(s, t):
    x, v_x, y, v_y = s

    res = mountain(x, L1, L2)

    f1 = res[0]
    f2 = res[1]

    #электрическая сила, действующая на частицу с зарядом q со стороны заряда Q,
    #движущегося из точки (a,b) со скоростью (Vx,Vt)
    Fe = q * Q / ((a + Vx * t - x)**2 + (b + Vy * t - y)**2)

    Fe_x = Fe*(a+Vx*t-x)/np.sqrt(((a+Vx*t-x)**2+(b+Vy*t-y)**2))
    Fe_y = Fe*(b+Vy*t-y)/np.sqrt(((a+Vx*t-x)**2+(b+Vy*t-y)**2))

    #множитель Лагранжа
    L = (f2*v_x**2 + g + f1*Fe_x-Fe_y)/(1+f1**2)

    #модуль скорости
    v = np.sqrt(v_x**2 + v_y**2)

    dxdt = v_x

    dv_xdt = -L*f1+Fe_x

    dydt = v_y

    dv_ydt = - g +Fe_y + L

    return dxdt, dv_xdt, dydt, dv_ydt
# Определяем начальные значения и параметры, входящие в систему диф. уравнений

x0 = 30
v_x0 = -7.89

y0 = 6
v_y0 = 0

g = 9.8

q = 1 #заряд движущейся частицы
Q = 100 #заряд частицы-источника поля

#первоначальные координаты заряда
a = 30
b = 4

#компоненты скорости источника
Vx = 0
Vy = 0

s0 = x0, v_x0, y0, v_y0

X = []
Y = []

N = 500
tau = np.linspace(-10,70,500)

Px=[]
Py=[]
PX=[]
PY=[]
T=np.linspace(0,30,N)

for i in range(N-1):
    t=[T[i],T[i+1]]
    sol = odeint(move_func, s0, t)

    PX.append(a+Vx*(t[1]))
    PY.append(b+Vy*(t[1]))
    Px.append(sol[1,0])
    Py.append(sol[1,2])

    x0 = sol[1,0]
    vx0 = sol[1,1]
    y0 = sol[1,2]
    vy0 = sol[1,3]

    s0 = x0, vx0, y0, vy0


# Построение фигуры
for i in range(500):
    X.append(tau[i])
    res = mountain(tau[i], L1, L2)
    Y.append(res[2])


fig = plt.figure()
body = []
step = int(N/100)

for i in range(0,len(T)-1,step):
    body1, = plt.plot(Px[i], Py[i], 'o', color='r')
    body1_line, = plt.plot(Px[:i], Py[:i],'-',color='r')
    body0, = plt.plot(PX[i], PY[i], 'o', color='g')

    body.append([body1, body1_line, body0])


ani = ArtistAnimation(fig,body,interval=5)

plt.plot(X,Y,color='b')
plt.xlim(-10,70)
plt.ylim(0,10)
plt.grid()

# plt.show()
ani.save('charge.gif')