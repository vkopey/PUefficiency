import matplotlib.pyplot as plt
import sympy as sp
import numpy as np
from scipy.optimize import curve_fit

# Визначаємо символьні змінні
t = sp.symbols('t', real=True, positive=True)
m, c, k, A, omega, v = sp.symbols('m c k A omega v', real=True, positive=True)
x = sp.Function('x')(t)

# 1. Записуємо диференціальне рівняння
# m*x''(t) + c*x'(t) + k*x(t) = k * A * sin(omega * t)
equation = sp.Eq(m * x.diff(t, t) + c * x.diff(t) + k * x, k * A * sp.sin(omega * t))

# 2. Розв'язуємо рівняння в символьному вигляді
solution = sp.dsolve(equation, x)
solution1=sp.simplify(solution.subs({'C1':0,'C2':0}))

fv={m:3961,k:44650,c:2121,A:2.1/2,omega:0.3}#6.4/60}
OMEG=0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0
AMP1=1.05, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.13, 1.15
fv={m:3961,k:44650/4,c:2121,A:2.1/2,omega:0.3}
AMP2=1.06, 1.07, 1.09, 1.11, 1.15, 1.20, 1.26, 1.33, 1.44, 1.56

def regr(x,a,b):
    return a*x**b

X=np.array(OMEG)
Y=np.array(AMP2)-np.array(AMP1)
popt,_=curve_fit(regr, OMEG, Y)
print(popt)
print("R^2:", np.corrcoef(Y, regr(X,*popt))[0,1]**2 )

plt.plot(OMEG,AMP1,'ks--',OMEG,AMP2,'ko-')
plt.xlabel('$\omega$, рад/c');plt.ylabel("A, м");plt.grid()
plt.show()



x_t=solution1.rhs.subs(fv)
print((44650/3961)**0.5)
print((0.25*44650/3961)**0.5)
sp.plot(x_t, (t,0,500))#, show=False
#force=k * A * sp.sin(omega * t) # ?
u_t = A * sp.sin(omega * t) # Визначаємо функцію переміщення верхньої точки
F_dyn = k * (u_t - x_t) # Сила пружності (динамічна складова)
sp.plot(F_dyn.subs(fv), (t,0,500))
sp.plot_parametric(x_t, F_dyn.subs(fv), xlim=[-1.1, 1.1], ylim=[-800, 800], line_color='red') # динамограма




"""
# 2. Визначаємо початкові умови
# x(0) = 0
# x'(0) = v
initial_conditions = {x.subs(t, 0): 0, x.diff(t).subs(t, 0): v}

# 3. Розв'язуємо рівняння з урахуванням початкових умов
solution = sp.dsolve(equation, x, ics=initial_conditions)

print("Розв'язок задачі Коші (з початковими умовами):")
sp.pprint(solution)

# Додатково: виведемо тільки праву частину розв'язку для подальшої роботи
x_t = solution.rhs
print(sp.simplify(x_t))
sp.plot(x_t.subs({k:1, A:1, c:0.1, v:1, omega:0.01, m:1}))
"""