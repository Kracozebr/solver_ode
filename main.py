import numpy as np
from rk.rk4 import RungeKutta4

x1 = 0

def M1(n):
    a0 = 28.44 
    a1 = -0.27138
    a2 =  1.00783 * 10 **(-3)
    a3 = -1.56667 * 10 **(-6)
    a4 = 8.66667 * 10 **(-10)

    return a0 + a1 * n + a2 * (n**2) + a3 * (n**3) + a4 * (n**4)

fi0 = 28 * np.pi /180 # угол канавки
d01 = 0.055 # м начальный диаметр ведущего шкива
d02 = 0.179 # м начальный диаметр ведомого шкива


a = 0.2

# Захардокили L иначе закольцовыввются переменные
L = 2 * a + np.pi / 2 * (d01+d02)+((d02-d01) **2) / 4 * a

def get_dx1(x1):
    return d01+x1/np.tan(fi0/2) # диаметр расопложения ремня на ведущем шкиве


def R(x1):
    # print(get_dx1(x1))
    return np.sqrt(np.pi * a * (np.pi * a)/4 + L/np.pi - (2 * a)/np.pi - get_dx1(x1))


def get_x2(x1):
    return (np.pi * a + d02 - d01 - 2 * R(x1)) * np.tan(fi0 / 2) - x1


def get_dx2(x1):
    return d02-get_x2(x1)/np.tan(fi0/2) # диаметр расопложения ремня на ведомом шкиве


G_sum = 52 + 75 # кгс  

r_k = 0.265 # м радиус колеса
i_dp =  4 # передаточное число доп передачи
niu_dp = 0.9 # кпд доп передачи

f_c = 0.02 # коэфф сопротивления качению
alfpha_n = 0 # угол наклона дороги к горизонту


k_b = 0.07 # коэфф сопротивления воздуха
s = 0.5 # м2 площадь проекции транспортного средства


A = 0.8 # безразмерный коэффициент нажимных механизмов
J1 = 0.0005161 # ??кгс м2  момент инерции приведенный к ведущему валу вариатора
J02 = 0.00050342 # Собственный момент инерции ведомого шкива
J_k = 0.0191098 # Момент инерции 2 колес
G_m = 52 #  вес транспортного средства
G_n = 75 # вес водителя
g = 9.8
J2 = J02 + 1 / (i_dp**2 * niu_dp) * ((G_m / g + G_n/ g) *r_k**2 + 2 * J_k)
r0 = 0.0325 #м  начальный радиус центра тяжести шариков 
m_sh = 0.0056 # кг масса шариков

def get_DIAM(x1):
    dx1 = get_dx1(x1)
    x2 = get_x2(x1)
    dx2 = get_dx2(x2)

    return dx1/dx2

GR=G_sum * (r_k / (i_dp * niu_dp)) * (f_c * np.cos(alfpha_n) + np.sin(alfpha_n))
KR=((k_b * s * r_k**3)/((i_dp**3) * niu_dp))


def get_DP(x1):
    return ((1 / get_dx2(x1)) + (get_dx1(x1) / get_dx2(x1)**2) * (((np.pi * a) / R(x1)) - 1))


def get_RA(x1):
    return r0 + A * x1


def get_RA(x1):
    return (r0 + A * x1)

m1 = 0.03 # кгс вес подвижного диска ведущего шкива
m2 = 0.023 # кгс вес подвижного дискка ведомого шкива
c_pr = 840  # суммарная жесткость пружин ведомого шкива
x02 = 0.026 # начальнаяа деформация пружин
D = np.pi * a + d02 - d01 
dc1 = 0.035 # диаметр ступицы ведущего шкива
dc2 = 0.108 # диаметр ступицы ведомого шкива


def get_gamma(x1):
    return np.tan(((get_dx2(x1) - get_dx1(x1)) / (2 * a)) / (np.sqrt (1 - (((get_dx2(x1) - get_dx2(x1)) / 2 * a) **2))))  # угол между линией центров и ветвями ремня


alpha1 = np.pi - 2 * get_gamma(x1) # угол обхвата меньшего шкива
alpha2 = np.pi + 2 * get_gamma(x1) # угол обхвата большего шкива
alpha_c = 0.7 * alpha1 # угол скольжения от 30 вроде до 80
f_prime = 0.3 # f * sin(fio) приведенный коэффициент трения ремня и шкива 0.041
m = np.exp(alpha_c*f_prime) # m= 3.82 ... 1.77
alpha_n1 = 30 * np.pi /180 # угол покоя ведущего шкива
ro = 10 * np.pi /180 # угол трения
ro_R = 60 *  np.pi /180 # угол радиальной составляющей силы трения
alpha_n2 = 30 * np.pi /180 # угол покоя ведомого шкива
psi = m-1/m+1

lc1 = 0.022 # длина ступицы ведущего шкива
lc2 = 0.037 # длина ступицы ведомого шкива

f = 0.319 # коэффициент трения пары ремень-шкива
f_tr = 0.015 # коэффициент трения подвижной пары

def get_P(fi, x1):
    return 2 * M1(fi) / get_dx1(x1) # кгс окружная сила



Y1 = ((np.cos(fi0/2))/2*f) + (m/(m-1))*(alpha_n1 / (2 * np.tan(fi0/2+ro))) # относительная осевая сила ведущей ветви
Y2=(1 / 2 * f_prime * np.tan((fi0/2)+ro_R)) + (1/ m - 1) * (alpha_n2 / (2 * np.tan((fi0 / 2) + ro))) # относительная осевая сила ведомой ветви
Tetta = Y1/Y2 # отношение потребных относительных осевых сил вариатора

# np
def get_MP(x1):
    return (1 / 2) * m2 * ((np.pi**2 * a**2) / R(x1)**3) * (((np.pi * a) / R(x1)) - 1) * 1/np.tan(fi0 / 2) # исправил 1/tg на 1/cot

def get_CP(x1):
    return c_pr * (x02 + ((D - 2 * R(x1)) * np.tan(fi0 / 2)) - x1)


def get_TD(x1):
    return (get_dx1(x1) / dc1) + Tetta * (get_dx2(x1) / dc2) + ((np.cos(get_gamma(x1))) / psi) * (1 + Tetta)


def get_YD(x1):
    return ((get_dx1(x1) / lc1) * abs((np.sin(alpha1) / alpha1)) + (get_dx2(x1) / lc2) * abs((np.sin(alpha2) / alpha2)))


def get_MS(x1):
    return (m1 + m2 * ((((np.pi * a)) / R(x1)) - 1)**2 + m_sh * A**2)


def coupled_system_of_odes(t, vars):
        omega, fi, v, x1 = vars

        omega_dot = (M1(fi) - (get_DIAM(x1)) * GR - (get_DIAM(x1))**3 * KR * omega**2 - J2 * (fi0 / 2) * \
                (get_DIAM(x1)) * (get_DP(x1)) * omega * v - 2 * m_sh * A * (get_RA(x1)) * omega * v) / \
                (J1 + J2 * (get_DIAM(x1))**2 + m_sh * (get_RA(x1))**2) 

        fi_dot = omega

        v_dot = (m_sh * A * (get_RA(x1)) * omega**2 - (get_MP(x1)) * v - Tetta * (get_CP(x1)) - np.sign(v) * \
                get_P(fi, x1) * f_tr * ((get_TD(x1)) + Y1 * (get_YD(x1)))) / (get_MS(x1))

        x1_dot = v
        print(f"t={t} x1={x1} x2={get_x2(x1)} dx1={get_dx1(x1)} dx2={get_dx2(x1)} omega_dot={omega_dot} omega={omega} v_dot={v_dot} x1_dot={x1_dot}")

        return np.array([omega_dot, fi_dot, v_dot, x1_dot])


if __name__ == "__main__":
    # Initial conditions
    initial_conditions = [360.0, 0.0, 0.0, 0.0] #omega, fi, v, x1
    t0 = 0
    tf = 1
    h = 0.001
    print("hello")

    # Create an instance of the RungeKutta4 class
    rk4_solver = RungeKutta4(coupled_system_of_odes, initial_conditions, t0, tf, h)

    # Solve the system
    t_values, y_values = rk4_solver.solve()