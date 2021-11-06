import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math as ma


class black_hole():
    def __init__(self, alpha0, eta0, nu0):
        self.__c = 6.407e4
        self.__G = 39.309
        self.alpha = alpha0
        self.eta = eta0
        self.nu = nu0

    def masa(self):
        if self.alpha < 0:
            '''calcular masa branch positivo'''
            M = -self.__c ** 2 / self.__G * \
                (3 * self.eta ** 2 + self.alpha) / self.eta ** 3 / 6
            return M
        elif self.alpha > 0:
            # calcular masa branch negativo
            M = self.__c ** 2 / self.__G * \
                (3 * self.eta ** 2 + self.alpha) / self.eta ** 3 / 6
            return M
        else:
            print("alpha no puede ser 0")

    def horizonte_hairy_x(self):
        '''calcular x horizonte'''
        if self.alpha > 0:
            x_inicial = 0.00004
            return fsolve(self.f, x_inicial)
        elif self.alpha < 0:
            x_inicial = 26
            return fsolve(self.f, x_inicial)

    def horizonte_hairy(self):
        if self.alpha > 0:
            x_inicial = 0.00004
            x_h = fsolve(self.f, x_inicial)
            return self.omega(x_h)**0.5
        elif self.alpha < 0:
            x_inicial = 26
            x_h = fsolve(self.f, x_inicial)
            return self.omega(x_h)**0.5

    def horizonte(self):
        r_h = 2*self.__G*self.masa()/self.__c**2
        return r_h

    def omega(self, x):
        O = self.nu ** 2 * x ** (self.nu - 1) / \
            self.eta ** 2 / (x ** self.nu - 1) ** 2
        return O

    def f(self, x):
        f = self.alpha * (1 / (self.nu ** 2 - 4) - x ** 2 * (1 + x ** (-self.nu) / (self.nu - 2) - x ** self.nu / (
            self.nu + 2)) / self.nu ** 2) + x * self.eta ** 2 * (x ** self.nu - 1) ** 2 / self.nu ** 2 / x ** (self.nu - 1)
        return f

    def diff_omega(self, x):
        dO = self.nu ** 2 * x ** (self.nu - 1) * (self.nu - 1) / x / self.eta ** 2 / (x ** self.nu - 1) ** 2 - \
            2 * self.nu ** 3 * \
            x ** (self.nu - 1) / self.eta ** 2 / \
            (x ** self.nu - 1) ** 3 * x ** self.nu / x
        return dO

    def diff_f(self, x):
        df = self.alpha * (-2 * x * (1 + x ** (-self.nu) / (self.nu - 2) - x ** self.nu / (self.nu + 2)) / self.nu ** 2 - x ** 2 * (-x ** (-self.nu) * self.nu / x / (self.nu - 2) - x ** self.nu * self.nu / x / (self.nu + 2)) / self.nu ** 2) + self.eta ** 2 * (
            x ** self.nu - 1) ** 2 / self.nu ** 2 / x ** (self.nu - 1) + 2 * self.eta ** 2 * (x ** self.nu - 1) / self.nu / x ** (self.nu - 1) * x ** self.nu - self.eta ** 2 * (x ** self.nu - 1) ** 2 / self.nu ** 2 / x ** (self.nu - 1) * (self.nu - 1)
        return df

    def estado(self):
        return "La masa del BH es: ", self.masa(), "\nEl horizonte hairy es: ", self.horizonte_hairy(), "\nEl horizonte tipo Schwarzschild es: ", self.horizonte()


class particula_time_like(black_hole):
    def __init__(self, alpha0, eta0, nu0, energia0, J0):
        super().__init__(alpha0, eta0, nu0)
        self.__c = 6.407e4
        self.__G = 39.309
        self.energia = energia0
        self.J = J0

    def U_potencial(self, x):
        U = self.omega(x)*self.f(x)*(1+self.J**2*self.__c**2/(self.omega(x)))-1
        return U

    def potencial(self):
        if self.alpha < 0:
            x = np.linspace(1, 1000, 100000)
            U = self.U_potencial(x)
            r = self.omega(x)**(0.5)
            return r, U
        elif self.alpha > 0:
            x = np.linspace(self.horizonte_hairy(), 1, 10000)
            U = self.U_potencial(x)
            r = self.omega(x)**(0.5)
            return r, U
        else:
            print("alpha no puede ser 0")

    def diff_U(self, x):
        dU = self.diff_f(x)*self.omega(x)+self.diff_omega(x) * \
            self.f(x)+(self.J*self.__c)**2*self.diff_f(x)
        return dU

    def maximo(self):
        x_inicial = float(self.horizonte_hairy_x())
        x_max = fsolve(self.diff_U, x_inicial)
        U_max = self.omega(x_max)*self.f(x_max)*(1+self.J **
                                                 2*self.__c**2/(self.omega(x_max)))-1
        return self.omega(x_max)**(0.5), U_max

    def minimo(self):
        x_inicial = float(self.horizonte_hairy_x())+0.8
        x_min = fsolve(self.diff_U, x_inicial)
        U_min = self.omega(x_min)*self.f(x_min)*(1+self.J **
                                                 2*self.__c**2/(self.omega(x_min)))-1
        return self.omega(x_min)**(0.5), U_min

    def U_nulo(self):
        x_inicial = float(self.horizonte_hairy_x())+0.6
        x_nulo = fsolve(self.U_potencial, x_inicial)
        return self.omega(x_nulo)**(0.5), self.U_potencial(x_nulo)

    def cond_init(self, r):
        if self.alpha > 0:
            x_initial_guess = 0.83
            def O_r(x): return self.nu ** 2 * x ** (self.nu - 1) / \
                self.eta ** 2 / (x ** self.nu - 1) ** 2-r**2
            x_radio = fsolve(O_r, x_initial_guess)
            x_prima = ma.sqrt((self.energia-self.U_potencial(x_radio)
                            )/(self.eta**2*self.J**2*self.__c**2))
            return x_radio, x_prima
        elif self.alpha < 0:
            x_initial_guess = 1.22
            def O_r(x): return self.nu ** 2 * x ** (self.nu - 1) / \
                self.eta ** 2 / (x ** self.nu - 1) ** 2-r**2
            x_radio = fsolve(O_r, x_initial_guess)
            x_prima = ma.sqrt((self.energia-self.U_potencial(x_radio)
                            )/(self.eta**2*self.J**2*self.__c**2))
            return x_radio, x_prima
        else:
            print("Alpha no puede ser 0")

    def funH(self, x):
        H = self.diff_U(x)/(2*self.J**2*self.__c**2*self.eta**2)
        return H

    def x_turn_r(self, x):
        return self.omega(x)**0.5


class particula_null(black_hole):
    def __init__(self, alpha0, eta0, nu0, b0):
        super().__init__(alpha0, eta0, nu0)
        self.__c = 6.407e4
        self.__G = 39.309
        self.b = b0

    def U_potencial(self, x):
        U = self.__c**2*self.f(x)
        return U

    def potencial(self):
        if self.alpha < 0:
            x = np.linspace(1, 1000000, 1000000)
            U = self.U_potencial(x)
            r = self.omega(x)**(0.5)
            return r, U
        elif self.alpha > 0:
            x = np.linspace(0, 1, 1000000)
            U = self.U_potencial(x)
            r = self.omega(x)**(0.5)
            return r, U
        else:
            print("alpha no puede ser 0")

    def maximo(self):
        #x_inicial = float(self.horizonte_hairy_x())
        x_inicial = 0.01  #momentaneo
        x_max = fsolve(self.diff_f, x_inicial)
        return self.omega(x_max)**(0.5), self.U_potencial(x_max)

    def cond_init(self, r_x):
        r = ma.sqrt(r_x**2+self.b**2)
        theta = ma.atan(self.b/r_x)
        x_initial_guess = 0.8
        def O_r(x): return self.nu ** 2 * x ** (self.nu - 1) / \
            self.eta ** 2 / (x ** self.nu - 1) ** 2-r**2
        x_radio = fsolve(O_r, x_initial_guess)
        x_prima = ma.sqrt((1/self.b**2-self.f(x_radio))/(self.eta**2))
        return x_radio, x_prima, theta

    def x_turn_r(self, x):
        return self.omega(x)**0.5

    def funH(self, x):
        H = self.diff_f(x)/(2*self.eta**2)
        return H


def RK(x_n, y_n, h, s, particula):
    r = np.array([])

    # punto de referencia
    for i in range(s):
        r_n = particula.omega(x_n)**0.5
        r = np.append(r, r_n)

        k_0 = h*y_n
        l_0 = h*(-particula.funH(x_n))

        k_1 = h*(y_n+l_0/2)
        l_1 = h*(-particula.funH(x_n+0.5*k_0))

        k_2 = h*(y_n+l_1/2)
        l_2 = h*(-particula.funH(x_n+0.5*k_1))

        k_3 = h*(y_n+l_2)
        l_3 = h*(-particula.funH(x_n+k_2))

        y_n1 = y_n+(l_0+2*l_1+2*l_2+l_3)/6
        x_n1 = x_n+(k_0+2*k_1+2*k_2+k_3)/6

        y_n = y_n1
        x_n = x_n1
    return(r)   