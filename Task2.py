import numpy as np
import requests
import matplotlib.pyplot as plt
from scipy.constants import pi, speed_of_light
from scipy.special import spherical_jn as jn
from scipy.special import spherical_yn as yn

#Создание txt
def txt(x,y):
   f = open('result.txt','w')
   prob = ' '
   len1 = prob*4 + 'f' + prob*17 + 'RCS' + '\n'
   f.write(len1)
   for n in range(len(x)):
       f.write(str(x[n])+prob*4+str(y[n]) + '\n')
   f.close()

#Сферическая фуекция Бесселя 3-го рода
def hn(n, x):
    return jn(n, x) + 1j * yn(n, x)

def bn(n, x):
    return (x * jn(n - 1, x) - n * jn(n, x)) / (x * hn(n - 1, x) - n * hn(n, x))

def an(n, x):
   return jn(n, x) / hn(n, x)


#Выгрузка задания
url = requests.get('https://jenyay.net/uploads/Student/Modelling/task_02_01.txt')
par_str = url.text
number_str = 4
my_str = par_str.split('\n')[number_str].split()
D=float(my_str[1][my_str[1].find("=")+1:-1])
fmin=float(my_str[2][my_str[2].find("=")+1:-1])
fmax=float(my_str[3][my_str[3].find("=")+1:])
fstep = 1e+6
R = D / 2
f = np.arange(fmin, fmax, fstep)
lambd = speed_of_light / f
k = 2 * pi / lambd

#Подсчет ЭПР
arr_sum = [((-1) ** n) * (n + 0.5) * (bn(n, k * R) - an(n, k * R)) for n in range(1, 20)]
summa = np.sum(arr_sum, axis=0)
rcs = (lambd ** 2) / pi * (np.abs(summa) ** 2)

txt(f, rcs)

#График
plt.plot(f / 10e6, rcs)
plt.xlabel("$f, МГц$")
plt.ylabel(r"$\sigma, м^2$")
plt.grid()
plt.show()
