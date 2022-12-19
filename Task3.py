'''
Размер области моделирования 2,5 метров
Диэлектрическая постонная 6
Гармонический сигнал
Слева: ABC 1 порядка
Справа: PEC
'''
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft, fftshift


def Signal(m,q, name):
    if str(name) == 'sin':
        # Параметры сигнала синуса
        N_la = 50
        phi = -2 * np.pi / N_la
        return np.sin(2*np.pi / N_la * (Sc*q - np.sqrt(mu * eps) * m) + phi)
    if str(name) == 'gauss':
        return np.exp(-(((q - m * np.sqrt(eps * mu) / Sc) - 100.0) / 30.0) ** 2)
    if str(name) == 'modgauss':
        N_la = 50
        phi = -2 * np.pi / N_la
        return np.sin(2*np.pi / N_la * (Sc*q - np.sqrt(mu * eps) * m) + phi) * np.exp(-(((q - m * np.sqrt(eps * mu) / Sc) - 100.0) / 30.0) ** 2)
    if str(name) == 'rect':
        return 1 if 0 <= q - m <= 20 else 0

if __name__ == '__main__':
    # Название сигнала
    name_signal = 'sin'
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * np.pi

    # Число Куранта
    Sc = 1.0

    # Размер области в см
    l = 250

    # Время расчета в отсчетах
    maxTime = 1000

    # Размер области моделирования в отсчетах
    maxSize = 500

    # Шаг по пространству
    dx = 1e-3

    # Шаг по времени
    c = 3*10**8
    dt = Sc * dx / c

    # Положение датчика, регистрирующего поле
    probePos = 300

    # Положение источника
    sourcePos = 250

    # Параметры поля
    mu = 1.0
    eps = 6.0

    Ez = np.zeros(maxSize)
    Hy = np.zeros(maxSize - 1)

    oldEzLeft = Ez[1]

    tempLeft = Sc / np.sqrt(mu * eps)
    koeffABCLeft = (tempLeft - 1) / (tempLeft + 1)

    # Поле, зарегистрированное в датчике в зависимости от времени
    probeTimeEz = np.zeros(maxTime)
    probeTimeEz[0] = Ez[probePos]

    # Подготовка к отображению поля в пространстве
    xlist = np.arange(maxSize)

    # Включить интерактивный режим для анимации
    plt.ion()

    # Создание окна для графика
    fig, ax = plt.subplots()

    # Установка отображаемых интервалов по осям
    ax.set_xlim(0, maxSize * dx)
    ax.set_ylim(-1.1, 1.1)

    # Установка меток по осям
    ax.set_xlabel('x, см')
    ax.set_ylabel('Ez, В/м')

    # Включить сетку на графике
    ax.grid()

    # Отобразить поле в начальный момент времени
    line, = ax.plot(xlist * dx, Ez)

    # Отобразить положение источника
    ax.plot(sourcePos * dx, 0, 'ok')

    # Отобразить положение датчика
    ax.plot(probePos * dx, 0, 'xr')

    for q in range(1, maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)
        # Полное поле/ рассеянное поле
        Hy[sourcePos - 1] -= (Sc / W0) * Signal(0, q, name_signal)
        # Граничное условие
        Ez[0] = Ez[1]
        Ez[-1] = 0
        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:-1] = Ez[1:-1] + (Hy[1:] - Hy_shift) * Sc * W0 / eps
        # Полное поле/рассеянное поле
        Ez[sourcePos] += (Sc / np.sqrt(eps * mu)) * Signal(-0.5, q + 0.5, name_signal)
        # АВС первого порядка
        Ez[0] = oldEzLeft + koeffABCLeft * (Ez[1] - Ez[0])
        oldEzLeft = Ez[1]

        # Регистрация поля в точке
        probeTimeEz[q] = Ez[probePos]
        if q % 10 == 0:
            plt.title(format(q * dt * 1e9, '.3f') + ' нc')
            line.set_ydata(Ez)
            fig.canvas.draw()
            fig.canvas.flush_events()

    # Отключить интерактивный режим по завершению анимации
    plt.ioff()
    
    # Расчёт спектра сигнала
    EzSpec = fftshift(np.abs(fft(probeTimeEz)))

    # Рассчёт шага частоты
    df = 1.0 / (maxTime * dt)
    # Рассчёт частотной сетки
    freq = np.arange(-maxTime / 2 , maxTime / 2 )*df
    # Оформляем сетку
    tlist = np.arange(0, maxTime * dt, dt) 

    # Отображение сигнала, сохраненного в датчике
    tlist = np.arange(maxTime)
    fig, ax = plt.subplots()
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel('Время, мкс')
    ax.set_ylabel('Ez, В/м')
    ax.plot(tlist * dt * 10 ** 6, probeTimeEz)
    ax.grid()
    
    # Спектр сигнала
    fig, ax = plt.subplots()
    ax.plot(freq, EzSpec / np.max(EzSpec))
    ax.set_xlabel('Частота, ГГц')
    ax.set_ylabel('|S| / |Smax|, б/р')
    ax.set_xlim(0, 10e9)
    ax.set_ylim(0, 1.1)
    ax.grid()
    plt.show()
   

