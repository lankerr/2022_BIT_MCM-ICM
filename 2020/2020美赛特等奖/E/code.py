import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import lognorm
import time
import csv

AP = np.zeros(202)
APs = np.zeros(202)
APo = np.zeros(202)
AP1a = np.zeros(202)
APoa = np.zeros(202)
AWs = np.zeros(202)
AWo = np.zeros(202)
AW1a = np.zeros(202)
AWoa = np.zeros(202)
a = 104580.93
b = 409084324.46
c = 400055166835.44
D = [0.1, 3, 7, 2, 1, 3, 1.5, 1.5]  # plastic lifetime standard deviation
E = [0.5, 13, 35, 8, 3, 20, 5, 5]  # plastic lifetime expectation
ratio = [0.359, 0.066, 0.160, 0.044, 0.103, 0.007, 0.145, 0.115]  # plastic ratio
a1 = 0.1
ao = 0.05
C = 0
P = 0


def use_age(kind, year):
    d = D[kind] ** 2
    s = math.log((1 + d / (E[kind] ** 2))) ** 0.5
    mu = math.log(E[kind]) - 0.5 * math.log((1 + D[kind] / (E[kind] ** 2)))
    res = lognorm.cdf(year, s, scale=math.exp(mu)) - lognorm.cdf(year - 1, s, scale=math.exp(mu))
    return res


def AP_Func(t):
    # Forecast the future plastic production
    return a * (t ** 2) - b * t + c


def AP_Gene(t):
    for i in range(0, t):
        tmp = AP_Func(1950 + i)
    AP[i] = tmp
    for i in range(0, t):
        tmp = AP_s(i)
    APs[i] = tmp
    tmp = APo_a(i)
    APoa[i] = tmp
    tmp = AP_o(i)
    APo[i] = tmp
    tmp = AP1_a(i)
    AP1a[i] = tmp


def AW_Gene(t):
    for i in range(1, t):
        tmp = AW_s(i)
        AWs[i] = tmp
        tmp = AW_o(i)
        AWo[i] = tmp
        tmp = AW1_a(i)
        AW1a[i] = tmp
        tmp = AWo_a(i)
        AWoa[i] = tmp


def ratio_a1(t):
    # Alternative material ratio in single-use plastic
    return min(a1 * P_Func(t), 0.8)


def ratio_a2(t):
    # Alternative material ratio in other plastic
    return min(ao * P_Func(t), 0.5)


def AP_s(t):
    # Single-use plastic production (alternated)
    return AP[t] * ratio[0] * (1 - ratio_a1(t))


def AP1_a(t):
    # Alternative material production in packaging
    if t < 67:
        return 0
    return AP[t] * ratio[0] * ratio_a1(t)


def AP_o(t):
# Other plastic production
    return AP[t] * (1 - ratio[0]) * (1 - ratio_a2(t))


def APo_a(t):
# Alternative material production in other plastic
    if t < 67:
        return 0
    return AP[t] * (1 - ratio[0]) * ratio_a2(t)


def AW_s(t):
# Annual waste in single-use plastic
    res = 0
    if t < 67:
        res += AP[t - 1] * ratio[0] * use_age(0, 1)
        res += AP[t - 2] * ratio[0] * use_age(0, 2)
    else:
        res += AP_s(t - 1) * use_age(0, 1)
        res += AP_s(t - 2) * use_age(0, 2)
    return res


def AW_o(t):
# Annual waste in other plastic
# t = year - 1950
    res = 0
    if t < 67:
        for i in range(1, t):
            for j in range(1, 8):
                res += AP[t - i] * ratio[j] * use_age(j, i)
    else:
        for i in range(1, 67):
            for j in range(1, 8):
                res += AP[67 - i] * ratio[j] * use_age(j, i)
    for i in range(67, t):
        for j in range(1, 8):
            tmp = AP_o(t - i) * ratio[j] * use_age(j, i - 66)
            res += tmp
            return res


def AW1_a(t):
# Alternative material’s waste in packaging
    res = 0
    if t < 67:
        pass
    elif t == 67:
        res += AP1_a(t - 1) * use_age(0, 1)
    else:
        res += AP1_a(t - 1) * use_age(0, 1)
        res += AP1_a(t - 2) * use_age(0, 2)
    return res


def AWo_a(t):
# Alternative material’s waste in other plastic
    if t < 67:
        return 0
    else:
        res = 0
    for i in range(67, t):
        for j in range(1, 8):
            res += APo_a(t - i) * ratio[j] * use_age(j, i - 66)
    return res


def C_Func(t):
    if t < 66:
        return 1
    else:
        return math.exp(C * (t - 65))

def P_Func(t):
    if t < 66:
        return 1
    else:
        return math.exp(P * (t - 65))

def i_o(t):
# incineration ratio in other plastic
    year = t + 1950
    return min((0.7045 * year - 1394) / 100, (0.8 - r_o(t)), 0.5)


def i_s(t):
# incineration ratio in single-use plastic
    year = t + 1950
    return min((0.7045 * year - 1394) / 100, 0.5)


def r_o(t):
# recycling ratio in other plastic
    year = t + 1950
    return min((0.7 * year - 1391) / 100 * C_Func(t), 0.5)


def J(t):
# human’s recollection of the plastic waste
    year = t + 1950
    if year < 2013:
        return 0
    res = 500000 * (year - 2013) - 7500000
    return res


def Hp(t):
# harm
    return AWs[t] * (1 - i_s(t)) + AWo[t] * (1 - i_o(t) - r_o(t))

def Hp_origin(t):
# origin harm
    res1 = 0
    res1 += AP[t - 1] * 0.359 * use_age(0, 1)
    res1 += AP[t - 2] * 0.359 * use_age(0, 2)
    res2 = 0
    for i in range(1, t + 1):
        for j in range(1, 8):
            res2 += AP[t - i] * ratio[j] * use_age(j, i)
    year = t + 1950
    r10 = (0.7045 * year - 1394) / 100  # origin incineration ratio
    r20 = (0.7 * year - 1391) / 100  # origin recycling ratio
    r2 = min(r20, 0.5, 0.8 - r10)
    r1 = min(r10, 0.5)
    res = res1 * (1 - r1) + res2 * (1 - r1 - r2)
    return res


def max_H(Pollu):
    maxh = 0
    max_h_year = 1950
    t = 0
    for h in Pollu:
        t += 1
    if h > maxh:
        maxh = h
    max_h_year = 1950 + t
    print(str(max_h_year) + ": " + "Max H is " + format(maxh, ’.2f’))

def figure(Original=None):
    plt.figure(1)
    plt.subplot(121)
    plt.plot(Year, Origin, label=’Original’, color =’b’)
    plt.plot(Year, Pollution, label=’Alter’, color =’r’)
    plt.plot(Year, H_J, label=’H - J’, color =’g’)
    plt.title(’Pollution’)
    plt.legend()
    plt.subplot(122)
    plt.plot(Year, Accum, label=’Alter’, color =’r’)
    plt.plot(Year, Acc_origin, label=’Origin’, color =’b’)
    plt.title(’Accumulation’)
    plt.legend()
    plt.figure(3)
    plt.plot(Year, AWs[:200], label=’AW_s’)
    plt.legend()
    plt.figure(3)
    plt.subplot(121)
    plt.plot(Year, Origin, label=’Original’, color =’b’)
    plt.plot(Year, Pollution, label=’Alter’, color =’r’)
    plt.plot(Year, H_J, label=’H - J’, color =’g’)
    plt.title(’Pollution’)
    plt.legend()
    plt.subplot(122)
    plt.plot(Year, Accum, label=’Alter’, color =’r’)
    plt.plot(Year, Acc_origin, label=’Origin’, color =’b’)
    plt.title(’Accumulation’)
    plt.legend()
    plt.figure(2)
    plt.plot(Year, AWs[:200], label=’AW_s’)
    plt.legend()
    plt.show()


def write():
    with open(’full_data_1951-2150.csv’, ’w’, newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Year", "AP", "AP_o", "AP_s", "APo_a",
                 "AP1_a", "AW_o", "AW_s", "AWo_a", "AW1_a", "H", "H-J"])
    for i in T:
        writer.writerow([i + 1950, AP[i], APo[i], APs[i], APoa[i], AP1a[i], AWo[i],AWs[i], AWoa[i], AW1a[i], Pollution[i],H_J[i]])
if __name__ ==  ’__main__’:
    with open(’pc_data.csv’, ’w’, newline="") as csvfile:
        writer = csv.writer(csvfile)
    writer.writerow([’a1’, ’ao’, ’H(2120)’])
    while C < 0.020:
        C += 0.002
        P = 0
    while P < 0.03:
        P += 0.003
    AP_Gene(130)
    AW_Gene(130)
    res = Hp(100) - J(100)
    print("p: " + str(P) + ", c: " + str(C) + ", H: " +format(res, ’.2f’) + " AW_s: "+ format(AWs[100], ’.2f’) + " AW_o: " +
    format(AWo[100], ’.2f’))writer.writerow([P, C, format(res, ’.2f’)])
    AP_Gene(201)
    AW_Gene(201)
    Year = np.arange(1951, 2151)
    T = Year - 1950
    Pollution = []
    H_J = []
    Origin = []
    Accum = []
    Acc_origin = []
    for i in T:
        Origin.append(Hp_origin(i))
    H0 = Hp(i)
    Pollution.append(H0)
    HJ = Hp(i) - J(i)
    H_J.append(HJ)
    res = 0
    for item in Pollution:
        res += item
    Accum.append(res)
    res = 0
    for item in Origin:
        res += item
    Acc_origin.append(res)
    AWs_origin = []
    for i in T:
        res = AP[i - 1] * ratio[0] * use_age(0, 1) + AP[i - 2] * ratio[0] *
    use_age(0, 2)
    AWs_origin.append(res)
    max_H(Pollution)
    max_H(H_J)
    with open(’714_1951-2150.csv’, ’w’, newline="") as csvfile:
        writer = csv.writer(csvfile)
    writer.writerow(["Year", "AP", "AP_o", "AP_s", "APo_a",
                 "AP1_a", "AW_o", "AW_s", "AWo_a", "AW1_a", "H",
                 "H-J", "H_origin", "Accum", "Accum_origin"])
    for t in T:
        i = t - 1
    writer.writerow(
        [i + 1950, AP[i], APo[i], APs[i], APoa[i], AP1a[i], AWo[i],
        AWs[i], AWoa[i], AW1a[i], Pollution[i],
        H_J[i], Origin[i], Accum[i], Acc_origin[i]])
