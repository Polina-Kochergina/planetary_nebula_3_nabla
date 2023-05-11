from matplotlib import pyplot as plt
import numpy as np
from scipy import integrate
import scipy.optimize as opt
import pandas as pd

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def Planet_Neb(file):

    data = pd.read_csv(f"data\{file}", delimiter=",", names=["wave", "flux"], skiprows=1)
    waves = data["wave"]
    flux = data["flux"]
    print(data.shape)
    n = data.shape[0]//700
    n += 1

    # for i in range(n):
    fig, axes = plt.subplots()
        # axes.plot(waves[i*700 : (i+1)*750], flux[i*700 : (i+1)*750], "k")
        # fig.savefig(f"img3{i}.png")

    axes.plot(waves[700 : 750], flux[700 : 750], "k")

    b = 4862.5
    print(data.sort_values("flux", ascending = False)[:20])
    l = 4861
    dl = b - l
    # in_x0 = [3870, 3972, 4110, 4340, 4690, 4856, 4957, 5003, 5878, 6550] #
    z = dl/l
    print("z =", z)

    fig, axes = plt.subplots()
    
    axes.plot(waves[700 : 900], flux[700 : 900], "k")
    

    poptHa, _ = opt.curve_fit(gauss, waves, flux, p0 = [1, 1*10**(-13), 6550, 5])
    poptHb, _ = opt.curve_fit(gauss, waves, flux, p0 = [1, 1*10**(-13), 4861, 5])
    poptS1, _ = opt.curve_fit(gauss, waves, flux, p0 = [1, 1*10**(-14), 6716, 0.1])
    poptS2, _ = opt.curve_fit(gauss, waves, flux, p0 = [1, 1*10**(-14), 6731, 0.1])

    print(poptHa)
    print(poptHb)
    print(poptS1)
    print(poptS2)


    intens_Ha = integrate.quad(gauss, 6584-30, 6584+30, args=tuple(poptHa))
    intens_Hb = integrate.quad(gauss, 4861-25, 4861+25, args=tuple(poptHb))
    intens_S1 = integrate.quad(gauss, 6716-15, 6716+15, args=tuple(poptS1))
    intens_S2 = integrate.quad(gauss, 6731-15, 6731+15, args=tuple(poptS2))
    # print(intens_Ha[0])
    # print(intens_Hb[0])
    # print(intens_S1[0])
    # print(intens_S2[0])
    # plt.plot(waves[2080:2800], flux[2080:2800])
    # plt.plot(waves[2080:2800], gauss(waves[2080:2800], *poptS1))



    # plt.plot(waves, flux)

    # plt.plot(waves, gauss(waves, *poptHa))
    # plt.plot(waves, gauss(waves, *poptHb))
    # plt.plot(waves, gauss(waves, *poptS1))
    # plt.plot(waves, gauss(waves, *poptS2))

    print(intens_Ha[0]/intens_Hb[0], "balmer")
    print(intens_S1[0]/intens_S2[0], "sera")

    return intens_Ha, intens_Hb, intens_S1, intens_S2


# intens = Planet_Neb("A2.csv")
# intens = Planet_Neb("NGC6537.csv")
intens = Planet_Neb("IC289.csv")

# T_eff < 50 000 K получилось, линий из табл нет, но есть неионизов гелий и линии водорода
# print(data)

plt.show()