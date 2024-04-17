# %%
import numpy as np
import scipy.stats as stats
import statistics
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import make_interp_spline
import scipy.integrate as integrate

# %%
# Läs in data
data = pd.read_excel("pk.xlsx")
# Data från olika kolumner


# Sammanfattar data
'''
print(data)

print(data.describe)
'''

# Exponentiell modell


def exp_decay(t, k, C0):
    return C0 * np.exp(-k * t)


fig, (ax1, ax2) = plt.subplots(2, 1)

# Tar ut datan för varje individ
individuals = data['Person'].unique()


for person in individuals:
    # Delar upp datan i person, tid samt koncentration
    person_data = data[data['Person'] == person]
    time = person_data['Time']
    concentration = person_data['Conc']
    symptom = person_data['Symptom']

    # Anpassa modellen till data för varje individ
    optimala, kovariansen = curve_fit(exp_decay, time, concentration)
    k, C0 = optimala

    # Halveringstiden
    half_life = np.log(2) / k

    # Plotta data oför varje person
    # plt.scatter(time, concentration, label='Person ' + str(person))

    spline = UnivariateSpline(time, concentration)
    x_fit = np.linspace(min(time), max(time), 1000)
    spl = make_interp_spline(time, concentration, k=2)
    y_fit = spl(x_fit)

    ax1.plot(time, exp_decay(time, k, C0), '--',
             label='Person ' + str(person-100))
    ax1.scatter(time, concentration, label='Person ' + str(person-100))
    ax1.plot(x_fit, y_fit)

    # ax1.curve_fit(time, concentration,label='Person ' + str(person-100))
    ax1.set_ylabel('Koncentrationen (mg/liter)')
    # ax1.legend()
    ax1.grid()

    ax2.plot(time, symptom, label='person' + str(person-100))
    ax2.set_ylabel('Antalet Symptom')
    ax2.set_xlabel('Tid (h)')
    # ax2.legend()

    # skriver ut parametrarna från modellen för varje individ

    print(f'Person nr {person-100}:')
    print('Lambda:', k)
    print('Halveringstid:', half_life, 'timmar')
    print()

plt.grid()
plt.show()


# %%
person_data1 = data[data['Person'] == person]
time1 = person_data['Time']
concentration1 = person_data['Conc']
symptom1 = person_data['Symptom']

optimala1, kovariansen1 = curve_fit(exp_decay, time, concentration)
k1, C01 = optimala1

half_life1 = np.log(2) / k1

spline = UnivariateSpline(time, concentration)
x_fit = np.linspace(0, 100, 1000)
spl = make_interp_spline(time, concentration, k=2)
y_fit = spl(x_fit)

plt.scatter(time, concentration)
plt.plot(x_fit, y_fit)
plt.xlabel('tid (h)')
plt.ylabel('Koncentrationen (mg/l)')


print('Uppskattade gemensam parameter')
print('lambda', k1)
print('Halveringstid:', half_life1, 'timmar')
print()

plt.grid(True)
plt.show()

D = 50*3  # total dos

integrated_C = integrate.quad(lambda x: exp_decay(x, k1, C01), 0, 10)[0]
integrated_tC = integrate.quad(lambda x: exp_decay(x, k1, C01) * x, 0, 10)[0]

half_life = np.log(2)/k1

clearence = D / integrated_C

mrt = integrated_tC/integrated_C

v_d = clearence/k1

v_ss = clearence * mrt


print("Halflife: ", half_life)
print("Clearence: ", clearence)
print("MRT: ", mrt)
print("V_d: ", v_d)
print("V_ss: ", v_ss)
# %%
