# %%
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt
from operator import truediv
import scipy.integrate as integrate

# %%
data = pd.read_excel("pk.xlsx")

# Def modell


def modelnr2(t, FA, FB, ka, lam, mu):
    return (FA * ka * (np.exp(-lam * t) - np.exp(-ka * t)) / (ka - lam) +
            FB * ka * (np.exp(-mu * t) - np.exp(-ka * t)) / (ka - mu))


fittedCurve = []
vec_guess = [2, 2, 0.5, 0.09, 0.3]
individual_params = []

# Itererar igenom varje individ
for person in data['Person'].unique():
    fittedCurve.append([])

    person_data = data[data['Person'] == person]
    time = person_data['Time'].to_numpy()
    concentration = person_data['Conc'].to_numpy()

    # Tar maxvärdet av konc. för varje invdivid och startar därifrån
    maxindex = concentration.argmax()
    time = time[maxindex:]
    concentration = concentration[maxindex:]
    time_space = np.linspace(time[0], 96, 10000)
    params1, params1_cov = curve_fit(modelnr2, time, concentration, vec_guess)

    # Lägg till i lista med individuella parametrar
    individual_params.append(params1)

    # appendar varje inviduella funktion och plottar
    for t in time_space:
        fittedCurve[person-101].append(modelnr2(t, *params1))
    plt.plot(time_space, fittedCurve[person-101],
             label='Person ' + str(person-100))
    plt.scatter(time, concentration)
    # print(
    #    f'Person {person-100}: FA = {params1[0]} FB = {params1[1]} ka = {params1[2]} lamda = {params1[3]} mu = {params1[4]} \n')


Halflife = []
AUC = []
MRT = []
Clearence = []
Vss = []

# Beräkning av konstanter
for person in individual_params:
    Halflife.append(person.max())
    integrated_C = integrate.quad(lambda x: modelnr2(x, *person), 0, 1000)[0]
    integrated_tC = integrate.quad(
        lambda x: modelnr2(x, *person) * x, 0, 1000)[0]
    Clearence.append(150/integrated_C)
    AUC.append(integrated_C)
    MRT.append(integrated_tC/integrated_C)
    Vss.append(np.multiply((150/integrated_C), (integrated_tC/integrated_C)))


'''
# Inte klar med konstanterna
for person in fittedCurve:
    AUC.append(np.trapz(person, time_space))
    Clearence.append(150)
    MRT.append(np.trapz(time_space*person, time_space))

Clearence = list(map(truediv, Clearence, AUC))
Vss = np.multiply(Clearence, MRT)
'''

for i in range(10):
    print(
        f'Person nr {i+1}: Clearence = {Clearence[i]} AUC = {AUC[i]} MRT = {MRT[i]} Vss = {Vss[i]} Halflife= {Halflife[i]}\n')

Results = [Halflife, AUC, MRT, Clearence, Vss]
Results_names = ["Halflife", "AUC", "MRT", "Clearence", "Vss"]

for i in range(5):
    print(
        f'{Results_names[i]}:   Min: {min(Results[i])}  Max: {max(Results[i])} Mean: {np.mean(Results[i])}  STD: {np.std(Results[i])}')

table_data = []
colors = []

for i in range(10):
    table_data.append(["Person nr. " + str(i+1), round(Clearence[i], 4),
                      round(AUC[i], 4), round(MRT[i], 4), round(Vss[i], 4)])
    person_color = plt.gca().get_lines()[i].get_color()
    colors.append(person_color)


plt.xlabel('Tid(h)')
plt.ylabel('Koncentration (mg/l)')
plt.legend()
plt.grid()

table = plt.table(cellText=table_data, loc='top', cellLoc='center', colLabels=[
                  'Person', 'Clearence', 'AUC', 'MRT', 'Vss'], bbox=[1, 0, 1.1, 1])
'''
for i, cell in enumerate(table.get_celld().values()): 
    if i % 5 == 0:
        continue
    cell.set_facecolor(colors[(i)//5])
'''
plt.show()
# %%
# Sammma som innan med jag börjar från 0 och inte maxkoncentrationen
fittedCurve1 = []
for person in data['Person'].unique():
    fittedCurve1.append([])

    person_data = data[data['Person'] == person]
    time = person_data['Time'].to_numpy()
    concentration = person_data['Conc'].to_numpy()

    time_space1 = np.linspace(0, 96, 10000)
    params, params_cov = curve_fit(modelnr2, time, concentration, vec_guess)

    for t in time_space1:
        fittedCurve1[person-101].append(modelnr2(t, *params))
    plt.plot(time_space, fittedCurve1[person-101])
    # print(
    #    f'Person {person-100}: FA={params[0]} FB={params[1]} ka={params[2]} lambda={params[3]} mu={params[4]}')

plt.xlabel('Tid(h)')
plt.ylabel('Koncentration (mg/l)')
plt.grid()
plt.show()

# %%
