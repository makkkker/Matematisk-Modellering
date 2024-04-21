import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt

#%%
data = pd.read_excel("pk.xlsx")

#Def modell
def modelnr2(t, F, A, B, ka, lam, mu):
    return (F * ka * A * (np.exp(-lam * t) - np.exp(-ka * t))/ (ka - lam) + 
            F * ka * B  * (np.exp(-mu * t) - np.exp(-ka * t))/ (ka - mu))

fittedCurve = []
vec_guess = [0.4, 5, 5, 0.5, 0.09, 0.3]

#Itererar igenom varje individ
for person in data['Person'].unique(): 
    fittedCurve.append([])

    person_data = data[data['Person'] == person]
    time = person_data['Time'].to_numpy()
    concentration = person_data['Conc'].to_numpy()
    
    #Tar maxvärdet av konc. för varje invdivid och startar därifrån 
    maxindex = concentration.argmax()
    time = time[maxindex: ]
    concentration = concentration[maxindex: ]
    time_space = np.linspace(time[0], 96, 10000)
    params1, params1_cov = curve_fit(modelnr2, time, concentration, vec_guess)

    #appendar varje inviduella funktion och plottar 
    for t in time_space: 
        fittedCurve[person-101].append(modelnr2(t,*params1))
    plt.plot(time_space, fittedCurve[person-101])
    plt.scatter(time, concentration)
    print(f'Person {person-100}: F={params1[0]} A={params1[1]} B={params1[2]} ka={params1[3]} lamda={params1[4]} mu={params1[5]}')

AUC = []
MRT = []
Clearence = []
Vss = []

#Inte klar med konstanterna 
for person in fittedCurve:
    AUC.append(np.trapz(person,time_space))
    

print(AUC + '\n')

plt.xlabel('Tid(h)')
plt.ylabel('Koncentration (mg/l)')
plt.grid()
plt.show()
# %%
#Sammma som innan med jag börjar från 0 och inte maxkoncentrationen
fittedCurve1 = []
for person in data['Person'].unique(): 
    fittedCurve1.append([])

    person_data = data[data['Person'] == person]
    time = person_data['Time'].to_numpy()
    concentration = person_data['Conc'].to_numpy()

    time_space1 = np.linspace(0, 96, 10000)
    params, params_cov = curve_fit(modelnr2, time, concentration, vec_guess)

    for t in time_space1: 
        fittedCurve1[person-101].append(modelnr2(t,*params))
    plt.plot(time_space, fittedCurve1[person-101])
    print(f'Person {person-100}: F={params[0]} A={params[1]} B={params[2]} ka={params[3]} lamda={params[4]} mu={params[5]}')

plt.xlabel('Tid(h)')
plt.ylabel('Koncentration (mg/l)')
plt.grid()
plt.show()

# %%


