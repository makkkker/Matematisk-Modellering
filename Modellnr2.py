import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad
import matplotlib.pyplot as plt

#%%
data = pd.read_excel("pk.xlsx")

#Definera modell
def model(t, vec):
    A, B, F, ka, lam, mu = vec
    return (F * ka * A * (np.exp(-lam * t) - np.exp(-ka * t))/ (ka - lam) + 
            F * ka * B  * (np.exp(-mu * t) - np.exp(-ka * t))/ (ka - mu))

#Minsta kvadratmetoden
def minimera(vec, t, y):
    return np.mean((y - model(t, vec))**2)

#Lista på parametrar för varje person
results = {}

time_vec = []
conc_vec = []

#Gissar parametervärden
vec_guess = [16, 0.5, 0, 0.1, 2, 3]

#Itererar igenom varje person i datan
for person in data['Person'].unique(): 
    #Plockar ut data för varje individ
    person_data = data[data['Person'] == person]
    time = person_data['Time']
    time_vec.extend(time)
    concentration = person_data['Conc']
    conc_vec.extend(concentration)
    #Här används minsta kvadratmetoden för att hitta parametrar för varje individ
    res = minimize(minimera, vec_guess, args=(time, concentration))
    #Här läggs parametrarna till i results 
    results[person] = res.x
    #Skriver ut parametrar för varje person
    print(f'Person {person-100}: {res.x}')

#%%
#Skapar ett tidsspann
time_space = np.linspace(0, 100, 1000)
#Använder dessa parametrar (Man ska nog använda andra) för att plotta funktionen
params = [24.41072424, -23.71260049, 26.63106566, 0.12670922, 1.15457906, 1.18363399]
concentration_plot = model(time_space, params)

#Plottar datapunkter 
plt.scatter(time_vec, conc_vec)
#Plottar funktionen efter de valda parametrarna 
plt.plot(time_space, concentration_plot)
plt.xlabel('Tid')
plt.ylabel('Koncentration')
plt.grid()
plt.show()

#Arean
AUC= np.trapz(concentration_plot, time_space)
print('AUC:', AUC)

#Clearance
CL = 150 / AUC
print('Clearance:', CL)

#Mean residence time
MRT = np.trapz(time_space*concentration_plot, time_space)/AUC
print('MRT: ', MRT)

#Volym i steady state
Vss = CL*MRT
print('Vss: ', Vss)


# %%
