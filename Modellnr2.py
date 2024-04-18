import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#%%
data = pd.read_excel("pk.xlsx")

def model(t, vec):
    A, B, F, ka, lam, mu = vec
    return (F * ka * A * (np.exp(-lam * t) - np.exp(-ka * t))/ (ka - lam) + 
            F * ka * B  * (np.exp(-mu * t) - np.exp(-ka * t))/ (ka - mu))


def minimera(vec, t, y):
    return np.mean((y - model(t, vec))**2)

results = {}
vec_guess = [16, 0.5, 0, 0.1, 2, 3]
time_vec = []
conc_vec = []

for person in data['Person'].unique(): 
    person_data = data[data['Person'] == person]
    time = person_data['Time']
    time_vec.extend(time)
    concentration = person_data['Conc']
    conc_vec.extend(concentration)
    res = minimize(minimera, vec_guess, args=(time, concentration))
    results[person] = res.x
    print(f'Person {person}: {res.x}')

#%%
time_space = np.linspace(0, 100, 1000)
params = [24.41072424, -23.71260049, 26.63106566, 0.12670922, 1.15457906, 1.18363399]
concenctration_plot = model(time_space, params)

plt.scatter(time_vec, conc_vec)
plt.plot(time_space, concenctration_plot)
plt.xlabel('Tid')
plt.ylabel('Koncentration')
plt.grid()
plt.show()
