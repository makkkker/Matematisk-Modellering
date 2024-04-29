import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt
from operator import truediv
from scipy.stats import linregress
from scipy.linalg import norm

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#%%
#läser data
data = pd.read_excel("pk.xlsx")
data_withoutperson = data.drop(columns=['Person'])
print(data_withoutperson.describe())

#Def model
def model(t, FA, FB, ka, lam, mu):
    return (FA * ka  * (np.exp(-lam * t) - np.exp(-ka * t))/ (ka - lam) + 
            FB * ka * (np.exp(-mu * t) - np.exp(-ka * t))/ (ka - mu))


#Massa listor 
fitted_curves = []
half_life_vec = []
terminal_velocity = []
time_vector = []
conc_vector = []
params_list = []
params_noF = []
vec_guess = [4, 4, 0.5, 0.09, 0.3]
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Itererar igenom varje individ
for person in data['Person'].unique(): 
    fitted_curves.append([])

    #läser individuell data
    person_data = data[data['Person'] == person]
    time = person_data['Time'].to_numpy()
    concentration = person_data['Conc'].to_numpy()

    #konkatinerar data
    time_vector.extend(time)
    conc_vector.extend(concentration)
    
    #Tar maxvärdet av konc. för varje invdivid och startar därifrån 
    maxindex = concentration.argmax()
    time = time[maxindex: ]
    concentration = concentration[maxindex: ]    

    #tar def mängd samt parametervärde
    time_space = np.linspace(time[0], 96, 10000)
    params, params_cov = curve_fit(model, time, concentration, vec_guess)
    FA, FB, ka, lam, mu = params

    #Beräknar halveringstid och eliminationshastighet
    terminal_velocity.append(lam)
    half_life = np.log(2)/lam
    half_life_vec.append(half_life)
    params_list.append(params)

    #appendar varje inviduella funktion och plottar 
    for t in time_space: 
        fitted_curves[person-101].append(model(t,*params))

    #Plottar och printar varje värde 
    plt.plot(time_space, fitted_curves[person-101], label='Person ' + str(person-100))
    plt.scatter(time, concentration)
    print(f'Person {person-100}: FA= {params[0]} FB= {params[1]} ka= {params[2]} lamda= {params[3]} mu= {params[4]}\n')


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Beräknar konstanterna för varje individ
AUC = []
MRT = []
Clearence = []
Vss = []
Kel = []

#Beräkning
for person in fitted_curves:
    AUC.append(np.trapz(person,time_space))
    Clearence.append(150)
    MRT.append(np.trapz(time_space*person,time_space))

Clearence = list(map(truediv, Clearence, AUC))
MRT = list(map(truediv,MRT, AUC))
Vss = np.multiply(Clearence,MRT)

#printar
for i in range(10):
   print(f'Person nr {i+1}: Clearence = {Clearence[i]} AUC = {AUC[i]} MRT = {MRT[i]} Vss = {Vss[i]} Halveringstid = {half_life_vec[i]} Kel: {terminal_velocity[i]} \n')

#skapar tabell och plottar
table_data = []
for i in range(10):
    table_data.append(["Person nr. " + str(i+1), round(Clearence[i],4), round(AUC[i],4), round(MRT[i],4), round(Vss[i],4)])

plt.xlabel('Tid(h)')
plt.ylabel('Koncentration (mg/l)')
plt.legend()
plt.grid()

table = plt.table(cellText=table_data,loc='top', cellLoc='center', colLabels=['Person', 'Clearence', 'AUC', 'MRT', 'Vss', 'Half-life'], bbox=[1,0,1.1,1])
plt.show()
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


# %%
#Sammma som innan med jag börjar från 0 och inte maxkoncentrationen
half_life_vec0 = []
fitted_curves0 = []
for person in data['Person'].unique(): 
    fitted_curves0.append([])

    person_data = data[data['Person'] == person]
    time = person_data['Time'].to_numpy()
    concentration = person_data['Conc'].to_numpy()

    time_space1 = np.linspace(0, 96, 10000)
    params, params_cov = curve_fit(model, time, concentration, vec_guess)

    lam = params[3]
    half_life_vec0.append(np.log(2)/lam)

    for t in time_space1: 
        fitted_curves0[person-101].append(model(t,*params))
    plt.plot(time_space, fitted_curves0[person-101])

    print(f'Person {person-100}: FA={params[0]} FB={params[1]} ka={params[2]} lambda={params[3]} mu={params[4]}')

plt.xlabel('Tid(h)')
plt.ylabel('Koncentration (mg/l)')
plt.grid()
plt.show()


# %%
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Kombinerar data 
time_vector = np.array(time_vector)
conc_vector = np.array(conc_vector)
time_space = np.linspace(0,100,10000)

#Får ut gemensamma parametrar 
comb_par, comb_cov = curve_fit(model, time_vector, conc_vector, vec_guess)
fitted_curve = model(time_space, *comb_par)

#Plottar
plt.scatter(time_vector, conc_vector)
plt.plot(time_space, fitted_curve, color='red')

#Beräknar medelvärde samt standardavvikelsen av parametrarna 
average_params = np.mean(params_list, axis=0)
params_std = np.std(params_list, axis=0)

#skriver ut gemensamma 
table_data = []
for i, param_name in enumerate(['FA', 'FB', 'ka', 'lam', 'mu']):
    table_data.append([param_name, round(average_params[i], 4), round(params_std[i], 4)])

#plottar
table = plt.table(cellText=table_data,loc='top', cellLoc='center', colLabels=['Parameter', 'Genomsnitt', 'Standardavvikelse'], bbox=[1, 0, 1.1, 1])
plt.xlabel('Time (h)')
plt.ylabel('Concentration (mg/l)')
plt.legend()
plt.grid()
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Konstanterna gemensamt
AUC_comb = np.trapz(fitted_curve, time_space)
MRT_comb = np.trapz(time_space * fitted_curve, time_space) / AUC_comb
Clearance_comb = 150/AUC_comb
Vss_comb = Clearance_comb*MRT_comb
half_life_comb = np.log(2)/comb_par[3]
elimination = comb_par[3]

print(f'FA={comb_par[0]} FB={comb_par[1]} ka={comb_par[2]} lambda(eliminationshastigen)={comb_par[3]} mu={comb_par[4]} halvering= {half_life_comb}')

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# %%
#Simulerar en given dosering( hur mycket, hur ofta, hur många samt parametrar från model)
def simulation(dose, dose_interval, params_list, time_points):
    concentrations = np.zeros_like(time_points)  
    for params in params_list:
        FA, FB, ka, lam, mu = params
        for i in range(96//dose_interval):
            dose_time = i * dose_interval
            time_after_dose = time_points - dose_time
            time_after_dose[time_after_dose < 0] = 0
            dose_concentration = model(time_after_dose, FA, FB, ka, lam, mu)*(dose/150)
            concentrations += dose_concentration  
    return time_points, concentrations


if __name__ == "__main__":
    time_points = np.linspace(0, 96, 1000)
    doses_to_simulate = [5, 10, 15, 20, 40] 
    plt.figure(figsize=(10, 6))
    for dose in doses_to_simulate:
        _, concentrations = simulation(dose, dose_interval=12, params_list=params_list, time_points=time_points)
        plt.plot(time_points, concentrations, label=f'Dos: {dose} mg')

    plt.xlabel('Tid (h)')
    plt.ylabel('Koncentration (mg/L)')
    plt.title('Koncetration av profylax med dosering varje 12 timmar')
    plt.legend()
    plt.grid()
    plt.show()
    #----------------------------------------------------------------------------------------------------------------------------------------------------

# %%
def simulation(dose, dose_interval, params_list, time_points):
    concentrations = np.zeros_like(time_points)  
    for params in params_list:
        FA, FB, ka, lam, mu = params
        for i in range(96//dose_interval):
            dose_time = i * dose_interval
            time_after_dose = time_points - dose_time
            time_after_dose[time_after_dose < 0] = 0
            dose_concentration = model(time_after_dose, FA, FB, ka, lam, mu)*(dose/150)
            concentrations += dose_concentration  
    return time_points, concentrations

if __name__ == "__main__":
    time_points = np.linspace(0, 96, 1000)
    doses_to_simulate = [5, 10, 15, 20, 40] 
    plt.figure(figsize=(10, 6))
    for dose in doses_to_simulate:
        _, concentrations = simulation(dose, dose_interval=12, params_list=params_list, time_points=time_points)
        plt.plot(time_points, concentrations, label=f'Dos: {dose} mg')

    plt.xlabel('Tid (h)')
    plt.ylabel('Koncentration (mg/L)')
    plt.title('Koncetration av profylax med dosering varje 12 timmar')
    plt.legend()
    plt.grid()
    plt.show()

#%% 
def simulation2(dose1, dose2, dose_interval, params_list, time_points):
    concentrations = np.zeros_like(time_points)  
    for params in params_list:
        FA, FB, ka, lam, mu = params
        for i in range(96//dose_interval):
            dose_time = i * dose_interval
            time_after_dose = time_points - dose_time
            time_after_dose[time_after_dose < 0] = 0
            if i ==0: 
                dose = dose1
            else: 
                dose = dose2
            dose_concentration = model(time_after_dose, FA, FB, ka, lam, mu)*(dose/150)
            concentrations += dose_concentration  
    return time_points, concentrations

if __name__ == "__main__":
    time_points = np.linspace(0, 96, 1000)
    dose1 = [5, 10, 20, 40]
    dose2 = 11
    plt.figure(figsize=(10, 6))
    for d in dose1:
        _, concentrations = simulation2(d, dose2, dose_interval=12, params_list=params_list, time_points=time_points)
        plt.plot(time_points, concentrations, label=f'Dos1: {d} mg, dos2: {dose2} mg')

    plt.xlabel('Tid (h)')
    plt.ylabel('Koncentration (mg/L)')
    plt.title('Dosering var 12:te timme för en initialdos samt en koninuerlig dosering efter')
    plt.legend()
    plt.grid()
    plt.show()
# %%
