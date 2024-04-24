# %%
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt
from operator import truediv
from scipy.stats import linregress
import statsmodels.api as sm
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# %%
# läser data
data = pd.read_excel("pk.xlsx")
data_withoutperson = data.drop(columns=['Person'])
print(data_withoutperson.describe())

# Def model


def model(t, FA, FB, ka, lam, mu):
    return (FA * ka * (np.exp(-lam * t) - np.exp(-ka * t)) / (ka - lam) +
            FB * ka * (np.exp(-mu * t) - np.exp(-ka * t)) / (ka - mu))


def linear(t, a, b):
    return a*t + b


# Massa listor
fitted_curves = []
half_life_vec = []
terminal_velocity = []
time_vector = []
conc_vector = []
params_list = []
vec_guess = [4, 4, 0.5, 0.09, 0.3]
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Itererar igenom varje individ
for person in data['Person'].unique():
    fitted_curves.append([])

    # läser individuell data
    person_data = data[data['Person'] == person]
    time = person_data['Time'].to_numpy()
    concentration = person_data['Conc'].to_numpy()

    # konkatinerar data
    time_vector.extend(time)
    conc_vector.extend(concentration)

    # Tar maxvärdet av konc. för varje invdivid och startar därifrån
    maxindex = concentration.argmax()
    time = time[maxindex:]
    concentration = concentration[maxindex:]

    # tar def mängd samt parametervärde
    time_space = np.linspace(time[0], 96, 10000)
    params, params_cov = curve_fit(model, time, concentration, vec_guess)
    FA, FB, ka, lam, mu = params

    # Beräknar halveringstid och eliminationshastighet
    terminal_velocity.append(lam)
    half_life = np.log(2)/lam
    half_life_vec.append(half_life)
    params_list.append(params)

    # appendar varje inviduella funktion och plottar
    for t in time_space:
        fitted_curves[person-101].append(model(t, *params))

    # Plottar och printar varje värde
    plt.plot(time_space, fitted_curves[person-101],
             label='Person ' + str(person-100))
    plt.scatter(time, concentration)
    print(
        f'Person {person-100}: FA= {params[0]} FB= {params[1]} ka= {params[2]} lamda= {params[3]} mu= {params[4]}\n')


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Beräknar konstanterna för varje individ
AUC = []
MRT = []
Clearence = []
Vss = []
Kel = []

# Beräkning
for person in fitted_curves:
    AUC.append(np.trapz(person, time_space))
    Clearence.append(150)
    MRT.append(np.trapz(time_space*person, time_space))

Clearence = list(map(truediv, Clearence, AUC))
MRT = list(map(truediv, MRT, AUC))
Vss = np.multiply(Clearence, MRT)

# printar
for i in range(10):
    print(
        f'Person nr {i+1}: Clearence = {Clearence[i]} AUC = {AUC[i]} MRT = {MRT[i]} Vss = {Vss[i]} Halveringstid = {half_life_vec[i]} Kel: {terminal_velocity[i]} \n')

# skapar tabell och plottar
table_data = []
for i in range(10):
    table_data.append(["Person nr. " + str(i+1), round(Clearence[i], 4),
                      round(AUC[i], 4), round(MRT[i], 4), round(Vss[i], 4)])

plt.xlabel('Tid (h)')
plt.ylabel('Koncentration (mg/l)')
plt.legend()
plt.grid()

table = plt.table(cellText=table_data, loc='top', cellLoc='center', colLabels=[
                  'Person', 'Clearence', 'AUC', 'MRT', 'Vss', 'Half-life'], bbox=[1, 0, 1.1, 1])
plt.show()
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# %%
# Skriv ut min, max, mean och std för alla parametrar


Results = [half_life_vec, AUC, MRT, Clearence, Vss]
Results_names = ["Halflife", "AUC", "MRT", "Clearence", "Vss"]

for i in range(5):
    print(
        f'{Results_names[i]}:   Min: {min(Results[i])}  Max: {max(Results[i])} Mean: {np.mean(Results[i])}  STD: {np.std(Results[i])}')


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# %%
# Sammma som innan med jag börjar från 0 och inte maxkoncentrationen
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
        fitted_curves0[person-101].append(model(t, *params))
    plt.plot(time_space, fitted_curves0[person-101])

    print(
        f'Person {person-100}: FA={params[0]} FB={params[1]} ka={params[2]} lambda={params[3]} mu={params[4]}')

plt.xlabel('Tid(h)')
plt.ylabel('Koncentration (mg/l)')
plt.grid()
plt.show()


# %%
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Kombinerar data
time_vector = np.array(time_vector)
conc_vector = np.array(conc_vector)
time_space = np.linspace(0, 100, 10000)

# Får ut gemensamma parametrar
comb_par, comb_cov = curve_fit(model, time_vector, conc_vector, vec_guess)
fitted_curve = model(time_space, *comb_par)

# Plottar
plt.scatter(time_vector, conc_vector)
plt.plot(time_space, fitted_curve, color='red')

# Beräknar medelvärde samt standardavvikelsen av parametrarna
average_params = np.mean(params_list, axis=0)
params_std = np.std(params_list, axis=0)

# skriver ut gemensamma
table_data = []
for i, param_name in enumerate(['FA', 'FB', 'ka', 'lam', 'mu']):
    table_data.append(
        [param_name, round(average_params[i], 4), round(params_std[i], 4)])

# plottar
table = plt.table(cellText=table_data, loc='top', cellLoc='center', colLabels=[
                  'Parameter', 'Genomsnitt', 'Standardavvikelse'], bbox=[1, 0, 1.1, 1])
plt.xlabel('Time (h)')
plt.ylabel('Concentration (mg/l)')
plt.legend()
plt.grid()
plt.show()

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Konstanterna gemensamt
AUC_comb = np.trapz(fitted_curve, time_space)
MRT_comb = np.trapz(time_space * fitted_curve, time_space) / AUC_comb
Clearance_comb = 150/AUC_comb
Vss_comb = Clearance_comb*MRT_comb
half_life_comb = np.log(2)/comb_par[3]
elimination = comb_par[3]

print(f'FA={comb_par[0]} FB={comb_par[1]} ka={comb_par[2]} lambda(eliminationshastigen)={comb_par[3]} mu={comb_par[4]} halvering= {half_life_comb}')

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# %%
# Simulerar en given dosering (hur mycket, hur ofta, hur många samt parametrar från model)


def simulation(dose_amount, dose_interval, num_doses, model_params):
    # Skapar Df och Vf
    time_points = np.linspace(0, dose_interval * num_doses, 1000)
    concentrations = np.zeros_like(time_points)

    # För varje dos
    for i in range(num_doses):
        dose_time = i * dose_interval
        # Tid efter dos
        time_after_dose = time_points - dose_time
        # Ser till att koncentrationen håller sig vid ett visst värde innan till dos
        time_after_dose[time_after_dose < 0] = 0
        # Beräknar koncentrationen efter dosering
        dose_concentration = model(time_after_dose, *model_params)
        # Totala concentration av dos
        concentrations += dose_concentration * dose_amount
    return time_points, concentrations


# Gissning
model_params = [0.7, 4.7, 0.5, 0.02, 0.2]
# mg
dose_amount = 50
# Timmar
dose_interval = 24
# Antal doser
num_doses = 3
# Får ut värdena
time_points, concentrations = simulation(
    dose_amount, dose_interval, num_doses, model_params)

# Plottar
plt.plot(time_points, concentrations)
plt.xlabel('Tid (h)')
plt.ylabel('Koncentration (mg/l)')
plt.grid(True)
plt.show()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# %%
individual_conc = data['Conc']
individual_symp = data['Symptom']

concentration_vector = []
symp_vector = []

symp_vector.extend(individual_symp)
concentration_vector.extend(individual_conc)

optimala_symp, kovariansen_symp = curve_fit(
    linear, symp_vector, concentration_vector)
a, b = optimala_symp

conc_lin = np.linspace(0, max(concentration_vector), 1000)


plt.scatter(individual_conc, individual_symp, color='red')
plt.plot(linear(conc_lin, a, b), conc_lin)
plt.xlim(0, max(concentration_vector))
plt.ylim(0, max(symp_vector)+0.1)
plt.xlabel('Koncentration (mg/l)')
plt.ylabel('Symptomgrad')

plt.show()

# %%
