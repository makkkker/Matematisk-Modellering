import numpy as np
import scipy.stats as stats
import statistics
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit

#Läs in data 
data = pd.read_excel("pk.xlsx")
#Data från olika kolumner


#Sammanfattar data
'''
print(data)

print(data.describe)
'''

#Exponentiell modell
def exp_decay(t, k, C0):
    return C0 * np.exp(-k * t)


fig, (ax1, ax2) = plt.subplots(2,1)

#Tar ut datan för varje individ
individuals = data['Person'].unique()

for person in individuals:
    #Delar upp datan i person, tid samt koncentration
    person_data = data[data['Person'] == person]
    time = person_data['Time']
    concentration = person_data['Conc']
    symptom = person_data['Symptom']

    # Anpassa modellen till data för varje individ
    optimala, kovariansen = curve_fit(exp_decay, time, concentration)
    k, C0 = optimala

    #Halveringstiden 
    half_life = np.log(2) / k

    # Plotta data oför varje person
    #plt.scatter(time, concentration, label='Person ' + str(person))

    ax1.plot(time, exp_decay(time, k, C0), '--', label='Person ' + str(person-100))
    ax1.scatter(time, concentration,label='Person ' + str(person-100))
    #ax1.curve_fit(time, concentration,label='Person ' + str(person-100))
    ax1.set_ylabel('Koncentrationen (mg/liter)')
    #ax1.legend()
    ax1.grid()
    

    ax2.plot(time, symptom, label ='person' + str(person-100))
    ax2.set_ylabel('Antalet Symptom')
    ax2.set_xlabel('Tid (h)')
    #ax2.legend()




    # skriver ut parametrarna från modellen för varje individ
    '''
    print(f'Person nr {person-100}:')
    print('k:', k)
    print('Halveringstid:', half_life, 'timmar')
    print()
    '''
plt.grid()
plt.show()