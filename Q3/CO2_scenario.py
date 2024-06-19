import shutil
import sys
import numpy as np

#the simulation year from my_most_plasim_run.sh
year = int(sys.argv[1])

#year to start the simulation from, in simulation year 001
startyear = 1850

#namelist file to edit
infile = 'radmod_namelist'

#namelist parameter to edit
itemname = 'CO2' 

#RCPscenario choose
scenario = 'RCP85'  #choose from: 'RCP85' , 'RCP6' , 'RCP45', 'RCP3' 

#open the namelist file
f=open(infile, 'r')
lines = f.readlines()
f.close()

newfile = infile + '.new'
f=open(newfile, 'w')

CO2_data = np.loadtxt('CO2_scenarios.csv',skiprows=5,delimiter=',',usecols=(0,1,2,3,4))
years = CO2_data[:,0]

#check if year is in the dataset range
if 1764 < startyear+year < 2501:
    if scenario == 'RCP3':
        CO2_value = CO2_data[np.where(years == startyear+year),1][0,0]
    if scenario == 'RCP45':
        CO2_value = CO2_data[np.where(years == startyear+year),2][0,0]
    if scenario == 'RCP6':
        CO2_value = CO2_data[np.where(years == startyear+year),3][0,0]
    if scenario == 'RCP85':
        CO2_value = CO2_data[np.where(years == startyear+year),4][0,0]
    else:
        print('Error: scenario not found, choose from RCP85 , RCP6 , RCP45, RCP3')
else:
    print('Error: Simulation year should be between 1765 and 2500.')

for line in lines:  
    if 'CO2' in line:
        oldvalue = float(line.split('=')[1])
        newline = ' %s = %s\n'% (itemname, str(CO2_value))
        
    else:
        newline=line
    f.write(newline)

f.close()
print('Changed the parameter '+itemname+' in '+infile+' into '+str(CO2_value))
shutil.copy(newfile, infile)
