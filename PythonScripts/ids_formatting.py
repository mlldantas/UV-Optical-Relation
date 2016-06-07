import numpy as np
import os

mydata  = 'ids.csv'

plate   = np.loadtxt(mydata, dtype=str, usecols=[0], delimiter=',', skiprows=1)
mjd     = np.loadtxt(mydata, dtype=str, usecols=[1], delimiter=',', skiprows=1)
fiberid = np.loadtxt(mydata, dtype=str, usecols=[2], delimiter=',', skiprows=1)

new_plate = []
for i in range (plate.size):
    if (len(plate[i]) == 3):
        new_plate_i = os.path.join('0'+plate[i])
        new_plate.append(new_plate_i)
    elif (len(plate[i]) == 4):
        new_plate.append(plate[i])

new_fiberid = []
for j in range (fiberid.size):
    if (len(fiberid[j]) == 1):
        new_fiberid_j = os.path.join('00'+fiberid[j])
        new_fiberid.append(new_fiberid_j)
    elif (len(fiberid[j]) == 2):
        new_fiberid_j = os.path.join('0'+fiberid[j])
        new_fiberid.append(new_fiberid_j)
    elif (len(fiberid[j]) == 3):
        new_fiberid.append(fiberid[j])

formatted_ids = open('new_ids.csv', 'w')

print >> formatted_ids, '#identifiers'

for i in range (plate.size):
    line = os.path.join(new_plate[i] + '.' + mjd[i] + '.' + new_fiberid[i])
    print >> formatted_ids, line


formatted_ids.close()
