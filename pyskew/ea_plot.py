# Plot effective remanent inclination against latitude
#
# Created 10 May 2019 by Daniel Woodworth.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from geographiclib.geodesic import Geodesic

max_file = pd.read_csv('/home/dtw6/Code/PySkew/24r/data/C24r.max', sep=',', skiprows=3)

isochron = pd.read_csv('/home/dtw6/Code/24r_isochron.txt', sep='\t')

profiles = max_file.iloc[::2]
data = max_file.iloc[1::2]

prf = profiles.iloc[:-1,0]
e_a = data.iloc[:-1,0].astype(float)
std = data.iloc[:-1,1].astype(float)
lat = data.iloc[:-1,2].astype(float)
lon = data.iloc[:-1,3].astype(float)
str = data.iloc[:-1,4].astype(float)

fig = plt.figure(figsize=(4,9),dpi=200)
ax = fig.add_subplot(111)

col = []

# Plot real data
real_data = False
if real_data:
	for i,row in prf.iteritems():
		if "d.lp" in row:
			ax.scatter(e_a[i+1], lat[i+1], color='blue', s=8)
		else:
			ax.scatter(e_a[i+1], lat[i+1], color='red', s=8)

# Calculate predicted data
plat = 80.7
plon = 1.6
d2r = 3.14/180
strk = 168
er = []
for i,row in isochron.iterrows():
	gdsc = Geodesic.WGS84.Inverse(row['Lat'], row['Lon'], plat, plon)
	decl = gdsc['azi1']
	
	er.append(np.absolute(np.arctan2(np.tan((90 - row['Lat'])*d2r), np.sin((strk - decl)*d2r)))*80)

ax.plot(er, isochron['Lat'])

ax.plot([0, 0], [lat.max(), lat.min()], 'k--', linewidth=0.5)
#ax.axis('equal')

plt.show()
