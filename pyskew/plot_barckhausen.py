from pmagpy import demag_gui_utilities as dgu
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from track_processing import *
from re import findall

chron_to_analyse = (27,'c')
chrons_info = [(24,'green'),(25,'pink'),(26,'r'),(27,'c'),(28,'orange'),(29,'aqua'),(30,'grey')]
results_directory = '27_results'

ndf = open('GSFML.Barckhausen++_2013_MGR.picks.gmt','r')

#lonlat = np.array([list(map(float,line.strip('\n').split('\t'))) for line in ndf.readlines() if not line.startswith('#')])

meta_lonlat = ndf.readlines()[7:]

m = (9-3)/(-110-(-145))
b = 3-m*(-145)+1
ccz_div_line = lambda x: m*x+b #yes I picked this line by hand as it was all I could think of quit judging me. -_-

ccz_chrons = {}
gcz_chrons = {}
for line in meta_lonlat:
    if line.startswith('#'):
        chron = line.lstrip('# @DC').split('|')[0].replace('n','').replace('.1','')
        if chron not in list(chrons.keys()): chrons[chron] = []
        else: continue
    else:
        nd = list(map(float,line.split('\t')))
        if nd[1]>ccz_div_line(nd[0]): ccz_chrons[chron].append(nd)
        else: gcz_chrons[chron].append(nd)

def sort_cmp(x,y=''):
    if y=='': return 0
    sort_by_num = int(findall(r'\d+', x)[0])-int(findall(r'\d+', y)[0])
    if sort_by_num!=0: return sort_by_num
    else:
        for c1,c2 in zip(x,y):
            if (ord(c1)-ord(c2))!=0: return ord(c1)-ord(c2)
    if len(x)>len(y): return 1
    elif len(y)>len(x): return -1
    else: return 0

sort_key = dgu.cmp_to_key(sort_cmp)

chrons = ccz_chrons
m = create_basic_map(projection='merc')
m.drawparallels(np.arange(-30,50,5),labels=[1,0,0,0],linewidth=0)
m.drawmeridians(np.arange(-190,-90,5),labels=[0,0,0,1],linewidth=0)
for i,chron in enumerate(sorted(chrons.keys(),key=sort_key)):
    lonlat = np.array(chrons[chron])

    if len(lonlat)==0: continue

    (lon,lat) = m(lonlat[:,0],lonlat[:,1])

    m.plot(lon,lat,label=chron)

plt.title("Barckhausen 2013 Magnetic Anomalies")
plt.legend()

plt.show()
