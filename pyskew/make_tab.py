import os,sys,glob
import pandas as pd
import numpy as np
from .utilities import open_mag_file,convert_to_180_180,check_dir

def make_tabs(directory):
    ship_cut_tracks = glob.glob(os.path.join(directory,'**','**','*.lp'))
    aero_cut_tracks = glob.glob(os.path.join(directory,'**','**','**','*.DAT'))
    print("------NOTE: all _tab files will be skipped to prevent recursion-------")
    print(ship_cut_tracks+aero_cut_tracks)

    for cut_track in (ship_cut_tracks+aero_cut_tracks):
       if '_tab.' in cut_track: continue
       df = open_mag_file(cut_track)
       df.to_csv(cut_track.split('.')[0]+'_tab.'+cut_track.split('.')[1],sep='\t',index=False)

    # for cut_track in aero_cut_tracks:
    #     if '_tab.' in cut_track: continue
    #     df = pd.read_csv(cut_track,sep=' ',names=['dis','lat','lon','u1','u2','u3','u4','u5','u6','u7','u8','u9'])
    #     df.to_csv(cut_track.split('.')[0]+'_tab.'+cut_track.split('.')[1],sep='\t',index=False)

def make_cande_tab(cande_path,out_dir,coord_0_360=True):
    check_dir(out_dir)
    chron_path,n = cande_path,1
    if not os.path.isfile(chron_path): print("no file %s"%(chron_path)); return
    fchron = open(chron_path,'r')
    lines = fchron.readlines()
    fchron.close()
    entries=[[],[]]
    for line in lines[1:]:
        entry = line.split()
        if entry[0]=='>':
            if len(entries[0])<2 or len(entries[1])<2: entries=[[],[]]; continue
            lats = entries[1]
            lons = entries[0]
            df = pd.DataFrame(np.array(entries).T,columns=['lon','lat'])
            df['RGB'] = ['0,255,255' for i in range(len(entries[0]))]
            df.to_csv(os.path.join(out_dir,"sz_%d.tsv"%n),sep='\t',index=False)
            n += 1
            entries=[[],[]]
        else:
            if coord_0_360:
                entries[0].append(float(entry[0])); entries[1].append(float(entry[1]))
            else:
                entries[0].append(convert_to_180_180(entry[0])); entries[1].append(float(entry[1]))

if __name__=="__main__":
    directory = sys.argv[1]
    make_tabs(directory)
