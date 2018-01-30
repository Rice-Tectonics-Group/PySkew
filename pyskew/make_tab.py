import os,sys,glob
import pandas as pd

directory = sys.argv[1]

ship_cut_tracks = glob.glob(os.path.join(directory,'**','**','**','*.c?'))
aero_cut_tracks = glob.glob(os.path.join(directory,'**','**','**','**','*.c?'))
print("------NOTE: all _tab files will be skipped to prevent recursion-------")
print(ship_cut_tracks+aero_cut_tracks)

for cut_track in ship_cut_tracks:
    if '_tab.' in cut_track: continue
    df = pd.read_csv(cut_track,sep=' ',names=['dis','dec_year','mag','lat','lon'])
    df.to_csv(cut_track.split('.')[0]+'_tab.'+cut_track.split('.')[1],sep='\t',index=False)

for cut_track in aero_cut_tracks:
    if '_tab.' in cut_track: continue
    df = pd.read_csv(cut_track,sep=' ',names=['dis','lat','lon','u1','u2','u3','u4','u5','u6','u7','u8','u9'])
    df.to_csv(cut_track.split('.')[0]+'_tab.'+cut_track.split('.')[1],sep='\t',index=False)
