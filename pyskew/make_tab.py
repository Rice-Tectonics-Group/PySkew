import os,sys,glob
import pandas as pd

directory = sys.argv[1]

cut_tracks = glob.glob(os.path.join(directory,'**','**','**','*.c?'))
print(cut_tracks)

for cut_track in cut_tracks:
    df = pd.read_csv(cut_track,sep=' ')
    df.to_csv(cut_track.split('.')[0]+'_tab.'+cut_track.split('.')[1],sep='\t',header=False,index=False)
