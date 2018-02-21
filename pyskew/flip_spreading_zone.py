import sys,os,shutil
import pandas as pd
from utilities import open_mag_file,write_mag_file_df

def flip_spreading_zone(deskew_path,spreading_zone):
    dskf = open_deskew_file(deskew_path)
    for i,row in dskf.iterrows():
        if row['sz_name']==spreading_zone:
            flip_data_file(os.path.join(row['data_dir'],row['comp_name']))

def flip_data_file(data_path):

    #backup old data file
    data_back = data_path+'.bak'
    print('backing up %s to %s'%(data_path,data_back))
    shutil.copyfile(data_path,data_back)

    #flip data file
    data_df = open_mag_file(data_path)
    print("flipping data in %s"%data_path)
    import pdb; pdb.set_trace()
    data_df = data_df.iloc[::-1]

    #write out data file
    print("overwriting %s"%data_path)
    write_mag_file_df(data_df,data_path)

if __name__=='__main__':
    if len(sys.argv)>2 and os.path.isfile(sys.argv[1]):
        deskew_path = sys.argv[1]
        spreading_zone = sys.argv[2]
        flip_spreading_zone(deskew_path,spreading_zone)
    else:
        raise ValueError("Inproper inputs to script. This script requires a deskew file to look up data files and a spreading zone for which you want the data flipped, in that order.")
