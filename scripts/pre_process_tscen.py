import sys,os
from functools import reduce

def pre_process_tscen(infilename,epoch=0):
    infile = open(infilename,'r')
    ages = list(map(float,infile.readlines())) #convert ages to floats
    ages = list(map(lambda x: x-float(epoch), ages)) #remove the epoch 0 if there isn't one
    out_str = ""
    for i,age in enumerate(reversed(ages)):
        out_str += "ts(%d)= %.3f\n"%(i+1,-age)
    final_neg_index=i+1
    for j,age in enumerate(ages):
        out_str += "ts(%d)= %.3f\n"%(j+final_neg_index+1,age)
    out_str = out_str.rstrip("\n")
    out_filename = os.path.join(os.path.dirname(infilename),reduce(lambda x,y: x+'.'+y,os.path.basename(infilename).split('.')[0:-1])) #same name minus an extension
    print("saving to %s"%out_filename)
    out_file = open(out_filename,'w+')
    out_file.write(out_str)
    out_file.close()

if __name__=="__main__":
    if len(sys.argv)<3: print("This script requires at least 2 inputs in this order: in file and epoch. Please rerun with these, aborting."); sys.exit()
    if '-h' in sys.argv: help(__name__)

    pre_process_tscen(sys.argv[1],sys.argv[2])
