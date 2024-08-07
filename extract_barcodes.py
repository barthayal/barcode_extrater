__author__ = 'espinozada'
import pandas as pd
import numpy as np
import glob
import pickle
import multiprocessing as mp
import sys
from functools import partial
 
 
'''
this function will create a dictionary
of the barcodes in a fastq file with
reads over 100.
'''
def extractBarcodes(fastq, libID, threshold):
    print ("Extracting " + str(fastq))
    file_stream = open(fastq,'r')
    nonthresh_dict = {}
    id_length = len(libID)
    reads = 0
    mapped = 0
    '''go through this'''
    for i, line in enumerate(file_stream):
        if (i % 4 == 1):
            reads += 1
            if line[:id_length] == libID:
                mapped += 1
                if line[:50] in nonthresh_dict:
                    nonthresh_dict[line[:50]] += 1
                else:
                    nonthresh_dict[line[:50]] = 1
    file_stream.close()
    readmeinfo = (fastq, mapped, reads, (100*float(mapped)/reads),  threshold, libID)
    thresh_dict = {}
    for key in nonthresh_dict:
        if nonthresh_dict[key] > threshold:
            thresh_dict[key] = nonthresh_dict[key]
    return (thresh_dict, fastq, readmeinfo)
 
 
 
if __name__ == "__main__":
 
    '''inputs by user'''
    previous_boolean = sys.argv[1] #has this been done for this before'
    old_file_name = sys.argv[2] #if so where'
    userlibID = sys.argv[3] #LID'
    barcode_length = int(sys.argv[4]) #'barcode length'
    input_thresh = int(sys.argv[5]) #'we usually use 100'
    newfile = sys.argv[6] #'name of file to create'
    directory = sys.argv[7]
     
     
    if(previous_boolean == "Y"):
        old_file = open(old_file_name,'rb')
        old_dictlist = pickle.load(old_file)
        previous_files = [x[1] for x in old_dictlist]
        fileList = []
        for file in glob.glob(directory+"/*.fastq"):
            if file not in previous_files:
                fileList.append(file)
    elif(previous_boolean == "N"):
        fileList = []
        for file in glob.glob(directory+"/*.fastq"):
                fileList.append(file)
    else:
        exit("Must choose Y or N.")
 
    dictpool = mp.Pool(processes=4)
    partial_extract = partial(extractBarcodes, libID = userlibID, threshold = input_thresh)
    new_dictlist = dictpool.map(partial_extract, fileList)
    if(previous_boolean == "Y"):
        new_dictlist = new_dictlist + old_dictlist
 
    file_dump = open(newfile,'wb')
    pickle.dump(new_dictlist, file_dump)
    file_dump.close()
