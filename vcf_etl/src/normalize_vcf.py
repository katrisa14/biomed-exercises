import logging
import pandas as pd
import io
import re
import itertools


## Turns the vcf into a pandas dataframe for ease of operations while preserving the meta-info lines
def get_vcf_as_table(vcf_in):
    
    #adding meta-info lines to one list and data lines to another
    with open(vcf_in, 'r') as f:
        dataLines = [] 
        infoLines = []
        for line in f:
            if not line.startswith('##'):
                dataLines.append(line.strip())
            else:
                infoLines.append(line.strip())
    
    #creating a table from the data lines
    table = pd.read_csv(
        io.StringIO('\n'.join(dataLines)),
        sep='\t'
    )
    
    return table, infoLines


## Removes records from the data table that don't have a FILTER value of "PASS" or "LowGQX"
def remove_no_pass(data):
    
    acceptable_statuses = ["PASS", "LowGQX"]
    data = data[data.FILTER.isin(acceptable_statuses)]
    return data


## Takes a list of strings, representing entries in the INFO column, and returns a 
## list of strings containing only the allele frequency data. Also changes the label 
## to be "AF" instead of "AF1000G"
def extractAlleleFrequencies(entry):

    search_info = re.search("AF1000G=", entry) #finding allele frequency in the string

    if search_info is None:
        clean_entry = "." #if not found, append vcf missing data symbol (.)

    else:
        start = search_info.start()
        current = start + 8 #adding eight to get to the end of the label searched above (aka start of frequency)
        while current < len(entry):   #If it is found, traverse the string until the end of the AF frequency
            if entry[current] == ';': # ";" indicates the end of a piece of info, so we stop there
                break
            current += 1
        freqency = entry[start + 8 : current]
        clean_entry = "AF=" + freqency #Changes the label to AF from AF1000G

    return clean_entry


## Cleans the info column so only allele frequency values are included
def trim_info(data):

    cleaned_column = []
    for entry in data.INFO.tolist():
        cleaned_column.append(extractAlleleFrequencies(entry))
    data.INFO = cleaned_column
    return data


## Takes a list of labels from the format column, and a list of corresponding information
## from the sample column, extracts only the pieces of information wanted, and adds them all
## to a string separated by colons, as they appear in a vcf.
def extractDesiredSampleInfo(formats, sampleDxs):

    new_format_entry = ""
    new_sample_entry = ""

    desired_fmts = ['GT', 'AD', 'DP']

    for j in range(len(formats)):
        if formats[j] in desired_fmts:
            new_format_entry += formats[j] + ":" #Only keeps the format labels from the approved list
            new_sample_entry += sampleDxs[j] + ":" # and it's corresponding sample descriptions 
    
    if new_format_entry == "":
        new_format_entry = "." #If no approved format labels are found, appends 
        new_sample_entry = "."
    else:
        new_format_entry = new_format_entry[0:-1]  #Remove trailing colon (:)
        new_sample_entry = new_sample_entry[0:-1]
    
    return new_format_entry, new_sample_entry


## Removes all sample descriptions that aren't of the type GT, AD, or DP from the Sample 
## (last) column, and its labels from the FORMAT column. 
def trim_format(data):
   
    fmt = data.FORMAT        #grabbing the format column
    sample = data.iloc[:,-1] #grabbing the last (sample) column

    #creating nested lists of each sample description for each entry 
    formats = [x.split(":") for x in fmt]
    samples = [x.split(":") for x in sample]

    new_format_list = []
    new_sample_list = []

    #traversing each entry 
    for (current_format, current_sample) in zip(formats, samples):

        new_format, new_sample = extractDesiredSampleInfo(current_format, current_sample)
        
        new_format_list.append(new_format)
        new_sample_list.append(new_sample)

    data.FORMAT = new_format_list
    data.iloc[:,-1] = new_sample_list

    return data 
                

## Removes records from the data table that have a missing (.) value for the ALT column
def remove_no_alternate(data):
    data = data[data.ALT != "."]
    return data


## Removes all meta info lines describing information now removed from the INFO column
def filter_meta_info(meta_info):

    filtered = [line for line in meta_info if not line.startswith("##INFO=")]
    filtered.append('##INFO=<ID=AF,Number=A,Type=Float,Description="The allele frequency">')

    return filtered


## Calls all other functions and outputs resulting data into a new, 'normalized' vcf
def normalize_vcf(vcf_in: str, vcf_out: str):
    
    data, meta_info = get_vcf_as_table(vcf_in)

    data = remove_no_pass(data)
    data = trim_info(data)
    data = trim_format(data)
    data = remove_no_alternate(data)

    meta_info = filter_meta_info(meta_info)

    newVcfString = data.to_csv(sep="\t", index=False) #Converts data table back to a string
    newVcfString = "\n".join(meta_info + [newVcfString]) #Adds meta info lines back on top

    outfile = open(vcf_out, 'w') 
    outfile.write(newVcfString)  #and write the whole thing out to a new vcf  
  