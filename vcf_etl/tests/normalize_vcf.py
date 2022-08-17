import logging
import pandas as pd
import io
import re

def get_vcf_as_table(vcf_in):
    
    with open(vcf_in, 'r') as f:
        infoLines = [] 
        headers = []
        for line in f:
            if not line.startswith('##'):
                infoLines.append(line.strip())
            else:
                headers.append(line.strip())
    
    table = pd.read_csv(
        io.StringIO('\n'.join(infoLines)),
        sep='\t'
    )
    
    return table, headers


def remove_no_pass(data):

    acceptable_statuses = ["PASS", "LowGQX"]
    data = data[data.FILTER.isin(acceptable_statuses)]
    return data

def extractAlleleFrequencies(infoCol):
    afs_extracted = []
    for entry in infoCol:
        search_info = re.search("AF1000G=", entry)
        if search_info is None:
            afs_extracted.append(".")
        else:
            start = search_info.start()
            current = start + 8
            while current < len(entry):
                if entry[current] == ';':
                    break
                current += 1
            freqency = entry[start+8:current]
            afs_extracted.append("AF="+freqency)
    return afs_extracted

def trim_info(data):

    afs_formatted = extractAlleleFrequencies(data.INFO.tolist())
    data.INFO = afs_formatted
    return data
    
    
def trim_format(data):
   
    fmt = data.FORMAT
    sample = data.iloc[:,-1].tolist()

    formats = [x.split(":") for x in fmt]
    samples = [x.split(":") for x in sample]
    desired_fmts = ['GT', 'AD', 'DP']
    new_format_list = []
    new_sample_list = []
    for i in range(len(formats)):
        current_format = formats[i]
        current_sample = samples[i]
        new_format_entry = ""
        new_sample_entry = ""
        for j in range(len(current_format)):
            if current_format[j] in desired_fmts:
                new_format_entry += current_format[j] + ":"
                new_sample_entry += current_sample[j] + ":"
        if new_format_entry == "":
            new_format_entry = ". "
            new_sample_entry = ". "
        new_format_list.append(new_format_entry[0:-1])
        new_sample_list.append(new_sample_entry[0:-1])

    data.FORMAT = new_format_list
    data.iloc[:,-1] = new_sample_list

    return data #Do I need to adjust this to take in any number of sample columns??
                

def remove_no_alternate(data):
    data = data[data.ALT != "."]
    return data

def filterHeaders(headers):

    filtered = [head for head in headers if not head.startswith("##INFO=")]
    filtered.append('##INFO=<ID=AF,Number=A,Type=Float,Description="The allele frequency">')

    return filtered

def normalize_vcf(vcf_in: str, vcf_out: str):
    
    data, headers = get_vcf_as_table(vcf_in)

    data = remove_no_pass(data)
    data = trim_info(data)
    data = trim_format(data)
    data = remove_no_alternate(data)

    headers = filterHeaders(headers)

    newVcfString = data.to_csv(sep="\t", index=False)
    newVcfString = "\n".join(headers + [newVcfString])

    outfile = open(vcf_out, 'w')
    outfile.write(newVcfString)    
  