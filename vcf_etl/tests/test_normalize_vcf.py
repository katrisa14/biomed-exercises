from unittest import result
from src.normalize_vcf import normalize_vcf
from src.normalize_vcf import remove_no_pass
from src.normalize_vcf import trim_info
from src.normalize_vcf import trim_format
from src.normalize_vcf import remove_no_alternate
from src.normalize_vcf import filter_meta_info
from src.normalize_vcf import get_vcf_as_table
import filecmp
import pandas as pd


def test_get_vcf_as_table():
    
    result_table, result_mi = get_vcf_as_table("data/source_test1.vcf")

    table_after = pd.read_csv("data/table_test1.tsv", sep="\t")
    mi_after = ['##fileformat=VCFv4.0',
        '##fileDate=20090805',
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FILTER=<ID=q10,Description="Quality below 10">',
        '##FILTER=<ID=s50,Description="Less than 50% of samples have data">',
        '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">']
    
    assert result_table.equals(table_after)
    assert result_mi == mi_after


def test_remove_no_pass():
    data = ["PASS", "LowGQX", "pass", "PASS", "Whatever", "BLue", "Red", "42"]
    df = pd.DataFrame(data, columns=["FILTER"])
    result = remove_no_pass(df)

    assert len(result.index) == 3

def test_trim_info():
    data_before = ["AF1000G=0.134;Filler=foo;Extra=Bar", "High=100;low=75;AFT=0.28;AF1000G=0.17;time=12:23", "src=refseq;lorem=ipsum/dolor|sitamit;", " Sed=aliquet/urna|vitae:dui;AF1000G=0.89;vestibulum=mattis.Ut", "a=turpis;Praesent=sagittis,AF1000G=0.9", "."]
    df_before = pd.DataFrame(data_before, columns=["INFO"])
    result = trim_info(df_before)

    data_after = ["AF=0.134", "AF=0.17", ".", "AF=0.89", "AF=0.9", "."]
    df_after = pd.DataFrame(data_after, columns=["INFO"])

    assert result.equals(df_after)

def test_trim_format():
    data_before = [["GT:AD:DP", "0:78:8"],["GT:DP:TY", "90:0/9:89"], ["AD:GT:DP", "hi:hey:howdy"], ["DP:SR:TV:QR", "89:90:i:j"], ["LOL:GG:RN", "b:b:b"],["BB:GT:HG:RN:ILY:SYA", "56:88:t:8:8:8"], ["GH:TP:AD:YMCA:PCP:KW:DP","23:22|67:0/0:0.01:hi:3#8:90" ]]
    df_before = pd.DataFrame(data_before, columns=['FORMAT', 'SAMPLE'])
    result = trim_format(df_before)

    data_after = [["GT:AD:DP", "0:78:8"],["GT:DP", "90:0/9"], ["AD:GT:DP", "hi:hey:howdy"], ["DP", "89"], [".", "."], ["GT", "88"],  ["AD:DP","0/0:90" ] ]
    df_after = pd.DataFrame(data_after, columns=['FORMAT', 'SAMPLE'])

    assert result.equals(df_after)

def test_remove_no_alternates():
    data_before = {"ALT":["G", "C", "A", ".", ".", "C", "T", "A", "."], "TEST":[1,2,3,4,5,6,7,8,9]}
    df_before = pd.DataFrame(data_before)
    result = remove_no_alternate(df_before).reset_index(drop=True)
    

    data_after = {"ALT":["G", "C", "A", "C", "T", "A"], "TEST":[1,2,3,6,7,8]}
    df_after = pd.DataFrame(data_after).reset_index(drop=True)

    assert result.equals(df_after)

def test_filter_meta_info():
    mi_before = ['##fileformat=VCFv4.0\n##fileDate=20090805',
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FILTER=<ID=q10,Description="Quality below 10">',
        '##FILTER=<ID=s50,Description="Less than 50% of samples have data">',
        '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">']
    
    mi_after = ['##fileformat=VCFv4.0\n##fileDate=20090805',
        '##FILTER=<ID=q10,Description="Quality below 10">',
        '##FILTER=<ID=s50,Description="Less than 50% of samples have data">',
        '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="The allele frequency">']
    
    result = filter_meta_info(mi_before)

    assert result == mi_after

def test_normalize_vcf():
    normalize_vcf("data/source_test1.vcf", "data/result_test1.vcf")
    assert filecmp.cmp("data/out_test1.vcf", "data/result_test1.vcf")

test_remove_no_alternates()
