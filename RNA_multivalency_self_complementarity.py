import os
import sys
from Bio import SeqIO
import math
from scipy.stats import pearsonr
import datetime

def find_complementary(subseq:str):
    complentary_dict={'A':'U','G':'C','C':'G','U':'A'}
    comp_seq=""
    for char in subseq:
        for key, val in complentary_dict.items():
            if char == key:
                comp_seq+=val
    return comp_seq


'''
Function to get the complementary sequences if the other strand has a T/U
Basically U can base pair with both A and G. So, if a complementary sequence has A, 
we need all combinations of by replacing A with G in the complementary sequence.

PMID : 11256617
'''
def get_combinations(original_str, replacements_dict, results_list, current_combination):
    if not original_str:
        results_list.append(current_combination)
    else:
        first_char = original_str[0]
        remaining_str = original_str[1:]
        if first_char in replacements_dict:
            for replacement in replacements_dict[first_char]:
                get_combinations(remaining_str, replacements_dict,
                                 results_list, current_combination + replacement)
        else:
            get_combinations(remaining_str, replacements_dict,
                             results_list, current_combination + first_char)

def scan_complement(subseq:str,complement_seq:str):
    c=0
    for i in range (len(subseq)):
        char_search = subseq[i:i+(len(complement_seq))]
        if char_search ==  complement_seq:
            c+=1
    return c

# Check for correct number of arguments
if len(sys.argv) != 3:
    print("Usage: python script.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

o1=open(output_file,'w')
line0=['Sequence-id','Length_sequence','Count','Stickiness_by_length','Stickiness_per_1000nts']
o1.write("\t".join(line0)+"\n")

now = datetime.datetime.now()
print (now)
stickiness_score_list=[]
sticky_count_list=[]
input_length=[]
icount=0
for data in SeqIO.parse(input_file,"fasta"):
    header=str(data.id)
    #print (header)
    input_seq_dna=str(data.seq)
    input_seq_rna=input_seq_dna.replace('T','U')
    #print (len(utr3_seq_rna))
    c=0
    window_5mer=0
    sticky_count=0
    while (c < len(input_seq_rna)-5):
        subseq=input_seq_rna[c:c+5]
        window_5mer+=1
        complement_seq_rev=find_complementary(subseq)[::-1]
        #print (subseq+"----"+complement_seq_rev+"\n")
        #upstream_start=c-3
        downstream_start=(c+5)+2
        #count_upstream=0
        count_downstream=0

        '''Only calculate complementary regions using watson-crick pairs'''
        '''if downstream_start <= (len(utr3_seq_rna)-5):
            downstream_seq_region = utr3_seq_rna[downstream_start:]
            count_downstream=scan_complement(downstream_seq_region,complement_seq_rev)'''
        
        '''Below if-else block is calculate complementary regions considering wobble bases
        U can pair with G along with normal watson-crick pairs'''
        
        if 'A' not in complement_seq_rev:
            if downstream_start <= (len(input_seq_rna)-5):
                downstream_seq_region=input_seq_rna[downstream_start:]
                count_downstream=scan_complement(downstream_seq_region,complement_seq_rev)

            sticky_count+=count_downstream            
        else:
            wobble_dict={'A':['A','G']}
            res=[]
            get_combinations(complement_seq_rev, wobble_dict, res, "")
            if downstream_start <= (len(input_seq_rna)-5):
                count_downstream=0
                for i1 in range (len(res)):
                    downstream_seq_region=input_seq_rna[downstream_start:]
                    count_downstream+=scan_complement(downstream_seq_region,res[i1])
            
        sticky_count+=count_downstream    
        c+=5  #Shift by 1 nucleotide,overlapping window ; NB - also run a non-overlapping ie c+=5

    try:
        stickiness_score_by_window_freq=sticky_count/window_5mer
        stickiness_score_by_utr_length=sticky_count/len(input_seq_rna)
        stickiness_score_per_1000_nt=(sticky_count/len(input_seq_rna))*1000
        #print (stickiness_score)
        #print (window_5mer)
        stickiness_score_list.append(stickiness_score_by_window_freq)
        sticky_count_list.append(sticky_count)
        input_length.append(len(input_seq_rna))
        o1.write(header+"\t"+str(len(input_seq_rna))+"\t"+str(sticky_count)+
                 "\t"+str(stickiness_score_by_utr_length)+"\t"+str(stickiness_score_per_1000_nt)+"\n")

    except ZeroDivisionError:
        o1.write(header+"\t"+str(len(input_seq_rna))+"\t"+"input length less than 5"+"\n")
    
    icount+=1
    if icount % 1000 == 0:
        now = datetime.datetime.now()
        print (str(now)+"\t"+str(icount))

o1.close()
pcorr=pearsonr(input_length,stickiness_score_list)
print ("Pearsonr : "+str(pcorr.statistic)+"with p-value of : "+str(pcorr.pvalue))

pcorr1=pearsonr(input_length,sticky_count_list)
print ("Pearsonr : "+str(pcorr1.statistic)+"with p-value of : "+str(pcorr1.pvalue))