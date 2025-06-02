import sys
from Bio import SeqIO

def count_stretch(nt: str, seq: str, repeat: int):
    c=0
    counter=0
    stretch=nt*repeat
    while (c <= len(seq)):
        sub_seq=seq[c:c+repeat]
        if sub_seq != stretch:
            c+=1
        else:
            start_point=c+repeat
            counter+=1
            while (start_point < (len(seq))):
                if seq[start_point] != nt:
                    #print (seq[c:start_point])
                    c=start_point+1
                    break
                else:
                    if start_point != (len(seq)-1):
                        start_point+=1
                    else:
                        #print ((seq[c:start_point+1]))
                        c=start_point+1
                        break
            c=start_point+1
    return counter

# Check for correct number of arguments
if len(sys.argv) != 4:
    print("Usage: python script.py <input_file> <output_file> <Repeat_number>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
repeat_number = int(sys.argv[3])


o1=open(output_file,'w')
line0=['ID','Input-sequence-length','A','T','G','C']
o1.write("\t".join(line0)+"\n")

c=0
for data in SeqIO.parse(input_file,'fasta'):
    c+=1
    header=str(data.description)
    enst_id=(header.strip().split()[0]).strip().split('_')[-1]
    #print (enst_id)
    s=str(data.seq)
    nts=['A','T','G','C']
    nts_count=[]
    for char in nts:
        repeat_count=count_stretch(char,s,repeat_number)
        nts_count.append(repeat_count)
    o1.write(enst_id+"\t"+str(len(s))+"\t"+"\t".join(map(str,nts_count))+"\n")
    if c% 1000 == 0:
        print ('Checkpoint:' +str(c))
o1.close()