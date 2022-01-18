from Bio import SeqIO
import pandas as pd
import numpy as np
import re
from Bio.Seq import Seq


def mut_start_stop_codon(guide):
    if guide=='empty':
        return 'empty'
    else:
        stop_codon=['taa','tag','tga']
        start_codon='atg'
        pos = {}
        for seq in stop_codon:
         pos[seq] = [m.start() for m in re.finditer(seq,guide)]

        element=list(guide)
        for seq in stop_codon:
            for index in pos[seq]:
                if (index % 3 == 0):
                    if seq!='taa':
                        element[index:index + 3] = ['t','g','g']
                    else:
                        element[index:index + 3] = ['t', 'a', 'c']

        stop_pos=guide.find('A')
        atg_pos=[m.start() for m in re.finditer(start_codon,guide)]
        for index in atg_pos:
            if (index % 3 == 0 and index>stop_pos):
                element[index:index + 3] = ['a','g','g']
        element=''.join([str(item) for item in element])
        return element

ms2_seq1='agacatgaggatcacccatgt'
ms2_seq2='aagggtggaggaacaccccaccct'
ms2_seq3='acagaagcaccatcagggcttctg'
ms2_seq4='gtgcgtggagcatcagcccacgca'
ms2_seq5='tcgacgcaggaccaccgcgtc'
ms2_seq6='agcgcagaggaacaccctgcg'
ms2_seq7='acgggtggaggatcaccccacccg'
ms2_seq8='tcgcgaagagcatcagccttcgcg'


def make_sensor(Genename,length):
    filename='Endo DNA Files/{}.gb'.format(Genename)
    record = SeqIO.read(filename,"gb")
    target=record.seq
    target=str(target)
    target=target.lower()
    cca= [m.start() for m in re.finditer('cca',target)]
    Middle=[]
    target=list(target)
    for index in cca:
        beg=int(index-((length-3)/2))
        end_index=beg+length
        #print(beg)
        if beg <1 or (beg+length)>len(target):
            res='empty'
        else:
            frac=target[beg:beg+length]
            frac=''.join([item for item in frac])
            rev=Seq(frac).reverse_complement()
            rev=list(str(rev))
            a_pos=int((length-1)/2)
            rev[a_pos]='A'
            rev = ''.join([item for item in rev])
            #Here spacer length between avidity region is 5bp
            #Avidity region is 27bp
            right1=target[end_index+5:beg+length+32]
            right2=target[end_index+37:end_index+64]
            left1 = target[beg-32 : beg-5]
            left2 = target[beg - 64 : beg - 37]
            right1= str(Seq(''.join(right1)).reverse_complement())
            right2 = str(Seq(''.join(right2)).reverse_complement())
            left1 = str(Seq(''.join(left1)).reverse_complement())
            left2 = str(Seq(''.join(left2)).reverse_complement())
            res=right2+ms2_seq1+right1+ms2_seq2+rev+ms2_seq3+left1+ms2_seq4+left2

        Middle.append(res)

    df=pd.DataFrame()
    df['Middle']=pd.Series(Middle)
    df_new=df.copy()
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            df_new.iloc[i,j]=mut_start_stop_codon(df.iloc[i,j])
    #df_new.to_csv('sensor.csv',index=None)

    name_list=[]
    for index in range(df_new.shape[0]):
        name = Genename + '_CCA_' + str(index + 1) + '_ms2_five_avidity_region'
        name_list.append(name)
    # a=df_new.values
    # print(a.shape)

    res=pd.DataFrame(name_list,columns=['Name'])
    res=pd.concat([res,df_new],axis=1)
    outputname="Endo DNA Files/{}_ms2_sensor_guide.csv".format(Genename)
    res.to_csv(outputname,index=None)
    return res

gene_list=['hsp70','IFNb']
for i in gene_list:
    res=make_sensor(i,51)