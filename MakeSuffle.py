# This is a sample Python script.
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import random
from Bio.Seq import Seq
from numpy.random import choice
from Bio import motifs

def suff_between_strs(prefices):
    df_prefix = pd.DataFrame(prefices, columns=[str(x) for x in np.arange(len(prefices[0]))])
    for col in df_prefix.columns:
        x = df_prefix[col].tolist()
        np.random.shuffle(x)
        df_prefix[col] = x
        print(col)
    prefices_shuff = [''.join(map(str, x)) for x in df_prefix.values.tolist()]
    return prefices_shuff
def seqs_suffle_speed(seqs: object) -> object:
    prefices= [seq[:len(seq)//2].split("")  for seq in seqs]
    suffices = [seq[len(seq) // 2:].split("")for seq in seqs]

    prefices_shuff=suff_between_strs(prefices)
    suffix_shuff = suff_between_strs(suffices)
    res = [i + j for i, j in zip(prefices_shuff, suffix_shuff)]
    return res
    i=5

def seqs_suffle(seqs: object) -> object:
    seqs_cpy = seqs.copy()
    seq_res=[]
    string_len=len(seqs[0])
    seqs_chars = [Seq(seq) for seq in seqs]
    m_1 = motifs.create(seqs_chars)
 #   for index in np.arange(string_len):
    for index in np.arange(string_len):
        indices =np.arange(len(seqs))
        random.shuffle(indices )
        print(index)

        for i, string in enumerate(seqs):
            if i ==6:
                print(seqs_cpy[i])
            x =  [t for t in seqs[i]]

            x[index]= seqs_cpy[indices[i]][index]
            seqs[i]="".join(x)
            if i ==0:
                print(seqs[i])

    #df = pd.DataFrame(seqs_chars, columns=[str(index) for index in np.arange(len(seqs_chars[0]))])
    #df_shuffled = df.copy()
    seqs_chars = [Seq(seq) for seq in seqs]
    m = motifs.create(seqs_chars)
    return seqs
    i = 8
if __name__ == '__main__':
    input_file='C:/Users/User/PycharmProjects/extract_afinities/RFormat/R0.thermo.txt'

    df=pd.read_csv(input_file, index_col=False)
    df.reset_index(drop=True, inplace=False)
    df_shuff=df.copy()
    x=df_shuff[df_shuff['+'] != '+'].index.values
    seqs= df_shuff[df_shuff['+'] != '+']['+'].tolist()

    seq_shuff=seqs_suffle_speed(seqs)
   # df_shuff.loc[x, '+'] = pd.Series([seq_shuff])
    for i, index in enumerate(x):
        df_shuff['+'][index]=seq_shuff[i]

    df_shuff['+'].to_csv('C:/SELEX_workspace/Suffle_seqs/thermoshafle-1.fastq2',index=False)
    df['+'].to_csv('C:/SELEX_workspace/Suffle_seqs/thermoOrig-1.fastq2', index=False)
    seq_shuff = seqs_suffle(seqs)
    for i, index in enumerate(x):
        df_shuff['+'][index] = seq_shuff[i]

    df_shuff['+'].to_csv('C:/SELEX_workspace/Suffle_seqs/thermoshafle-2.fastq2',index=False)
    i=5