import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import random




if __name__ == '__main__':
    mm_order = 4
    input_dir='C:/Users/User/PycharmProjects/HMMs/R0biasCorrection'
    output_dir = "R0biasCorrection/Thermo-bias-correction/"
    counts_file = 'C:/Users/User/PycharmProjects/HMMs/R0biasCorrection/Thermo-bias-correction/R26nM-Counts.cvs'
    df_counts_ = pd.read_csv(counts_file, sep="\t")
    df_counts=  df_counts_[df_counts_['ObservedCount-R26nM'] > 0].sort_values(by=['ObservedCount-R26nM'],ascending=True)
    experimen_name= 'R26nM'
    trail_name=['observed','non_homogenous','homogenous']
    description = ['R2:'+experimen_name+' R0-Thermo observed', 'R2:'+experimen_name+' R0-Thermo expected-non-homogenous mm='+str(mm_order-1),'R2:'+experimen_name+' R0--expected-homogenous mm='+str(11)]
    R2_=['ObservedCount-R26nM','ObservedCount-R26nM','ObservedCount-R26nM']
    R0_=['observed R0','Expected-non-homogenous-R0 mm=3','Expected-homogenous-R0 mm=11']
    for name,R2,R0, trail in zip(description,R2_,R0_,trail_name):
        sum_r2=df_counts[R2].sum()
        sum_R0=df_counts[R0].sum()
        print(R0)
        print(sum_r2)
        print(sum_R0)
        df_counts[trail+"-enrichment"]=[((row[R2]/sum_r2)/(row[R0]/sum_R0))for index, row in df_counts.iterrows() ]


        plt.scatter( df_counts[R0],df_counts[trail+"-enrichment"])
        plt.ylabel('Enrichment (f(R2)/f(R0))')
        plt.xlabel('R0')
        #pearson_coef, p_value = stats.pearsonr( df_counts[name+"-enrichment"], df_counts[R0])

        plt.title(name)
        plt.savefig(output_dir+" "+ trail +" "+experimen_name+'.png')
        plt.tight_layout()
        plt.close()
    df_counts.to_csv(counts_file, sep="\t")
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
if __name__ == '__main__':
