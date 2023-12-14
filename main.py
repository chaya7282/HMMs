# This is a sample Python script.
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import random

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
def readMultipleKmersTables(inputPath,excell_file):
    KmersTable={}
    with pd.ExcelWriter(excell_file) as writer:
        for file in os.listdir(inputPath):
            index=os.path.splitext(os.path.basename(file))[0]
            print(index)
            KmersTable[str(index)]= pd.read_csv(os.path.join(inputPath, file), sep='\t', index_col=False)

            KmersTable[str(index)]['frequency']=KmersTable[index]['ObservedCount']/sum(KmersTable[index]['ObservedCount'])
            KmersTable[str(index)].to_excel(writer, sheet_name="ofsset="+str(index))

    return KmersTable

def compute_mmExtended(seq,win_size,kmers_tables,info_tables,columns_tables):

    kmers=[seq[i:i+win_size] for i in range(len(seq) - win_size + 1)]
    calc_dict = []

    info_array=[]
    column_array=[]
    info_array.append(seq)
    column_array.append("kmer")
    mm_score=1

    if win_size==1:
        df = kmers_tables[str(win_size)][str(0)]
        row = df[df['Kmer'] == seq]
        return row['frequency'].item()

    if 1< win_size:
        p_shorter=compute_mm( kmers[0][:-1],win_size-1,kmers_tables,info_tables,columns_tables)
        mm_score = p_shorter*mm_score
        info_array.append(kmers[0][:-1])
        column_array.append("short mm kmer")
        info_array.append(p_shorter)
        column_array.append('p_value mm='+ str(win_size-1))

    for index, kmer in enumerate(kmers):

        df_long = kmers_tables[str(win_size)][str(index)]
        df_short=kmers_tables[str(win_size-1)][str(index)]

        row_long = df_long[df_long['Kmer'] == kmer]
        row_short = df_short[df_short['Kmer'] == kmer[:-1]]

        curr_p=row_long['frequency'].item()/row_short['frequency'].item()
        mm_score=mm_score*curr_p

        info_array.append(kmer)
        info_array.append(kmer[:-1])
        info_array.append(row_long['frequency'].item())
        info_array.append(row_short['frequency'].item())
        info_array.append(curr_p)
        column_array.append("long-kmer pos="+str(index))
        column_array.append("short-kemr pos=" + str(index))
        column_array.append("long-f " + str(index))
        column_array.append("short-f " + str(index))
        column_array.append("ratio pos-" + str(index))

    info_array.append( mm_score)
    column_array.append('mm_score')
    info_tables[str(win_size)].append(info_array)
    columns_tables[str(win_size)]=column_array
    return mm_score


# Press the green button in the gutter to run the script.
def compute_mm(seq,win_size,kmers_tables,info_tables,columns_tables):
    mm_order=win_size-1
    kmers=[seq[i:i+win_size] for i in range(len(seq) - win_size + 1)]
    calc_dict = []

    info_array=[]
    column_array=[]
    info_array.append(seq)
    column_array.append("kmer")
    mm_score =0
    for index, kmer in enumerate(kmers):

        df_long = kmers_tables[str(win_size)][str(index)]
        df_short=kmers_tables[str(win_size-1)][str(index)]

        row_long = df_long[df_long['Kmer'] == kmer]
        row_short = df_short[df_short['Kmer'] == kmer[:-1]]
        if index==0:
            mm_score = row_short['frequency'].item()
            column_array.append("first kmer")
            info_array.append(kmer[:-1])
            column_array.append("f(first)")
            info_array.append(np.round(row_short['frequency'].item(),3))
        curr_p=row_long['frequency'].item()/row_short['frequency'].item()
        mm_score=mm_score*curr_p

        info_array.append(kmer)
        info_array.append(kmer[:-1])
        info_array.append(row_long['frequency'].item())
        info_array.append(row_short['frequency'].item())
        info_array.append(np.round(curr_p,3))
        column_array.append("seq" + str(index) + "(" + str(mm_order+1) + ")")
        column_array.append("seq" + str(index) + "(" + str(mm_order) + ")")
        column_array.append("f"+str(index)+"("+str(mm_order+1)+")")
        column_array.append("f"+str(index)+"("+str(mm_order)+")")
        column_array.append("ratio-" + str(index))

    info_array.append( mm_score)
    column_array.append('probabaility')
    info_tables[str(win_size)].append(info_array)
    columns_tables[str(win_size)]=column_array
    return mm_score




def computeExpectedR0(mm_order,experiment_desc,train_path,test_path,output_path):
    train="Thermo"
    test="Thermo"

    kmers_tables={}
    info_tables={}
    columns_tables={}
    for i in range(1,mm_order+1):
        kmers_tables[str(i)]= readMultipleKmersTables(train_path +str(i),'C:/SELEX_workspace/Thermo_kmer_counts-offset/all-'+str(i)+'-kmers-offsets.xlsx')
        info_tables[str(i)]=[]
        columns_tables[str(i)] =[]
    Test_kmers_ = readMultipleKmersTables(test_path,'C:/SELEX_workspace/Thermo_kmer_counts-offset/0-'+str(test)+'.xlsx')

    probs = []
    Test_kmers = Test_kmers_[str(0)].copy()

    array_column=[]
    #pos_suffle(Test_kmers['Kmer'])

    for index, row in Test_kmers.iterrows():
        compute_mm(row['Kmer'], mm_order, kmers_tables,info_tables,columns_tables)

    for i in range(1,mm_order+1):
        if info_tables[str(i)]:
            result=pd.DataFrame(info_tables[str(i)],columns=columns_tables[str(i)])
            if i== mm_order:
                result['ObservedCount']=Test_kmers['ObservedCount'].tolist()
                result['Expected'] = result['probabaility'] * sum(result['ObservedCount'])
                result=result.drop_duplicates(subset=["kmer"])

                result.to_csv(output_path+experiment_desc+'.cvs',sep="\t")

def check_kmer_count(kmer):
    df = pd.read_csv('R0.thermo.table',sep="\t")
    count=0
    for index, row in df.iterrows():
        if row['Kmer'][1:].startswith(kmer):
            count=count+row['Kmer_count']
    print(kmer,count)


if __name__ == '__main__':
    #check_kmer_count('AACA')
    mm_order=4
    test='kinetics'
    for mm_order in np.arange(4,7,1):
        print(mm_order)
        output_path="train-thermo-suffle- test-kinetics/"
        experiment_desc = "train-thermo--test-kinetics"+" mm-order = "+str(mm_order)
        train_path='C:/SELEX_workspace/Thermo_kmer_counts-offset/'
        test_path = 'C:/SELEX_workspace/Thermo_kmer_counts-offset/18kmers-kinetics/'
        computeExpectedR0(mm_order,experiment_desc,train_path,test_path,output_path)
        df=pd.read_csv(output_path+experiment_desc+'.cvs', sep="\t")


        plt.scatter(df['ObservedCount'],df['Expected'])
        plt.xlabel('ExpectedCount')
        plt.xlabel('ObservedCount')
        pearson_coef, p_value = stats.pearsonr(df['ObservedCount'],df['Expected'])
        plt.title(experiment_desc+' pearson='+str(np.round(pearson_coef,3)))
        plt.savefig(output_path+experiment_desc+str(mm_order)+'.png')


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
