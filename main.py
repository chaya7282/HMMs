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

                result.to_csv(output_path,sep="\t")

def check_kmer_count(kmer):
    df = pd.read_csv('R0.thermo.table',sep="\t")
    count=0
    for index, row in df.iterrows():
        if row['Kmer'][1:].startswith(kmer):
            count=count+row['Kmer_count']
    print(kmer,count)


def mergeTwoAfinnityRounds(file_1, file_2, file_1_indicator, file_2_indicator, output_file):
    df1 = pd.read_csv(file_1, sep="\t", index_col=None)
    df1.columns = [str(col) + file_1_indicator if str(col) != 'kmer' else str(col) for col in df1.columns]
    df2 = pd.read_csv(file_2, sep="\t", index_col=None)
    df2.columns = [str(col) + file_2_indicator if str(col) != 'kmer' else str(col) for col in df2.columns]
    df4 = pd.merge(df1, df2, on='kmer', indicator=True, how='outer')
    df4.replace('', np.nan, inplace=True)

   # df4 = df4.dropna()
    df4.to_csv(output_file,sep="\t")

def mergeR1WithR2(ile_1, file_2, file_1_indicator, file_2_indicator, output_file):

    mergeTwoAfinnityRounds(os.path.join(input_dir, 'affinities-R6nM_R1.csv'), os.path.join(input_dir, 'affinities-R6nM_R2.csv'), "R1","R2", os.path.join(output_dir, "6nMR1-R2.csv"))
    mergeTwoAfinnityRounds(os.path.join(input_dir, 'affinities-R20nM_R1.csv'), os.path.join(input_dir, 'affinities-20nM_R2.csv'),"R1", "R2", os.path.join(output_dir, "20nMR1-R2.csv"))
    mergeTwoAfinnityRounds(os.path.join(input_dir, 'affinities-60nM_R1.csv'),os.path.join(input_dir, 'affinities-60nM_R2.csv'), "R1", "R2",os.path.join(output_dir, "60nMR1-R2.csv"))


def corrct_bias():
    #check_kmer_count('AACA')
    mm_order=4
    test='kinetics'
    input_dir= 'raw_files/Thermo-tmp'
    output_dir = "R0biasCorrection/Thermo-bias-correction/"
    experiment_desc='mm based on Thermo'+" mm order="+str(mm_order-1)
    train_path = 'C:/SELEX_workspace/Thermo_kmer_counts-offset/'
    r0_path='C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R0.thermo.table'
    for experiment in  [ os.path.join(input_dir, file) for file in os.listdir(input_dir)]:

        print(mm_order)
        experimen_name=os.path.splitext(os.path.basename(experiment))[0]
        test_path = input_dir+"/"+ experimen_name+"/"
        output_path = output_dir + experiment_desc
        computeExpectedR0(mm_order,experiment_desc+" "+experimen_name,train_path,test_path,output_path+".cvs")

        df=pd.read_csv(output_path+".cvs", sep="\t")
        counts_file = output_dir + '' + experimen_name + "-Counts" + ".cvs"
        df[['kmer','ObservedCount','Expected']].to_csv(counts_file, sep="\t")

        mergeTwoAfinnityRounds(counts_file, r0_path,  "-" +experimen_name,"-" +'R0',counts_file )
        homogenous_mm_flie= 'raw_files/R26nMHomogenous/20.cvs'
        mergeTwoAfinnityRounds(counts_file, homogenous_mm_flie, "-non_homogenous", experimen_name+'-homogenous', counts_file)
        df_counts = pd.read_csv(counts_file, sep="\t")
        df_counts =df_counts[df_counts['_merge-non_homogenous']=='both']
        df_counts = df_counts[df_counts['_merge'] == 'both']
        trail_name=['observed','non_homogenous','homogenous']
        description = ['R2-div-R0-observed', 'R2-div-R0-expected-non-homogenous mm='+str(mm_order-1),'R2-div-R0-expected-homogenous mm='+str(11)]
        R2_=['ObservedCount-R26nM-non_homogenous','ObservedCount-R26nM-non_homogenous','ObservedCount-R26nM-non_homogenous']
        R0_=['ObservedCount-R0-non_homogenous','Expected-R26nM-non_homogenous','ExpectedCountR26nM-homogenous']
        for name,R2,R0 in zip(description,R2_,R0_):
            sum_r2=df_counts[R2].sum()
            sum_R0=df_counts[R0].sum()
            print(R0)
            print(sum_r2)
            print(sum_R0)
            df_counts[name+"-enrichment"]=[((row[R2])/(sum_r2))/((row[R0])/(sum_R0)) for index, row in df_counts.iterrows() ]


            plt.scatter( df_counts['ObservedCount-R0-non_homogenous'],df_counts[name+"-enrichment"])
            plt.ylabel('Enrichment')
            plt.xlabel('R0')
            pearson_coef, p_value = stats.pearsonr( df_counts[name+"-enrichment"], df_counts[R0])

            plt.title(experimen_name+" "+ name +"\n" +'pearson='+str(np.round(pearson_coef,3)))
            plt.savefig(output_dir+" "+ name +" "+experimen_name+'.png')
            plt.tight_layout()
            plt.close()
        df_counts.to_csv(counts_file, sep="\t")
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
if __name__ == '__main__':
    corrct_bias()