import pandas as pd
import numpy as np
import main
import matplotlib.pyplot as plt
import meme_orig
import logomaker as lm
from Bio import motifs as mtf
from numpy import log2
from scipy import stats
from sklearn.preprocessing import normalize
import functools
def get_logo_from_table_counts(df_table_count_,sort_column,experiment_column,row_index_thresh,ascending):
    df_table_count=df_table_count_.sort_values(by=[sort_column],ascending=ascending)
    count_mat=np.zeros((len(df_table_count['kmer'][0]),4), dtype='float')
    index=0
    sum_count=0
    for index2, row in df_table_count.iterrows():
        count_mat += meme_orig.makeCountMatrix(row['kmer']).T*row[sort_column]

        index +=1
        if  row_index_thresh < index:
            break
    pfm = {'A': count_mat[:,0] ,
           'C': count_mat[:,1],
           'G': count_mat[:,2],
           'T': count_mat[:,3]}

    m = mtf.Motif(counts=pfm)
    return m

def IC_matrix(m_bio):

    arr= np.array([m_bio['A'],m_bio['C'],m_bio['G'],m_bio['T']])

    logo_mat =arr
    arr = arr
    U=  arr * log2(arr)
    U= np.nan_to_num(U)
    IC= np.sum(U, axis=0)
    IC+=2


    for row in range(arr.shape[0]):
        for col in range(arr.shape[1]):
            logo_mat[row][col]=arr[row][col]*IC[col]


    return logo_mat,  sum(IC)

def plotPssmScores(df,ro_nmae,experiment_name,fig_save_path):
    input_path = 'C:/Users/User/PycharmProjects/HMMs/R0biasCorrection/Thermo-bias-correction/data/'


    motif_ = read_motif(input_path, experiment_name)
    enrichments = np.array(df['ObservedCount_' + experiment_name]) / np.array(df['ObservedCount_' + ro_nmae])

    motif_scores = [motif_.pssm.calculate(row['kmer']) for  index, row in  df .iterrows()]
    df_all['pssm ' + experiment_name] =motif_scores
    df_all['enrichment ' + experiment_name] = enrichments
    pearson_coef, p_value = stats.pearsonr( motif_scores,np.log10(enrichments))
    print(pearson_coef)
    plt.scatter(motif_scores,np.log10(enrichments))
    plt.ylabel('Log10 (Enrichment)')
    plt.xlabel('pssm-score')

    plt.title("Experiment= " +experiment_name+ " pearson="+str(np.round(pearson_coef,3)) )
    plt.savefig(fig_save_path)
    plt.tight_layout()
    plt.close()
    return df_all
def read_motif (input_path,seqs_experimen_name):
    save_xls_path = input_path + "pwm for " + seqs_experimen_name + ".xlsx.cvs"
    motif_counts= pd.read_csv( save_xls_path)
    pfm = {'A': motif_counts['A'],
           'C': motif_counts['C'],
           'G': motif_counts['G'],
           'T': motif_counts['T']}
    m1 = mtf.Motif(counts=pfm)
    return m1
def plotTwoPssmRatioScores(df_all,motif_first,motif_second,ro_nmae, first_exp_name,second_exp_name,fig_save_path):
    input_path='C:/Users/User/PycharmProjects/HMMs/R0biasCorrection/Thermo-bias-correction/data/'

    motif_first_=read_motif(input_path,motif_first)
    motif_second_= read_motif(input_path, motif_second)
    motif_scores = [motif_first_.pssm.calculate(row['kmer']) /motif_second_.pssm.calculate(row['kmer']) for  index, row in  df_filtered.iterrows()]
    enrichments= np.array(df_all[first_exp_name])/np.array(df_all[second_exp_name])

    plt.scatter(enrichments,motif_scores)

    plt.ylabel('Pssm-score')
    plt.xlabel('Log10 (Enrichment)')

    plt.title("Pssm-Score "+first_exp_name +"/Pssm-Score "+ second_exp_name  + "Enrichment "+ first_exp_name+"/"+second_exp_name  )
    plt.savefig(fig_save_path)
    plt.tight_layout()
    plt.close()
def create_save_logo(mot,save_fig_path,save_xls_path,title):
    logo_mat, IC = IC_matrix(mot.pwm)
    df = pd.DataFrame(mot.counts, columns=['A', 'C', 'G', 'T'])
    df.to_csv(save_xls_path+".cvs")
    df = pd.DataFrame(logo_mat.T, columns=['A', 'C', 'G', 'T'])
    fig, ax = plt.subplots()
    crp_logo = lm.Logo(df, vpad=0)
    crp_logo.ax.set_xticks(range(20))
    crp_logo.ax.set_xticklabels('%d' % x for x in range(20))
    plt.title(title)
    ax.set_xticks(range(20))
    ax.set_xticklabels(range(20))
    plt.savefig(save_fig_path+".png")
    plt.close()
def plotEnrichmentVsR0(r0_path ,r1_apth,experimen_name,R0_name):

    threshold=100
    output_dir = "R0biasCorrection/Thermo-bias-correction/Logos/"

    counts_file="R0biasCorrection/Thermo-bias-correction/"+experimen_name+'.cvs'
    main.mergeTwoAfinnityRounds(r1_apth, r0_path, "-" + experimen_name, "-" + 'R0', counts_file)
    df_counts_=pd.read_csv(counts_file, sep="\t")
    df_counts_ = pd.read_csv(counts_file, sep="\t")
    df_counts_=    df_counts_[df_counts_['_merge'] =='both']
    sum_r1=df_counts_['ObservedCount-' + experimen_name].sum()
    sum_r0 =df_counts_['ObservedCount-R0'].sum()
    df_counts_['R1vsR0']= [(row['ObservedCount-' + experimen_name]/sum_r1)/(row['ObservedCount-R0']/sum_r0) for index, row in df_counts_.iterrows() ]
    fig_save_name= "Sort by enrichment-"+experimen_name +" 1000-Top"
    save_fig_path = output_dir + fig_save_name + ".png"
    save_xls_path=output_dir + "pwm for " +experimen_name+ ".xlsx"
    m_biopythonR1 = get_logo_from_table_counts(df_counts_, 'R1vsR0','ObservedCount-' + experimen_name, 1000, False)
    create_save_logo( m_biopythonR1, save_fig_path, save_xls_path, fig_save_name)

    count_matrix_path=output_dir + "Counts for " +experimen_name+ ".cvs"
    df_counts_.to_csv(count_matrix_path,sep="\t", )
    fig_save_path= output_dir + "pssm vs Enrichment for "+experimen_name+ ".png"
    plotPssmScores(df_counts_[::100], m_biopythonR1, 'R1vsR0', experimen_name, fig_save_path)

    t=7

def analyzePssm(seqs_experimen_name,pssm_experiemt_name):
    output_dir = "R0biasCorrection/Thermo-bias-correction/Logos/"
    input_path="R0biasCorrection/Thermo-bias-correction/data/"

    df_counts_1['ObservedCount_'+experiment_name]=[row['R1vsR0-6nM']/row['R1vsR0-60nM'] for index, row in df_counts_1.iterrows()]




    save_xls_path = input_path + "pwm for " +pssm_experiemt_name + ".xlsx.cvs"
    motif_counts = pd.read_csv(save_xls_path)
    pfm = {'A': motif_counts['A'],
           'C': motif_counts['C'],
           'G': motif_counts['G'],
           'T': motif_counts['T']}

    m2 = mtf.Motif(counts=pfm)


    fig_save_path= output_dir + "score vs enrich divisions" +seqs_experimen_name+"vs" +pssm_experiemt_name

    plotTwoPssmRatioScores(df_counts_1[::100], m1, m2, 'R1vsR0', seqs_experimen_name,pssm_experiemt_name, fig_save_path)
def merge_all_files():
    r0_thermo_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R0.Thermo.table'
    r0_kinetics_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R0.kinetics.table'
    r1_60nM_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.60nM.table'
    r1_6nM_path= 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.6nM.table'
    r1_20nM_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.20nM.table'

    r1_mix5_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.mix5.table'
    r1_mix10_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.mix10.table'
    r1_mix25_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.mix25.table'
    r1_mix50_path= 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.mix50.table'
    file_list=[r0_thermo_path,r0_kinetics_path,r1_60nM_path,r1_6nM_path,r1_20nM_path,r1_mix5_path,r1_mix10_path,r1_mix25_path,r1_mix50_path]
    experiment_nemse= ['R0_thermo','R0_kinetics', 'R1_6nM','R1_60nM','R1_20nM','R1_mix5','R1_mix10','R1_mix25','R1_mix50']
    data_frames =[]
    for exp_name,file in zip(experiment_nemse,file_list):
        df=pd.read_csv(file, sep="\t")[['kmer',	'ObservedCount']]
        df.rename(columns={'ObservedCount':'ObservedCount_'+exp_name}, inplace=True)
        data_frames.append(df)
    df_merged = functools.reduce(lambda left, right: pd.merge(left, right, on='kmer', how='outer'), data_frames)
    df_merged.to_csv('Allconditions.cvs', sep="\t")
    i=8
if __name__ == '__main__':
    merge_all_files()
    df_all=pd.read_csv('Allconditions.cvs', sep="\t")
    df_all = df_all.replace( '',np.nan, regex=True)
    df_all=df_all.dropna()
    fig_save_path='tmp_fig.png'
    fig_save_path='pssmVsEnrichmentGraphs/R1.6nm.png'
    df_all=plotPssmScores(df_all, 'R0_thermo', 'R1_6nM', fig_save_path)
    fig_save_path='pssmVsEnrichmentGraphs/R1.20nm.png'
    df_all=plotPssmScores(df_all, 'R0_thermo', 'R1_20nM', fig_save_path)
    fig_save_path = 'pssmVsEnrichmentGraphs/R1.60nm.png'
    df_all=plotPssmScores(df_all, 'R0_thermo', 'R1_60nM', fig_save_path)
    df_all.to_csv('AllDatawithpssm.cvs', sep="\t")


