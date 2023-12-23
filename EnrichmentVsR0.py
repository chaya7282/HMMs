import pandas as pd
import numpy as np
import main
import matplotlib.pyplot as plt
import meme_orig
import logomaker as lm
from Bio import motifs as mtf
from numpy import log2
from sklearn.preprocessing import normalize
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

def plotPssmScores(df,motif_, enricment_col_name,seqs_experiment_name,pssm_experiment_name,fig_save_path):
    motif_scores = [motif_.pssm.calculate(row['kmer']) for  index, row in  df.iterrows()]
    scipy.pearsonr(list1, list2)
.
    plt.scatter( [np.log10(x) for x in df[enricment_col_name].tolist()],motif_scores)
    plt.ylabel('Pssm-score')
    plt.xlabel('Log10 (Enrichment)')

    plt.title("Pssm- "+ pssm_experiment_name + "Enrichment "+seqs_experiment_name )
    plt.savefig(fig_save_path)
    plt.tight_layout()
    plt.close()

def plotTwoPssmRatioScores(df,motif_seqs,motif_second, enricment_col_name,seqs_experiment_name,pssm_experiment_name,fig_save_path):
    motif_scores = [motif_seqs.pssm.calculate(row['kmer']) /motif_second.pssm.calculate(row['kmer']) for  index, row in  df.iterrows()]


    plt.scatter( [np.log10(x) for x in df[enricment_col_name].tolist()],motif_scores)
    earson_coef, p_value = stats.pearsonr(df_counts[name + "-enrichment"], df_counts[R0])
    plt.ylabel('Pssm-score')
    plt.xlabel('Log10 (Enrichment)')

    plt.title("Pssm- "+ pssm_experiment_name + "Enrichment "+seqs_experiment_name )
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
    count_matrix_path = input_path + "Counts for " + seqs_experimen_name + ".cvs"
    df_counts_=pd.read_csv(count_matrix_path, sep="\t")
    save_xls_path = input_path + "pwm for " + pssm_experiemt_name + ".xlsx.cvs"
    motif_counts= pd.read_csv( save_xls_path)
    pfm = {'A': motif_counts['A'],
           'C': motif_counts['C'],
           'G': motif_counts['G'],
           'T': motif_counts['T']}
    m = mtf.Motif(counts=pfm)
    fig_save_path= output_dir + "Scores for pssm of "+pssm_experiemt_name +"for enrichment of " +seqs_experimen_name
    plotPssmScores(df_counts_[::100], m,'R1vsR0', seqs_experimen_name,pssm_experiemt_name, fig_save_path)

if __name__ == '__main__':
    analyzePssm( 'R1-6nM','R1-60nM')

    r0_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R0.Thermo.table'
    r1_apth='C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.60nM.table'
    plotEnrichmentVsR0(r0_path ,r1_apth,'R1-60nM','R0-thermo')
    r1_apth='C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.6nM.table'
    plotEnrichmentVsR0(r0_path ,r1_apth,'R1-6nM','R0-thermo')
    r1_apth='C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.20nM.table'
    plotEnrichmentVsR0(r0_path ,r1_apth,'R1-20nM','R0-thermo')
    r0_path = 'C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R0.kinetics.table'
    r1_apth='C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.mix5.table'
    plotEnrichmentVsR0(r0_path ,r1_apth,'R1-mix5','R0-kinetics')
    r1_apth='C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.mix10.table'
    plotEnrichmentVsR0(r0_path ,r1_apth,'R1-mix10','R0-kinetics')
    r1_apth='C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.mix25.table'
    plotEnrichmentVsR0(r0_path ,r1_apth,'R1-mix25','R0-kinetics')
    r1_apth='C:/Users/User/PycharmProjects/extract_afinities/20mer_data_threshold 0/R1.mix50.table'
    plotEnrichmentVsR0(r0_path ,r1_apth,'R1-mix50','R0-kinetics')