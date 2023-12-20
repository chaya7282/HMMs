import pandas as pd
import numpy as np
import main
import matplotlib.pyplot as plt
import meme_orig
from Bio import motifs as mtf
def get_logo_from_table_counts(df_table_count_,measure_column,row_index_thresh,ascending):
    df_table_count=df_table_count_.sort_values(by=[measure_column],ascending=ascending)
    count_mat=np.zeros((len(df_table_count['kmer'][0]),4), dtype='float')
    index=0
    for index2, row in df_table_count.iterrows():
        count_mat += meme_orig.makeCountMatrix(row['kmer']).T
        index +=1
        if  row_index_thresh < index:
            break
    pfm = {'A': count_mat[:,0] ,
           'C': count_mat[:,1],
           'G': count_mat[:,2],
           'T': count_mat[:,3]}

    m = mtf.Motif(counts=pfm)
    return m

def plotEnrichmentVsR0(r0_path ,r1_apth,experimen_name,R0_name):

    threshold=100
    output_dir = "R0biasCorrection/Thermo-bias-correction/1InmillionLog/"

    counts_file="R0biasCorrection/Thermo-bias-correction/"+experimen_name+'.cvs'
    main.mergeTwoAfinnityRounds(r1_apth, r0_path, "-" + experimen_name, "-" + 'R0', counts_file)
    df_counts_=pd.read_csv(counts_file, sep="\t")
    df_counts_ = pd.read_csv(counts_file, sep="\t")
    df_counts_=    df_counts_[df_counts_['_merge'] =='both']

    threshold = df_counts_['ObservedCount-' + experimen_name].sum() / 500000
    sum_r1 = df_counts_['ObservedCount-' + experimen_name].sum()
    sum_R0 = df_counts_['ObservedCount-R0'].sum()
    df_counts =  df_counts_[(df_counts_['ObservedCount-' + experimen_name]) > threshold].sort_values(by=['ObservedCount-'+experimen_name],ascending=True)


    title_graph='R1:'+experimen_name+', R0:'+R0_name +'mean enrichment'

    df_counts[experimen_name + "-enrichment"] = [((row['ObservedCount-'+experimen_name])/(row['ObservedCount-R0'])) for index, row in df_counts.iterrows()]
    max_r1 = df_counts[experimen_name + "-enrichment"].max()
    #df_counts_bins=df_counts.copy()

    #df_counts_bins['bins'] = pd.cut(df_counts_bins['ObservedCount-R0'], bins=100)
    #df_counts_bins = df_counts_bins.groupby(['bins'])
    #bins_of_grous = [x.right for x in list(df_counts_bins.groups.keys())]
    plt.scatter(df_counts['ObservedCount-R0'], df_counts[experimen_name + "-enrichment"], label="count-thrshold="+str(np.round(threshold,2)))
    plt.title(title_graph, fontsize = 14)
    plt.ylim([0,30])
    plt.ylabel('#std enrichment ', fontsize = 12)
    plt.xlabel('R0-observed', fontsize = 12)
    plt.legend()
    # pearson_coef, p_value = stats.pearsonr( df_counts[name+"-enrichment"], df_counts[R0])

    plt.savefig(output_dir + 'Observed' + " " + experimen_name + '.png')
    plt.tight_layout()
    plt.close()

def draw_logos()
    save_fig_path = input_folder + "logo-1000 bottom " + fig_save_name + ".png"
    create_save_logo(m_biopython, save_fig_path, experiment_desc + " 1000 bottom")
    plt.close()

if __name__ == '__main__':
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