#!/usr/bin/env python
# Author: Ioannis Anastopoulos
# Last Modified: 02/07/2019

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import argparse

# global colors for all plots
color_dict = {"FUSION" :"#7300E7",
       "TRUNCATING" : "#FF8C00", 
       "AMP" : "red",
       "DEL" : "blue",
       "DEEPDEL" : "#0000FF",
       "SHALLDEL" : "light blue",
       "MISSENSE":"#698B69",
       "PROMSNV":"#8B008B",
       "UTR3SNV":"#1A1A1A",
       "UTR5SNV":"#11E6EE",
       "SPLICE":"#00CED1",
       "CDS":"#008000",
       "SV":"#7300E7",
       "ENH":"#7F000F",
	"MISSENSE_SET1":"#698B69",
	"MISSENSE_SET2":"#698B69",
	"MISSENSE_SET3":"#698B69",
	"MISSENSE_SET4":"#698B69"}




def memoSort(M, p=None, p_encoder=None):
    # M will be a matrix with samples x genes pandas df
    # for the purposes of this plot the bait should be excluded from the matrix
    # p is a list of perturbations
    # if p is none, it is assumed that the df contains zero and ones
    # if no ValueError will be raised
    
    def scoreRow(x):
        score=0
        for i in range(len(x)):
            if x[i]!=0:
                score += 2**(len(x)-i)
        return score
    
    if p and p_encoder:
        M = M.replace(p_encoder)
    else:
        M = M.replace(dict(zip(p, np.arange(1,len(p)+1))))

    col_sum = M.sum(axis=0) #getting sum of columns(genes)
    col_sum_sort = col_sum.sort_values(ascending=False)
    
    M=M.T.loc[col_sum_sort.index.values].T #genes are sorted here
    
    sample_scores = {M.index[i]: scoreRow(sample) for i, sample in enumerate(M.values)}
    sample_scores = pd.Series(sample_scores).sort_values(ascending=False)
    
    M=M.loc[sample_scores.index] #sorting by sample score
    return M, col_sum_sort
    

    
def oncoprint(classifier_score, output_matrix, out, dpi, cmap='tab20b'):
    
    #colors to be used for the different annotations
    #USE same colors for legend of alterations
    
    #make a color wheel instead of hardcoding these
    
    # get colors
    
    classifier_score = classifier_score.sort_values(by='x', ascending=False) #sorting classifier scores
    output_matrix=output_matrix.T
    output_matrix=output_matrix.loc[classifier_score.index] #sorting output matrix in the same order as the sorted classifier score file
    output_matrix = output_matrix.fillna(0)
    
    perturbations=set()
    for c in output_matrix.columns:

        u=(output_matrix[c].unique())
        for perturbation in u:
            if 0 != perturbation:
                perturbations.add(str(perturbation))
    perturbations = sorted(list(perturbations))
    perturbation_encoder = dict(zip(perturbations, np.arange(1,len(perturbations)+1)))
    
    output_matrix = output_matrix.replace(perturbation_encoder)
    
    bar_width=1
    total_panels=len(output_matrix.columns)+1
    scores=classifier_score.x.values.tolist()

    if total_panels >12:
        total_panels=12

    ########----------------------Canvas--------------------##########
    width=10
    height=total_panels
    plt.figure(figsize=(width,height))  #units are always inches
    ########----------------------Canvas--------------------##########


    ########----------------------Main panel dimensions--------------------##########
    panel_width=3/width 
    panel_height=0.7/height
    classifier_height=0.85
    classifier_panel=plt.axes([0.1,classifier_height,panel_width,panel_height]) #this adds axis labels

    ########----------------------Main panel dimensions--------------------##########


    ########----------------------classifier panel--------------------##########
    for i,e in enumerate(scores):
        left=i
        bottom=0
        height=e
        width=bar_width
        rectangle=mplpatches.Rectangle((left,bottom),width,height,
                                               facecolor='black',
                                                edgecolor='black',
                                                linewidth=0) 
        classifier_panel.add_patch(rectangle)

    classifier_panel.set_xlim(0,len(scores))
    classifier_panel.set_ylim(0,1)
    classifier_panel.tick_params(axis='both', which='both',bottom=False, labelbottom=False, left=True, labelleft=True, 
            right=False, labelright=False, top=False, labeltop=False)
    classifier_panel.set_ylabel('Classifier'+'\n'+'Score', fontsize=10)

    classifier_panel.plot([0,len(scores)], [0.5, 0.5], linewidth=2, color='red', zorder=3)

    ########----------------------classifier panel--------------------##########
    c=0
    legend = set() 
    last_panel=None
    
    #reversing the encoder dict so its {int:annotation}
    encode_pert = dict()
    for k,v in perturbation_encoder.items():
        encode_pert[int(v)]=k
        
    for i,attribute in enumerate(output_matrix.columns):
        if c!=11:
            c+=1
        else:
            break

        if last_panel is None:
            last_panel=classifier_height

        next_panel=last_panel-1.05*panel_height
        panel=plt.axes([0.1,next_panel,panel_width,panel_height], frameon=False) #this adds axis labels
        last_panel=next_panel
        h=1

        panel.set_xlim(0,len(scores))
        panel.set_ylim(0,h)
        panel.tick_params(axis='both', which='both',bottom=False, labelbottom=False, left=False, labelleft=False, 
        right=False, labelright=False, top=False, labeltop=False)

        if i==0:
            panel.text(-0.8, h/2,'Bait', rotation = 0, fontsize=15, ha='right', va='top' )
        else:
            panel.text(-0.8, h/2,'Catch{}'.format(i), rotation = 0, fontsize=15, ha='right',va='top')


        count=0
        for idx,e in enumerate(output_matrix[attribute].values):
            #e is the perturbation

            left=idx
            bottom=0
            height=h
            width=bar_width

            if e!=0:
                color = color_dict[encode_pert[int(e)]]

                if ';' in encode_pert[int(e)]:
                    label=encode_pert[int(e)].replace(';','')
                else:
                    label = encode_pert[int(e)]
                    
                legend.add((color, label))
                
                count+=1
                
                rectangle=mplpatches.Rectangle((left,bottom),width,height,
                                                   facecolor=color,
                                                    edgecolor='black',alpha=1,
                                                    linewidth=0)
                
            else:
                rectangle=mplpatches.Rectangle((left,bottom),width,height,
                                                   facecolor='lightgrey',alpha=1,
                                                    edgecolor='black',
                                                    linewidth=0) 
            panel.add_patch(rectangle)



        panel.yaxis.set_label_position("right")
        panel.set_ylabel('{} ({})'.format(attribute, count), rotation = 0, fontsize=15, horizontalalignment='left')   

    legend_elements=[]
    for element in legend:
        legend_elements+=[Patch(facecolor=element[0], label=element[1])]
#     plt.subplots_adjust(wspace=2, hspace=0.1)

    plt.rcParams['legend.handlelength'] = 1.05
    plt.rcParams['legend.handleheight'] = 1.125

    legend=plt.legend(handles=list(legend_elements),bbox_to_anchor=(2.6,total_panels), fontsize=10, title="Alterations", frameon=False,labelspacing=0.1,loc='center right',
                     columnspacing=0)

    legend.get_title().set_fontsize('15')
    legend._legend_box.align = "left"

    plt.savefig(out, dpi=dpi,bbox_inches='tight')
#     plt.show()

    
def memo_oncoprint(classifier_score, output_matrix, out, dpi, cmap='tab20b'):
    output_matrix = output_matrix.T
    output_matrix = output_matrix.fillna(0)
    bait = output_matrix.columns[0]
    
    # get colors
    perturbations=set()
    for c in output_matrix.columns:

        u=(output_matrix[c].unique())
        for perturbation in u:
            if 0 != perturbation:
                perturbations.add(str(perturbation))
                
    perturbations = sorted(list(perturbations))
    perturbation_encoder = dict(zip(perturbations, np.arange(1,len(perturbations)+1)))
    
    #separate bait from rest 
    new_out = output_matrix.T.loc[output_matrix.columns[1:]].T
    #get a bait df
    bait_df = output_matrix.loc[output_matrix[output_matrix.columns[0]]!=0] #first row is always the bait
    #remove bait indeces from new_out
    new_out = new_out.drop(index=bait_df.index)
    
    # memo sort the new_output matrix
    m, gene_order = memoSort(new_out, perturbations,perturbation_encoder)
    
    bait_df = bait_df.replace(perturbation_encoder)
    bait_df = bait_df.T.loc[[bait]+list(gene_order.index)].T #updating bait_df with gene order determined by memo sort

    # sorting samples in classifer_df that match samples in bait df based on classfier scores
    bait_classifier_scores = classifier_score.loc[bait_df.index].sort_values(by='x', ascending=False)
    # Dropping bait samples, because we do not want to sort based on score
    # sorting rest classifier scores based on the sorting from memo sort
    rest_classifier_scores = classifier_score.drop(index=bait_classifier_scores.index).loc[m.index]
    # updating bait_df to reflect sorting of bait samples based on classifier scores
    bait_df= bait_df.loc[bait_classifier_scores.index]
    
    ###########---------resorting classifier scores for each catch----------###################
    sorted_rest_samples=[]
    for gene in m.columns:
        #get samples with gene perturbation
        gene_samples = m.loc[m[gene]!=0].index

        c = rest_classifier_scores.loc[gene_samples].sort_values(by='x', ascending=False)
        for samples in c.index:
            if samples not in sorted_rest_samples:
                sorted_rest_samples+=[samples]

    zero_gene_samples = list(set(m.index.values.tolist()).difference(sorted_rest_samples))
    zero_samples_sorted = rest_classifier_scores.loc[zero_gene_samples].sort_values(by='x', ascending=False).index.values.tolist()

    final_sample_sorting = sorted_rest_samples+zero_samples_sorted
    ###########---------resorting classifier scores for each catch----------###################
    
    rest_classifier_scores = rest_classifier_scores.loc[final_sample_sorting]
    m = m.loc[final_sample_sorting]  
    
    #final df describing perturbations of bait and everything after memo sort
    final_df = pd.concat([bait_df, m], axis=0, sort=False).fillna(0)
    #final classifier scores reflecting classifgier sorting for bait only samples, and memo-informed sorting for everything else
    final_scores=pd.concat([bait_classifier_scores,rest_classifier_scores], axis=0)
    
    #plotting
    bar_width=1
    total_panels=len(final_df.columns)+1
    scores=final_scores.x.values.tolist()
    encode_pert = dict()

    for k,v in perturbation_encoder.items():
        encode_pert[int(v)]=k

    if total_panels >12:
        total_panels=12

    ########----------------------Canvas--------------------##########
    width=10
    height=total_panels
    plt.figure(figsize=(width,height))  #units are always inches
    ########----------------------Canvas--------------------##########


    ########----------------------Main panel dimensions--------------------##########
    panel_width=3/width 
    panel_height=0.7/height
    classifier_height=0.85
    classifier_panel=plt.axes([0.1,classifier_height,panel_width,panel_height]) #this adds axis labels

    ########----------------------Main panel dimensions--------------------##########


    ########----------------------classifier panel--------------------##########
    for i,e in enumerate(scores):
        left=i
        bottom=0
        height=e
        width=bar_width
        rectangle=mplpatches.Rectangle((left,bottom),width,height,
                                               facecolor='black',
                                                edgecolor='black',
                                                linewidth=0) 
        classifier_panel.add_patch(rectangle)

    classifier_panel.set_xlim(0,len(scores))
    classifier_panel.set_ylim(0,1)
    classifier_panel.tick_params(axis='both', which='both',bottom=False, labelbottom=False, left=True, labelleft=True, 
            right=False, labelright=False, top=False, labeltop=False)
    classifier_panel.set_ylabel('Classifier'+'\n'+'Score', fontsize=10)

    classifier_panel.plot([0,len(scores)], [0.5, 0.5], linewidth=2, color='red', zorder=3)

    ########----------------------classifier panel--------------------##########

    c=0
    legend = set() 
    last_panel=None
    for i,attribute in enumerate(final_df.columns):
        if c!=11:
            c+=1
        else:
            break

        if last_panel is None:
            last_panel=classifier_height

        next_panel=last_panel-1.05*panel_height
        panel=plt.axes([0.1,next_panel,panel_width,panel_height], frameon=False) #this adds axis labels
        last_panel=next_panel
        h=1

        panel.set_xlim(0,len(scores))
        panel.set_ylim(0,h)
        panel.tick_params(axis='both', which='both',bottom=False, labelbottom=False, left=False, labelleft=False, 
        right=False, labelright=False, top=False, labeltop=False)

        if i==0: #always bait
            panel.text(-0.8, h/2,'Bait', rotation = 0, fontsize=15, ha='right', va='top')
        else:
            panel.text(-0.8, h/2,'Catch{}'.format(i), rotation = 0, fontsize=15, ha='right',  va='top')


        count=0
        for idx,e in enumerate(final_df[attribute].values):
            #e is the perturbation

            left=idx
            bottom=0
            height=h
            width=bar_width


            if e!=0:
                color = color_dict[encode_pert[int(e)]]

                if ';' in encode_pert[int(e)]:
                    label=encode_pert[int(e)].replace(';','')
                else:
                    label = encode_pert[int(e)]
                legend.add((color,label))

                count+=1
                rectangle=mplpatches.Rectangle((left,bottom),width,height,
                                                   facecolor=color_dict[encode_pert[int(e)]],
                                                    edgecolor='black',alpha=1,
                                                    linewidth=0)
            else:
                rectangle=mplpatches.Rectangle((left,bottom),width,height,
                                                   facecolor='lightgrey',alpha=1,
                                                    edgecolor='black',
                                                    linewidth=0) 
            panel.add_patch(rectangle)



        panel.yaxis.set_label_position("right")
        panel.set_ylabel('{} ({})'.format(attribute, count), rotation = 0, fontsize=15, horizontalalignment='left')   

    legend_elements=[]
    for element in legend:
        legend_elements+=[Patch(facecolor=element[0], label=element[1])]

    
    plt.rcParams['legend.handlelength'] = 1.05
    plt.rcParams['legend.handleheight'] = 1.125

    legend=plt.legend(handles=list(legend_elements),loc='center left', bbox_to_anchor=(2,total_panels/2 ), fontsize=10, title="Alterations", frameon=False,labelspacing=0.1,
                     columnspacing=0)


    legend.get_title().set_fontsize('15')
    legend._legend_box.align = "left"

#     plt.subplots_adjust(wspace=2, hspace=0.1)
    plt.savefig(out, dpi=dpi, bbox_inches='tight')
#     plt.show()

    
    
    
    
def main():
    ########----------------------Command line arguments--------------------##########
    parser = argparse.ArgumentParser(description="Arguments for preranked an single sample GSEA")

    parser.add_argument('-c', '--classifier_score', default=None, type=str, required=True, help='Classifier score file')
    parser.add_argument('-m', '--output_matrix', default=None, type=str, required=True, help='Output matirx file')
    parser.add_argument('-o', '--out', default='./plots/out.png', required=False,type=str, help='out file for image')
    parser.add_argument('-cmap', '--cmap', default='tab20b',type=str, required=False, help='Set colormap')

    parser.add_argument('-dpi', '--dpi', default=300,type=int, required=False, help='set dpi for image')
    parser.add_argument('-memo', '--memo', required=False, action='store_true',help='whether to memosort or not')

    
    args=parser.parse_args()
    ########----------------------Command line arguments--------------------##########
    
    
    print('parsing and sorting files...')
    classifier_score=pd.read_csv(args.classifier_score, index_col=0)
    output_matrix=pd.read_csv(args.output_matrix, index_col=0)
    file_out=args.out
    cmap=args.cmap
    dpi=int(args.dpi)
    memo=args.memo
    

    
    print('plotting...')

    if memo:
        memo_oncoprint(classifier_score, output_matrix, args.out, args.dpi, cmap)
    else:
        oncoprint(classifier_score, output_matrix, args.out, args.dpi, cmap)
    print('DONE')

if __name__ == '__main__':
    main()
