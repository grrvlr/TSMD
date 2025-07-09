import json
import os 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#settings for the plots
methods_dict={'MatrixProfile': 'STOMP', 'PanMatrixProfile':'PanMP', 'LatentMotif':'LatentMotif','MDL':'MDL-Clust','Motiflets':'Motiflets','Motiflets_numba': 'Motiflets','BasePersistentPattern':'PEPA','Valmod':'Valmod','Baseline':'SetFinder','AdaptativeBasePersistentPattern':'A-PEPA','Grammarviz':'Grammarviz','LocoMotif':'LoCoMotif'}
color_palette={'MatrixProfile':'darkviolet','PanMatrixProfile':'mediumorchid','LocoMotif':'lightskyblue','LatentMotif':'darkorange','MDL':'cornflowerblue','Motiflets':'crimson','BasePersistentPattern':'deeppink','Valmod':'violet','Baseline':'orange','AdaptativeBasePersistentPattern':'hotpink','Grammarviz':'royalblue','Motiflets_numba':'crimson'}
marker_dict={'MatrixProfile':'H','PanMatrixProfile':'h','LocoMotif':'D','LatentMotif':'P','MDL':'^','Motiflets':'x','BasePersistentPattern':'s','Valmod':'p','Baseline':'+','AdaptativeBasePersistentPattern':'s','Grammarviz':'v' , 'Motiflets_numba':'x'}
size_dict={'MatrixProfile':10,'PanMatrixProfile':10,'LocoMotif':10,'LatentMotif':17,'MDL':10,'Motiflets':20,'BasePersistentPattern':10,'Valmod':10,'Baseline':23,'AdaptativeBasePersistentPattern':10,'Grammarviz':10, 'Motiflets_numba':20}
metrics_dict={'es_fscore_0.5_mean' : 'Accuracy (fscore) average','es_fscore_0.5_std' : 'Accuracy (fscore) standard deviation ','es_precision_0.5_mean' : 'Accuracy (precision) average','es_precision_0.5_std' : 'Accuracy (precision) standard deviation','es_recall_0.5_mean' : 'Accuracy (recall) average','es_recall_0.5_std' : 'Accuracy (recall) standard deviation',  'execution_time_mean':'execution time average ','execution_time_std':'execution time standard deviation',
              'es_fscore_0.25_mean' : 'Accuracy (fscore) average','es_fscore_0.25_std' : 'Accuracy (fscore) standard deviation ','es_precision_0.25_mean' : 'Accuracy (precision) average','es_precision_0.25_std' : 'Accuracy (precision) standard deviation','es_recall_0.25_mean' : 'Accuracy (recall) average','es_recall_0.25_std' : 'Accuracy (recall) standard deviation',
              'es_fscore_0.75_mean' : 'Accuracy (fscore) average','es_fscore_0.75_std' : 'Accuracy (fscore) standard deviation ','es_precision_0.75_mean' : 'Accuracy (precision) average','es_precision_0.75_std' : 'Accuracy (precision) standard deviation','es_recall_0.75_mean' : 'Accuracy (recall) average','es_recall_0.75_std' : 'Accuracy (recall) standard deviation'}

def extract_scores(results_path,algo,alpha=0.5):
    
    infos_path=results_path+'/Infos/'+str(algo)
    metrics_path=results_path+'/Metrics/'+str(algo)

    es_precision_list=[]
    es_recall_list=[]
    es_fscore_list=[]

    metrics_files=sorted(os.listdir(metrics_path))
    
    for file in metrics_files:
        df = pd.read_csv(metrics_path+'/'+file)
        es_precision_list.append(df.loc[df['metric']=='es-precision_'+str(alpha),'score'].iloc[0])
        es_recall_list.append(df.loc[df['metric']=='es-recall_'+str(alpha),'score'].iloc[0])
        es_fscore_list.append(df.loc[df['metric']=='es-fscore_'+str(alpha),'score'].iloc[0])
    
    execution_times=[]
    infos_files=sorted(os.listdir(infos_path))

    for file in infos_files:
        if file != '.DS_Store':
            infos_file = open(infos_path+'/'+file)
            infos=json.load(infos_file)
            execution_times.append(infos['execution_time'])

    return es_precision_list,es_recall_list,es_fscore_list,execution_times

def build_results_csv(general_path,write_results_path,dataset_list,algos_list,alpha=0.5):
    for dataset in dataset_list:
        df=pd.DataFrame()
        df['metrics']=['es_fscore_'+str(alpha)+'_mean','es_fscore_'+str(alpha)+'_std','es_precision_'+str(alpha)+'_mean','es_precision_'+str(alpha)+'_std','es_recall_'+str(alpha)+'_mean','es_recall_'+str(alpha)+'_std','execution_time_mean','execution_time_std']
        dataset_path=general_path+dataset+'/'
        #algos_list=[algo for algo in os.listdir(dataset_path+'Results/Metrics') if algo != '.DS_Store']
        for algo in algos_list:
            try:
                es_precision_list,es_recall_list,es_fscore_list,execution_times=extract_scores(dataset_path,algo,alpha)
                metrics_list=[round(np.mean(es_fscore_list),3),round(np.std(es_fscore_list),3),
                            round(np.mean(es_precision_list),3),round(np.std(es_precision_list),3)
                            ,round(np.mean(es_recall_list),3),round(np.std(es_recall_list),3),
                            round(np.mean(execution_times),3),round(np.std(execution_times),3)]
            except:
                metrics_list=[0]*8
            df[algo]=metrics_list
        if not os.path.exists(write_results_path):
            os.makedirs(write_results_path)
        #dataset_=dataset.replace('/','_')
        df.to_csv(write_results_path+ dataset+'_'+str(alpha)+'.csv')
    return df


def get_results_by_param_by_algo(path,param_list,param_name,alpha=0.5):
    results_dict={}
    #initialization of the dictionnary with the first element of the parameter list 
    param=param_list[0]
    file_name=path+param_name +'_'+str(param)+'_'+str(alpha)+'.csv'
    df=pd.read_csv(file_name,index_col=False)
    columns=[column for column in df if column!='metrics'and column!='Unnamed: 0']
    for column in columns:
        for metric in df['metrics']:
            results_dict[column+'_'+metric]=[df.loc[df['metrics'] == metric, column].iloc[0]]

    for param in param_list[1:]:
        file_name=path+param_name +'_'+str(param)+'_'+str(alpha)+'.csv'
        df=pd.read_csv(file_name,index_col=False)
        for column in columns:
            for metric in df['metrics']:
                results_dict[column+'_'+metric].append(df.loc[df['metrics'] == metric, column].iloc[0])
    return results_dict, df['metrics'] 


def plot_results_by_param_by_algo(results_dict, param_list, xlabel, metrics,alpha=0.5):
    fig, axes = plt.subplots(len(metrics), 1, figsize=(8, 5 * len(metrics)))
    fig.suptitle("Results for alpha="+str(alpha), fontsize=16, y=0.95)
    
    if len(metrics) == 1:
        axes = [axes]  # S'assurer que axes est toujours une liste
    
    for ax, metric in zip(axes, metrics):
        for key in results_dict.keys():
            if metric in key:
                algo = key.split('_' + metric)[0]
                ax.plot(param_list, results_dict[key], label=methods_dict[algo],
                        color=color_palette[algo], marker=marker_dict[algo],
                        markersize=size_dict[algo] / 1.7)
        
        if metric == 'execution_time_mean':
            ax.set_yscale('log')
            ax.set_xscale('log')
        
        ax.set_ylabel(metrics_dict[metric])
        ax.set_xlabel(xlabel)
        ax.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.2))
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

    
import matplotlib.pyplot as plt
import numpy as np

def get_results_2d(path,param_list,param_name,alpha=0.5):
    n=len(param_list)
    results_dict={}
    #initialization of the dictionnary with the first element
    param1,param2=param_list[0],param_list[0]
    file_name=path+param_name+'_'+str(param1)+'and'+str(param2)+'_'+str(alpha)+'.csv'
    df=pd.read_csv(file_name,index_col=False)
    columns=[column for column in df if column!='metrics'and column!='Unnamed: 0']
    for column in columns:
        for metric in df['metrics']:
            results_dict[column+'_'+metric]=np.zeros((n,n))
            results_dict[column+'_'+metric][0,0]=df.loc[df['metrics'] == metric, column].iloc[0]

    for i,param1 in enumerate(param_list):
        for j,param2 in enumerate(param_list):
            if param1<=param2 and (param1!=param_list[0] or param2!=param_list[0]):
                file_name=path+param_name+'_'+str(param1)+'and'+str(param2)+'_'+str(alpha)+'.csv'
                df=pd.read_csv(file_name,index_col=False)
                for column in columns:
                    for metric in df['metrics']:
                        results_dict[column+'_'+metric][i,j]=df.loc[df['metrics'] == metric, column].iloc[0]
                        results_dict[column+'_'+metric][j,i]=df.loc[df['metrics'] == metric, column].iloc[0]
    return results_dict, df['metrics'] 

def plot_results_2d(results_dict, param_list, labels,  metric):
    # Filtrer les algorithmes qui contiennent le metric
    filtered_keys = [key for key in results_dict.keys() if metric in key]

    # Créer une figure avec autant de sous-graphes que d'algorithmes
    n_plots = len(filtered_keys)
    fig, axes = plt.subplots(n_plots, 1, figsize=(10, n_plots * 6), squeeze=False)

    # Stocker les valeurs min et max pour homogénéiser l'échelle de la colorbar
    vmin = min(np.min(results_dict[key]) for key in filtered_keys)
    vmax = max(np.max(results_dict[key]) for key in filtered_keys)

    # Boucler sur les résultats pour tracer chaque heatmap
    for i, key in enumerate(filtered_keys):
        algo = key.split('_' + metric)[0]
        results = results_dict[key]
        
        # Tracer la heatmap dans le sous-graphique correspondant avec vmin et vmax
        im = axes[i, 0].imshow(results, cmap='hot', interpolation='nearest', vmin=vmin, vmax=vmax)

        # Ajouter les labels pour les axes
        axes[i, 0].set_xticks(range(len(param_list)))
        axes[i, 0].set_yticks(range(len(param_list)))
        axes[i, 0].set_xticklabels(param_list)
        axes[i, 0].set_yticklabels(param_list)
        axes[i, 0].set_title(f"fscore heatmap for {methods_dict[algo]}")
        axes[i, 0].set_xlabel(labels + " 2")
        axes[i, 0].set_ylabel(labels+ " 1")

    # Ajouter une colorbar commune, mais à côté de chaque heatmap
    fig.subplots_adjust(right=0.85)  # Laisser de la place à droite pour la colorbar
    cbar_ax = fig.add_axes([0.88, 0.1, 0.02, 0.08])  # Positionner la colorbar à droite de toutes les heatmaps
    fig.colorbar(im, cax=cbar_ax)

    # Ajuster l'espacement entre les sous-plots
    #plt.tight_layout()
    plt.show()


def build_exec_time_csv(general_path,write_results_path,dataset_list,algos_list):
    for dataset in dataset_list:
        df=pd.DataFrame()
        df['metrics']=['execution_time_mean']
        dataset_path=general_path+dataset+'/'
        for algo in algos_list:
            try:
                infos_path=dataset_path+'Infos/'+str(algo)
                execution_times=[]
                infos_files=sorted(os.listdir(infos_path))

                for file in infos_files:
                    if file != '.DS_Store':
                        infos_file = open(infos_path+'/'+file)
                        infos=json.load(infos_file)
                        execution_times.append(infos['execution_time'])
                metrics_list=[np.round(np.mean(execution_times),3)]
            except:
                metrics_list=[0]
            df[algo]=metrics_list
        if not os.path.exists(write_results_path):
            os.makedirs(write_results_path)
        dataset_=dataset.replace('/','_')
        df.to_csv(write_results_path+ dataset_+'.csv')
    return df

def get_results_fixed_motifs_var_length_by_algo(path,exact_ts_length_list):
    results_dict={}
    exact_ts_length=exact_ts_length_list[0]
    file_name=path+'exact_ts_length_'+str(exact_ts_length)+'.csv'
    df=pd.read_csv(file_name,index_col=False)
    columns=[column for column in df if column!='metrics'and column!='Unnamed: 0']
    for column in columns:
        for metric in df['metrics']:
            results_dict[column+'_'+metric]=[df.loc[df['metrics'] == metric, column].iloc[0]]

    for exact_ts_length in exact_ts_length_list[1:]:
        file_name=path+'exact_ts_length_'+str(exact_ts_length)+'.csv'
        df=pd.read_csv(file_name,index_col=False)
        for column in columns:
            for metric in df['metrics']:
                results_dict[column+'_'+metric].append(df.loc[df['metrics'] == metric, column].iloc[0])
    return results_dict, df['metrics'] 

def plot_exec_time_fixed_motifs_var_length(results_dict,exact_ts_length_list):
    for key in results_dict.keys():
        if 'execution_time_mean' in key:
            algo=key.split('_'+'execution_time_mean')[0]
            filtered_values=[x for x in results_dict[key] if (x>0 and x<=20000)]
            n=len(filtered_values)
            plt.plot(exact_ts_length_list[:n],filtered_values,label=methods_dict[algo],color=color_palette[algo],marker=marker_dict[algo],markersize=size_dict[algo]/1.7)
        
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Execution time (in seconds)')
    plt.axhline(y = 20000, color = 'r', linestyle = 'dashed') 

    plt.xlabel('Length of the time series')
    plt.legend(loc='upper center',ncol=4,bbox_to_anchor=(0.47,1.22))
    plt.show()


def build_fscore_by_ts(datasets,algos_list):
    big_df = pd.DataFrame(columns=['file_name']+algos_list)
    rq1_path=os.getcwd()+'/RQ1/'
    data_path=os.getcwd()+'/../data/processed_data/'
    for dataset in datasets:
        dataset_path=data_path+dataset+'/Data'
        n_ts=len(os.listdir(dataset_path))
        name_ex=os.listdir(dataset_path)[0].split('_')[0]
        results_path=rq1_path+dataset+'/Metrics/'
        for i in range(n_ts):
            file_name=name_ex+'_'+str(i) 
            current_line_df={'file_name':file_name}
            for algo in algos_list:
                specific_file=results_path+algo+'/' + algo+'_'+file_name+'.csv'
                if os.path.exists(specific_file):
                    df=pd.read_csv(results_path+algo+'/' + algo+'_'+file_name+'.csv')
                    current_line_df[algo]=round(df.loc[df['metric']=='es-fscore_0.5','score'].iloc[0],2)
                elif algo=='Motiflets' and os.path.exists(results_path+'Motiflets_numba'+'/' + 'Motiflets_numba_'+file_name+'.csv') :
                    df=pd.read_csv(results_path+'Motiflets_numba'+'/' + 'Motiflets_numba_'+file_name+'.csv')
                    current_line_df[algo]=round(df.loc[df['metric']=='es-fscore_0.5','score'].iloc[0],2)
                else:
                    current_line_df[algo]=0.00
            big_df.loc[len(big_df)]=current_line_df
    big_df.to_csv(rq1_path+'results_by_ts.csv')

