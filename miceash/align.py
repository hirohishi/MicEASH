import pandas as pd
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
import numpy as np
import os

def align(data_path, posi_file, clustering_method_name, threshold_spots_value=20):
    analysis = os.path.join(data_path, 'analysis')
    df_seqfish = pd.read_csv(os.path.join(analysis, posi_file), index_col=0)

    df_marked_list = []
    for fov_value in df_seqfish['fov'].unique():
        df_fov = df_seqfish[df_seqfish['fov'] == fov_value]
        for cellID_value in df_fov['cellID'].unique():
            for chr in df_seqfish['Chrom'].unique():
            #chr = "chr12"
                df_marked = spot_align_each(df_seqfish, cellID=cellID_value, fov=fov_value, chromosome=chr, clustering_method=clustering_method_name, threshold_spots=threshold_spots_value)
                df_marked_list.append(df_marked)

    # Concatenate the list of DataFrames into a single DataFrame
    df_marked = pd.concat(df_marked_list, ignore_index=True)
    #df_marked['finalcellID'] = df_marked.groupby(['cellID', 'fov', 'Allele']).ngroup()
    df_marked.to_csv(os.path.join(analysis, posi_file.replace('.csv', '_marked.csv')))

    return df_marked


def spot_align_each(df, cellID, fov, chromosome="chr1", clustering_method='ward', threshold_spots=20, sex = "male"):
    df_sample = df[(df['cellID'] == cellID) & (df['fov'] == fov) & (df['Chrom'] == chromosome)].copy()
    data = df_sample[['x', 'y', 'z']].values

    if len(data) < threshold_spots or (sex == "male" and chromosome == "chrX"):
        df_sample['Allele'] = 0
        df_sample['marked'] = 0
    else:
        if clustering_method == 'spectral':
            clustering_model = SpectralClustering(n_clusters=2, eigen_solver='arpack',affinity='nearest_neighbors',random_state=1)
            labels = clustering_model.fit_predict(data)
        elif clustering_method == 'ward':
            clustering_model = AgglomerativeClustering(n_clusters=2, linkage='ward')
            labels = clustering_model.fit_predict(data)
        elif clustering_method == 'gmm':
            clustering_model = GaussianMixture(n_components=2)
            clustering_model.fit(data)
            labels = clustering_model.predict(data)
        elif clustering_method == 'kmeans':
            clustering_model = KMeans(n_clusters=2, random_state=1, n_init=10)
            labels = clustering_model.fit_predict(data)
        else:
            raise ValueError("Invalid clustering method. Choose 'ward', 'spectral', 'kmeans' or 'gmm'.")
        
        # Assign labels to the original dataframe
        df_sample['Allele'] = labels
        df_sample['marked'] = 1

    df_marked = df_sample.copy()
    df_marked.reset_index(drop=True, inplace=True)
    
    return df_marked
