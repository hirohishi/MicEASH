import pandas as pd
import os
from scipy.spatial import distance
import seaborn as sns
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
import click
from sklearn.cluster import SpectralClustering, AgglomerativeClustering
from sklearn.mixture import GaussianMixture

sns.set(font_scale=1.8)
sns.set_style("whitegrid", {'grid.linestyle': '--'})

@click.group()
def cli():
    """ This is a command line tool for DNA-seqFISH+ analysis."""
    pass

@cli.command()
@click.option("--data_path", "-p", help="Full path for analysis", type=str, required=True)
@click.option("--posi_file","-f1", help="File name of spots coordinates from DNA-seqFISH", type=str, required=True)
@click.option("--anno_file", "-f2", help="File name of annoation for DNA-seqFISH", type=str, required=True)
@click.option("--voxel_size_xyz","-vs", default=(103, 103, 250), help="Voxel size for analysis, default is 103, 103, 250. Put x y z values like -vs 103 103 250.", type=(int, int, int), required=False)
def preprocess(data_path, posi_file, anno_file, voxel_size_xyz):
    # path of parental folder
    #data_path = "/nas/data/Microscope/Ohishi-2/analysis/240110-seqFISH_Takei_Nature/DNAseqFISH+"
    #posi_file = "DNAseqFISH+1Mbloci-E14-replicate2.csv"

    df_seqfish = pd.read_csv(os.path.join(data_path, posi_file), index_col=0)

    #anno_file = "science.abj1966_table_s1.xlsx"

    # Read the Excel file into a dataframe
    df = pd.read_excel(os.path.join(data_path, anno_file))

    # Display the dataframe
    df['hyb'] = df.groupby('chromID')['Start'].rank(method='first') - 1
    df = df.dropna(subset=['hyb']).astype({'hyb': int})

    df_seqfish = df_seqfish.merge(df[['geneID', 'regionID', 'chromID', 'Chrom', 'Start', 'End', 'hyb']], on='geneID', how='left')
    df_seqfish['z'] = df_seqfish['z'] * voxel_size_xyz[2]/voxel_size_xyz[0]

    # data_pathにanlysisという名前のフォルダを作成
    analysis = os.path.join(data_path, 'analysis')
    os.makedirs(analysis, exist_ok=True)
    df_seqfish.to_csv(os.path.join(analysis, posi_file.replace('.csv', '_preprocessed.csv')))

    return df_seqfish


@cli.command()
@click.option("--data_path", "-p", help="Full path for analysis", type=str, required=True)
@click.option("--posi_file","-f", help="File name of spots coordinates from DNA-seqFISH", type=str, required=True)
@click.option("--clustering_method_name","-cm", default='ward', help="Select methods for clustering to split two allelles. You can choose ward, gmm (GaussianMixutre)", type=str, required=False)
def align(data_path, posi_file, clustering_method_name):
    analysis = os.path.join(data_path, 'analysis')
    df_seqfish = pd.read_csv(os.path.join(analysis, posi_file), index_col=0)

    df_marked_list = []
    for fov_value in df_seqfish['fov'].unique():
        df_fov = df_seqfish[df_seqfish['fov'] == fov_value]
        for cellID_value in df_fov['cellID'].unique():
            for chr in df_seqfish['chromID'].unique():
            #chr = 12
                df_marked = spot_align_each(df_seqfish, cellID=cellID_value, fov=fov_value, chromosome=chr, clustering_method=clustering_method_name)
                df_marked_list.append(df_marked)

    # Concatenate the list of DataFrames into a single DataFrame
    df_marked = pd.concat(df_marked_list, ignore_index=True)
    #df_marked['finalcellID'] = df_marked.groupby(['cellID', 'fov', 'Allele']).ngroup()
    df_marked.to_csv(os.path.join(analysis, posi_file.replace('.csv', '_marked.csv')))

    return


def spot_align_each(df, cellID, fov, chromosome=1, clustering_method='ward'):
    df_sample = df[(df['cellID'] == cellID) & (df['fov'] == fov) & (df['chromID'] == chromosome)].copy()
    data = df_sample[['x', 'y', 'z']].values

    # Apply Spectral Clustering, Ward method, or Gaussian Mixture Model based on the parameter
    if clustering_method == 'spectral':
        clustering_model = SpectralClustering(n_clusters=2, eigen_solver='arpack',affinity='nearest_neighbors',random_state=1)
        # Fit the model and predict labels
        labels = clustering_model.fit_predict(data)
    elif clustering_method == 'ward':
        clustering_model = AgglomerativeClustering(n_clusters=2, linkage='ward')
        # Fit the model and predict labels
        labels = clustering_model.fit_predict(data)
    elif clustering_method == 'gmm':
        clustering_model = GaussianMixture(n_components=2)
        # Fit the model
        clustering_model.fit(data)
        # Predict labels
        labels = clustering_model.predict(data)
    else:
        raise ValueError("Invalid clustering method. Choose 'ward', or 'gmm'.")

    # Assign labels to the original dataframe
    df_sample['Allele'] = labels
    df_sample['marked'] = 1
    df_marked = df_sample.copy()
    df_marked.reset_index(drop=True, inplace=True)
    
    return df_marked


@cli.command()
@click.option("--data_path", "-p", help="Full path for analysis", type=str, required=True)
@click.option("--posi_file","-f", help="File name of spots coordinates from DNA-seqFISH", type=str, required=True)
@click.option("--chromosome_number","-cn", default=12, help="Define number of chromosome", type=int, required=False)
@click.option("--voxel_size_xyz","-vs", default=(103, 103, 250), help="Voxel size for analysis, default is 103, 103, 250. Put x y z values like -vs 103 103 250.", type=(int, int, int), required=False)
def matrix(data_path, posi_file, chromosome_number, voxel_size_xyz):

    analysis = os.path.join(data_path, 'analysis')
    os.makedirs(analysis, exist_ok=True)
    df_dna_new = pd.read_csv(os.path.join(analysis, posi_file), index_col=0)
    df_dna_new = df_dna_new[df_dna_new['chromID'] == chromosome_number]
    df_dna_new['finalcellID'] = df_dna_new.groupby(['cellID', 'fov', 'Allele']).ngroup()
    data_IDs, data_XYZ = data_in(df_dna_new, voxel_size_xyz)
    distance = dist(data_IDs, data_XYZ)
    mean = np_mean(distance,data_IDs)
    medi = np_median(distance, data_IDs)
    analysis = os.path.join(data_path, 'analysis', 'output_matrix')
    os.makedirs(analysis, exist_ok=True)
    # Save the distance matrix
    np.save(os.path.join(analysis, posi_file.replace('.csv', f'_chr{chromosome_number}_distance.npy')), distance)
    np.save(os.path.join(analysis, posi_file.replace('.csv', f'_chr{chromosome_number}_mean.npy')), mean)
    np.save(os.path.join(analysis, posi_file.replace('.csv', f'_chr{chromosome_number}_medi.npy')), medi)

    return distance, mean, medi

def data_in(df, voxel_size_xyz):

    pix_xy = voxel_size_xyz[0]  # defoult 103 nm/px
    #pix_z = voxel_size_xyz[2]  # defoult 250 nm/px
    df = df.copy()  # Make a copy to avoid modifying the original dataframe
    df['hyb'] = df['hyb'].astype(int)
    df[["x", "y", "z"]] *= pix_xy
    #df["z"] *= pix_z

    return df[["finalcellID", "hyb"]].values, df[["x", "y", "z"]].values

def dist(data_IDs, data_XYZ):

    SAMPLE_num = len(np.unique(data_IDs[:, 0])) + 1
    N = np.amax(data_IDs[:, 1]) + 1  # Number of hybs
    distance = np.full((N, N, SAMPLE_num), np.nan)  # Initialize the distance array

    for num, sample in enumerate(np.unique(data_IDs[:, 0])):
        indices = np.where(data_IDs[:, 0] == sample)[0]
        if len(indices) > 1:
            sample_XYZ = data_XYZ[indices]
            diff = sample_XYZ[:, np.newaxis, :] - sample_XYZ[np.newaxis, :, :]
            dists = np.sqrt(np.sum(diff**2, axis=-1))
            n, m = np.meshgrid(data_IDs[indices, 1], data_IDs[indices, 1])
            distance[n, m, num] = dists

    return distance

def np_median(distance, data_IDs):

    N = np.amax(data_IDs[:, 1]) + 1
    Sigma2 = np.zeros((N, N))
    for n in range(N):
        for m in range(n + 1, N):
            Sigma2[n, m] = Sigma2[m, n] = np.nanmedian(distance[n, m, :])

    return Sigma2

def np_mean(distance, data_IDs):

    N = np.amax(data_IDs[:, 1]) + 1
    Sigma2 = np.zeros((N, N))
    for n in range(N):
        for m in range(n + 1, N):
            Sigma2[n, m] = Sigma2[m, n] = np.nanmean(distance[n, m, :])

    return Sigma2

@cli.command()
@click.option("--data_path", "-p", help="Full path for analysis", type=str, required=True)
@click.option("--data_file","-f", help="File name of spots coordinates from DNA-seqFISH", type=str, required=True)
def plot(data_path, data_file, C_MIN=None, C_MAX=None, c = "seismic_r"):
    
    analysis = os.path.join(data_path, 'analysis', 'output_matrix')
    Sigma2 = np.load(os.path.join(analysis, data_file))
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 36
    cmap = plt.colormaps.get_cmap(f"{c}").copy()
    cmap.set_bad("gray")
    #C_MAX = 0.2
    if(C_MIN == None):
        C_MIN = np.nanmin(Sigma2)
    
    if(C_MAX == None):
        C_MAX = np.nanmax(Sigma2)

    analysis = os.path.join(data_path, 'analysis', 'output_images')
    os.makedirs(analysis, exist_ok=True)

    C = Sigma2
    plt.figure(figsize=(10, 10))
    plt.imshow(C, cmap=cmap, vmin=C_MIN, vmax=C_MAX)
    plt.colorbar(ticks=[C_MIN, C_MAX], shrink=0.5, orientation="vertical")
                # label="Normalized contact probability")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(analysis, data_file.replace('.npy', '.pdf')))
    plt.close()

    return

if __name__ == "__main__": 
    cli()