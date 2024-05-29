import pandas as pd
import numpy as np
import os

def matrix(data_path, posi_file, chromosome, voxel_size_xyz):

    analysis = os.path.join(data_path, 'analysis')
    df_dna_new = pd.read_csv(os.path.join(analysis, posi_file), index_col=0)
    df_dna_new = df_dna_new[df_dna_new['Chrom'] == chromosome]
    df_dna_new['finalcellID'] = df_dna_new.groupby(['cellID', 'fov', 'Allele']).ngroup()
    data_IDs, data_XYZ = data_in(df_dna_new, voxel_size_xyz)
    distance = dist(data_IDs, data_XYZ)
    mean = np_mean(distance,data_IDs)
    medi = np_median(distance, data_IDs)
    analysis = os.path.join(data_path, 'analysis', 'output_matrix')
    os.makedirs(analysis, exist_ok=True)
    # Save the distance matrix
    np.save(os.path.join(analysis, posi_file.replace('.csv', f'_{chromosome}_distance.npy')), distance)
    np.save(os.path.join(analysis, posi_file.replace('.csv', f'_{chromosome}_mean.npy')), mean)
    np.save(os.path.join(analysis, posi_file.replace('.csv', f'_{chromosome}_medi.npy')), medi)

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
