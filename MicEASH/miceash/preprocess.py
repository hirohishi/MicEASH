import pandas as pd
import os

def preprocess(data_path, posi_file, anno_file, voxel_size_xyz=(103, 103, 250)):
    """
    Preprocess the DNA-seqFISH+ data.

    Parameters:
    - data_path: str. Full path for analysis.
    - posi_file: str. File name of spots coordinates from DNA-seqFISH.
    - anno_file: str. File name of annotation for DNA-seqFISH.
    - voxel_size_xyz: tuple of int, default (103, 103, 250). Voxel size for analysis.

    Returns:
    - df_seqfish: DataFrame. Preprocessed DNA-seqFISH+ data.
    """

    # Read the positions file
    df_seqfish = pd.read_csv(os.path.join(data_path, posi_file))

    #anno_file = "science.abj1966_table_s1.xlsx"

    # Read the Excel file into a dataframe
    df = pd.read_excel(os.path.join(data_path, anno_file))

    # Display the dataframe
    df['hyb'] = df.groupby('chromID')['Start'].rank(method='first') - 1
    df = df.dropna(subset=['hyb']).astype({'hyb': int})

    df_seqfish = df_seqfish.merge(df[['geneID', 'regionID', 'chromID', 'Chrom', 'Start', 'End', 'hyb']], on='geneID', how='left')
    df_seqfish['z'] = df_seqfish['z'] * voxel_size_xyz[2]/voxel_size_xyz[0]

    if 'fov' not in df_seqfish.columns:
        df_seqfish['fov'] = 0

    analysis = os.path.join(data_path, 'analysis')
    os.makedirs(analysis, exist_ok=True)
    df_seqfish.to_csv(os.path.join(analysis, posi_file.replace('.csv', '_preprocessed.csv')))

    return df_seqfish
