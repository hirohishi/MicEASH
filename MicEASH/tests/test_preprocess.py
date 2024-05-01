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
    df_seqfish = pd.read_csv(os.path.join(data_path, posi_file), index_col=0)

    # Read the annotation Excel file
    df_anno = pd.read_excel(os.path.join(data_path, anno_file))
    df_anno['hyb'] = df_anno.groupby('chromID')['Start'].rank(method='first') - 1
    df_anno = df_anno.dropna(subset=['hyb']).astype({'hyb': int})

    # Merge the position and annotation data
    df_seqfish = df_seqfish.merge(df_anno[['geneID', 'regionID', 'chromID', 'Chrom', 'Start', 'End', 'hyb']], 
                                  on='geneID', how='left')
    df_seqfish['z'] = df_seqfish['z'] * voxel_size_xyz[2] / voxel_size_xyz[0]

    # Create an analysis folder if it doesn't exist
    analysis_path = os.path.join(data_path, 'analysis')
    os.makedirs(analysis_path, exist_ok=True)

    # Save the preprocessed data to a new CSV file
    output_file = posi_file.replace('.csv', '_preprocessed.csv')
    df_seqfish.to_csv(os.path.join(analysis_path, output_file))

    return df_seqfish
