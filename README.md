# MicEASH

Welcome to MicEASH. This Python package is dedicated to supporting molecular biology research, particularly focusing on seqFISH analysis to explore higher-order genome structures across the entire genome. By utilizing the spatial information of genomic coordinates provided by seqFISH, MicEASH allows researchers to calculate the actual distances between genomic regions within a single cell and analyze gene transcription and protein localization related to these regions.

## Visualizing genomic distances: Insights from seqFISH data

### Execution via a Personal Computer

For this analysis, a personal computer (PC) terminal is used. We have verified that PCs with 16–32 GB of memory and 500 GB to 2 TB of storage can effectively accommodate the analysis. Compatible operating systems include Mac OS, Windows OS (using WSL2), and Linux OS (Ubuntu).

### Setting up the Execution Environment

The MicEASH package requires Jupyter Notebook and other necessary packages, which can be managed within a conda environment. To set up the analysis environment, follow these steps:

1. Download and install [Anaconda](https://www.anaconda.com/products/distribution) and [Git](https://github.com/git-guides/install-git).
2. Clone the MicEASH repository:
    ```sh
    cd ~
    mkdir [desired folder name]
    cd [desired folder name]
    git clone https://github.com/hirohishi/MicEASH.git
    ```
3. Set up the conda environment:
    ```sh
    cd MicEASH
    conda create --name mice_env python=3.11
    conda activate mice_env
    pip install -e .
    pip install notebook
    python -m ipykernel install --user --name mice_env --display-name "mice_env"
    ```

### Downloading Test Data

Download the test data and supplementary files using the following URLs:

- [DNAseqFISH+ Test Data](https://zenodo.org/records/3735329/files/DNAseqFISH+.zip?download=1)
- [Supplementary Files](https://www.science.org/doi/suppl/10.1126/science.abj1966/suppl_file/science.abj1966_tables_s1_to_s4.zip)

Unzip the downloaded files and place them in the `test_data` folder.

## Methods

### Preprocessing Bright Spot Information Data

1. Open your terminal and navigate to the MicEASH folder:
    ```sh
    cd [desired folder name]/MicEASH/jupyter
    ```
2. Activate the conda environment:
    ```sh
    conda activate mice_env
    ```
3. Launch Jupyter Notebook:
    ```sh
    jupyter notebook
    ```
4. Open the `Spatial_distance_analysis.ipynb` file in the Jupyter Notebook interface.
5. Execute the notebook code to preprocess bright spot information data.

## Acknowledgements

We thank Hiroshi Ochiai (Kyushu University) and Soya Shinkai (RIKEN BDR) for their invaluable comments on the methods. This work was supported by JSPS KAKENHI to H.Oh (JP22K15084).


## Contact

If you have any questions or need support, please contact our team at [hirohishi@outlook.jp](mailto:hirohishi@outlook.jp).
