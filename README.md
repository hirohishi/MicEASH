Based on the detailed information you provided, I can create an English translation and format it in the style of a README.md. Let's start with translating and organizing this text.

---

# MicEASH

Welcome to MicEASH. This Python package is dedicated to supporting molecular biology research, particularly focusing on seqFISH analysis to explore higher-order genome structures across the entire genome. By utilizing the spatial information of genomic coordinates provided by seqFISH, MicEASH allows researchers to calculate the actual distances between genomic regions within a single cell and analyze gene transcription and protein localization related to these regions.

## Features

- **Data Preprocessing:** Prepare your DNA-seqFISH data for analysis by combining target region information with coordinate data, and adjusting Z coordinates based on voxel size.
- **Data Alignment:** Use preprocessed DNA-seqFISH data files for allelic differentiation across homologous chromosomes, employing clustering methods such as KMeans and spectral clustering.
- **Visualization of Nuclear Bright Spots:** Visualize chromosome 12 and other chromosomes to verify allelic distinction visually.
- **Distance Calculation:** Compute Euclidean distances between genomic regions within the same allele, and calculate median and average distance matrices stored as Numpy arrays.
- **Visualization of Distance Matrices:** Heatmaps display the features of higher-order genomic structures, highlighting regions with frequent contacts.
- **Genomic Distance Analysis:** Compare the physical distances within the nucleus to the genomic distances, providing insight into the folded structure of the genome.

## Installation

Install MicEASH using the following steps:

1. **Download MicEASH Package:**
   ```bash
   git clone https://github.com/hirohishi/MicEASH.git
   ```

2. **Set Up Environment:**
   Using Anaconda, manage packages required for analysis easily.
   ```bash
   cd MicEASH
   conda create --name mice_env python=3.11
   conda activate mice_env
   pip install hic-straw
   pip install -e .
   ```

3. **Download Test Data:**
   Download and unzip test data for seqFISH analysis.
   ```bash
   mkdir test_data
   cd test_data/
   unzip DNAseqFISH+.zip
   unzip science.abj1966_tables_s1_to_s4.zip
   ```

## Usage

Here is how to preprocess data, align DNA sequences, and visualize results with MicEASH:

1. **Preprocess Data:**
   Load the specified bright spot coordinate file and combine it with target region information.

2. **Data Alignment:**
   Classify bright spots into clusters, distinguishing alleles for each chromosome.

3. **Visualize and Analyze:**
   Generate heatmaps and calculate distances, using scripts provided in Jupyter notebooks, such as `240110-seqFISH_analysis_10.ipynb`.

## Documentation

For detailed documentation on using MicEASH, visit [MicEASH Documentation](https://example.com/miceash-docs).

## Contributing

Contributions are welcome! Please review our [contributing guidelines](https://example.com/miceash-contributing) before making a contribution.

## License

MicEASH is released under the MIT license. See [LICENSE](LICENSE.md) for more details.

## Contact

If you have any questions or need support, please contact our team at [your-email@example.com](mailto:your-email@example.com).

---

Would you like any changes or additional information included in this README.md file?
