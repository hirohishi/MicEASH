import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

plt.rcParams['pdf.fonttype'] = 42

class DataPlotter:
    def __init__(self):
        #plt.rcParams["font.family"] = "Arial"
        plt.rcParams["font.size"] = 14
    
    def plot2D(self, data_file, file_name, C_MIN=None, C_MAX=None, c="seismic_r", show=False, save=False):
        os.makedirs(os.path.join(data_file, "analysis", "output_images"), exist_ok=True)
        data_file_2 = os.path.join(data_file, "analysis", "output_matrix", file_name)
        Sigma2 = np.load(data_file_2)
        cmap = plt.get_cmap(c).copy()
        cmap.set_bad("gray")
        
        if C_MIN is None:
            C_MIN = np.nanmin(Sigma2)
        if C_MAX is None:
            C_MAX = np.nanmax(Sigma2)
        
        plt.figure(figsize=(8, 8))
        plt.imshow(Sigma2, cmap=cmap, vmin=C_MIN, vmax=C_MAX)
        plt.colorbar(ticks=[C_MIN, C_MAX], shrink=0.5, orientation="vertical")
        plt.axis("off")
        plt.tight_layout()
        
        if save:
            data_file = os.path.join(data_file, "analysis", "output_images", file_name)
            output_image_path = data_file.replace('.npy', '.pdf')
            plt.savefig(output_image_path)
        if show:
            plt.show()
        plt.close()

    def plot3D_onechr(self, df_filtered, data_file, file_name, show=False, save=False):
        os.makedirs(os.path.join(data_file, "analysis", "output_images"), exist_ok=True)
        plt.rcParams["font.size"] = 14
        df_marked = df_filtered.sort_values('hyb')
        grouped = df_marked.groupby('Allele')

        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        for allele, group in grouped:
            ax.plot(group['x'], group['y'], group['z'], label=f'Allele {allele}', marker='o', linestyle='-')
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Allele distinguish for DNA spots using Kmeans Clustering')
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

        if save:
            data_file = os.path.join(data_file, "analysis", "output_images", file_name)
            output_image_path = data_file + '.pdf'
            plt.savefig(output_image_path)
        if show:
            plt.show()
        plt.close()

    def plot3D_allchr(self, df_filtered, data_file, file_name, show=False, save=False):
        os.makedirs(os.path.join(data_file, "analysis", "output_images"), exist_ok=True)
        cmap = plt.get_cmap('rainbow', 20)
        df_filtered = df_filtered.copy()
        df_filtered['chromID'] = df_filtered['chromID'].astype(int)

        fig = plt.figure(figsize=(13, 9))
        ax = fig.add_subplot(111, projection='3d')
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

        marker_dict = {0: 'o', 1: '^'} 

        # アレルによる色のリスト
        colors = []
        for chromID in sorted(df_filtered['chromID'].unique()):
            df_subset = df_filtered[df_filtered['chromID'] == chromID]
            for allele in df_subset['Allele'].unique():
                allele_subset = df_subset[df_subset['Allele'] == allele]
                marker_style = marker_dict.get(allele, 'o') 
                color = cmap(chromID-1)
                colors.append(color)
                ax.scatter(allele_subset['x'], allele_subset['y'], allele_subset['z'], 
                        color=color, 
                        marker=marker_style, 
                        s=50, 
                        alpha=0.6, 
                        linewidths=0.5, 
                        edgecolor='black')

                # ラインの追加
                sorted_allele_subset = allele_subset.sort_values('hyb')
                ax.plot(sorted_allele_subset['x'], sorted_allele_subset['y'], sorted_allele_subset['z'],
                        color=color, linewidth=2)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=1, vmax=20))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.1)
        cbar.set_label('Chromosome ID')
        cbar.set_ticks(range(1, 21)) 
        cbar.set_ticklabels(range(1, 21))

        allele_legend = ax.legend(handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='grey', alpha=1),
                                        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='grey', alpha=1)],
                                labels=['Allele 0', 'Allele 1'], loc='upper left', title="Alleles")
        plt.gca().add_artist(allele_legend)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('DNA spots for all chromosomes in a cell')
        
        if save:
            data_file = os.path.join(data_file, "analysis", "output_images", file_name)
            output_image_path = data_file + '.pdf'
            plt.savefig(output_image_path)
        if show:
            plt.show()
        plt.close()

