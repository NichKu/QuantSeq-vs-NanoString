import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


df_seq = pd.read_csv("Data/star_alignment_plot.tsv", sep='\t', index_col=0)
df_info = pd.read_csv("Data/Corr_variables.txt", sep='\t', index_col=0)
df_MQC = pd.read_csv("Data/MultiQC.txt", sep='\t', index_col=0)

#df_seq["sum"] = df_seq.sum(axis=1)
#df_seq["perc_multi"] = df_seq["Mapped to multiple loci"] / df_seq["sum"] * 100


merged = pd.merge(df_MQC, df_info, left_on="Sample Name", right_index=True)



columns_perc = ['dupInt', 'Error rate', '% Mapped', '% Aligned', '% Dups', '% GC', '% BP Trimmed', '% Dups.1', '% GC.1']

columns_bp = ['Length.1','Length']

for col in columns_perc:
    merged[col] = merged[col].str.rstrip('%').astype(float) / 100

for col in columns_bp:
    merged[col] = merged[col].str.rstrip(' bp').astype(int)


#merged.to_csv("./Correlation_Matrix ")

merged.drop(columns=["dupInt"], inplace=True)

def heatMap(df):
    #Create Correlation df
    corr = df.corr(method="spearman")

    #Plot figsize
    fig, ax = plt.subplots(figsize=(11, 8))

    #Generate Color Map
    colormap = sns.diverging_palette(220, 10, as_cmap=True)

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=bool))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    #Generate Heat Map, allow annotations and place floats in map
    sns.heatmap(corr, mask=mask, cmap=colormap, center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5}, annot=True, fmt=".2f", annot_kws={"size": 5})

    #Apply xticks
    plt.xticks(range(len(corr.columns)), corr.columns);

    #Apply yticks
    plt.yticks(range(len(corr.columns)), corr.columns)
    plt.tight_layout()

    plt.savefig("/Users/nicholas/Desktop/CorrMatrixQCparam.png", dpi=500)
    #show plot
    plt.show()
heatMap(merged)



#sns.regplot(merged["DV200"], merged["perc_multi"])
#plt.show()