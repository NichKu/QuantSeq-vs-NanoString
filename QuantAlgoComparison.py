import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

df_ns = pd.read_csv("/Users/nicholas/Desktop/RNA-seq_Vifor/NS_COUNTS_NORMALIZED.csv", index_col=0)

df_seq = pd.read_csv("/Users/nicholas/Desktop/Results_quant_sf_DESeq/DESeq_output_new/DESeq2.counts_normalized.tsv"
                     , sep='\t', index_col=0)

df_seq_lex = pd.read_csv("/Users/nicholas/Desktop/RNA-seq_Vifor/Lexogen_Bluebee/DEG Analysis/"
                     "Quantification_Files/DESeq_output/DESeq2HTSseq.counts_normalized.tsv"
                     , sep='\t', index_col=0)

samples = pd.read_csv("Data/samples_annot.csv")
samples.reset_index(inplace=True)
samples.set_index("sample", inplace=True)
print(samples)

df_ns.drop(["15_post"], axis=1, inplace=True)

merged_ns_seq = pd.concat([df_ns, df_seq, df_seq_lex], axis=1, keys=['NanoString', 'RNAseq_sal', "RNAseq_lex"], join="inner").\
    swaplevel(0, 1, axis=1).sort_index(axis=1)


merged_seq = pd.concat([df_seq, df_seq_lex], axis=1, keys=['RNAseq_sal', "RNAseq_lex"], join="inner").\
    swaplevel(0, 1, axis=1).sort_index(axis=1)


count = 1

plt.subplots(figsize=(15, 10))

dict_lex = {}
dict_sal = {}

for col in merged_ns_seq.columns.levels[0]:
    slope, intercept, r_value, p_value, std_err = stats.linregress(merged_ns_seq[col]['RNAseq_lex'], merged_ns_seq[col]["RNAseq_sal"])
    Rsq = r_value ** 2
    ax = plt.subplot(3, 5, count)
    sns.regplot(x='RNAseq_sal', y="RNAseq_lex",
                        data=merged_ns_seq[col], fit_reg=False, scatter_kws={'s': 2})
    plt.axis('square')
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="black", linewidth=0.5)
    sub_merged = merged_ns_seq[col]
    kendall = sub_merged[['RNAseq_sal', "RNAseq_lex"]].corr(method="kendall").iloc[1]['RNAseq_sal']

    #ax.set(xscale="log", yscale="log")
    ax.text(0.3, 0.9, "kendall = {0:.2f}".format(kendall),
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes)

    plt.xlabel("Salmon")
    plt.ylabel("HTseq")
    plt.title(col)
    count += 1
    dict_lex[col] = r_value ** 2

plt.subplots_adjust(wspace=0.4, hspace=0.5, top=0.97, right=0.97, bottom=0.05, left=0.06)
plt.savefig("/Users/nicholas/Desktop/Salmon_vs_HTseq_NSgenes.png", dpi=500)
plt.show()


count = 1

plt.subplots(figsize=(15, 10))

for col in merged_seq.columns.levels[0]:
    slope, intercept, r_value, p_value, std_err = stats.linregress(merged_seq[col]['RNAseq_lex'], merged_seq[col]["RNAseq_sal"])
    sub_merged = merged_seq[col]
    kendall = sub_merged[['RNAseq_lex', "RNAseq_sal"]].corr(method="kendall").iloc[1]['RNAseq_lex']
    Rsq = r_value ** 2
    ax = plt.subplot(3, 5, count)
    sns.regplot(x="RNAseq_sal", y='RNAseq_lex',
                        data=merged_seq[col], fit_reg=False, scatter_kws={'s': 2})
    plt.axis('square')
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="black", linewidth=0.5)
    #ax.set(xscale="log", yscale="log")
    ax.text(0.3, 0.9, "kendall = {1:.2f}".format(r_value, kendall),
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes)
    plt.xlabel("Salmon")
    plt.ylabel("HTseq")

    plt.title(col)
    count += 1
    dict_sal[col] = r_value ** 2
    #print(col, r_value ** 2)

plt.subplots_adjust(wspace=0.4, hspace=0.5, top=0.97, right=0.97, bottom=0.05, left=0.06)
plt.savefig("/Users/nicholas/Desktop/Salmon_vs_HTseq_allgenes.png", dpi=500)
plt.show()
