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

df_ns.drop(["15_post"], axis=1, inplace=True)

merged_ns_seq = pd.concat([df_ns, df_seq, df_seq_lex], axis=1, keys=['NanoString', 'RNAseq_sal', "RNAseq_lex"], join="inner").\
    swaplevel(0, 1, axis=1).sort_index(axis=1)
merged_ns_lex = pd.concat([df_ns, df_seq_lex], axis=1, keys=['NanoString','RNAseq'], join="inner").\
    swaplevel(0, 1, axis=1).sort_index(axis=1)

count = 1



dict_lex = {}
dict_sal = {}

concat_lex = pd.DataFrame(columns=["NanoString", "RNAseq_lex"])
for col in merged_ns_seq.columns.levels[0]:
    slope, intercept, r_value, p_value, std_err = stats.linregress(merged_ns_seq[col]["NanoString"],
                                                                   merged_ns_seq[col]["RNAseq_lex"])
    Rsq = r_value ** 2

    sns.regplot(x="NanoString", y="RNAseq_lex",
                        data=merged_ns_seq[col], fit_reg=False, scatter_kws={'s': 2})
    sub_merged = merged_ns_seq[col]
    kendall = sub_merged[["NanoString", "RNAseq_lex"]].corr(method="kendall").iloc[1]["NanoString"]

    plt.yscale("log")
    plt.xscale("log")
    print(concat_lex)
    concat_lex = pd.concat([concat_lex, sub_merged[["NanoString", "RNAseq_lex"]]])
    print(concat_lex)
    plt.title(col)
    count += 1
    dict_lex[col] = r_value ** 2
    #print(col, r_value ** 2)
#plt.subplots_adjust(wspace=0.3, hspace=0.35, top=0.97, right=0.97, bottom=0.05, left=0.04)

slope, intercept, r_value, p_value, std_err = stats.linregress(concat_lex["NanoString"], concat_lex["RNAseq_lex"])
kendall = concat_lex[["NanoString", "RNAseq_lex"]].corr(method="kendall").iloc[1]["NanoString"]
print(kendall)
plt.text(200, 40000, "kendall = {0:.3f}, R-sq.{1:0.3f}".format(kendall, r_value**2),
            horizontalalignment='center',
            verticalalignment='center')

plt.title("Lexogen")
plt.savefig("/Users/nicholas/Desktop/NanoString_vs_RNAseq_lex_combined.png", dpi=500)
plt.show()


count = 1
concat_sal = pd.DataFrame(columns=["NanoString", "RNAseq_sal"])
for col in merged_ns_seq.columns.levels[0]:
    slope, intercept, r_value, p_value, std_err = stats.linregress(merged_ns_seq[col]["NanoString"],
                                                                   merged_ns_seq[col]["RNAseq_sal"])
    sub_merged = merged_ns_seq[col]
    kendall = sub_merged[["NanoString", "RNAseq_sal"]].corr(method="kendall").iloc[1]["NanoString"]
    Rsq = r_value ** 2

    sns.regplot(x="NanoString", y="RNAseq_sal", data=merged_ns_seq[col], fit_reg=False, scatter_kws={'s': 2})
    concat_sal = pd.concat([concat_sal, sub_merged[["NanoString", "RNAseq_sal"]]])
    plt.yscale("log")
    plt.xscale("log")
    plt.title("Salmon")
    plt.title(col)
    count += 1
    dict_sal[col] = r_value ** 2
    #print(col, r_value ** 2)

slope, intercept, r_value, p_value, std_err = stats.linregress(concat_sal["NanoString"], concat_sal["RNAseq_sal"])

kendall = concat_sal[["NanoString", "RNAseq_sal"]].corr(method="kendall").iloc[1]["NanoString"]
print(kendall)
plt.text(200, 40000, "kendall = {0:.3f}, R-sq.{1:0.3f}".format(kendall, r_value**2),
            horizontalalignment='center',
            verticalalignment='center')

plt.title("Salmon")
plt.savefig("/Users/nicholas/Desktop/NanoString_vs_RNAseq_sal_combined.png", dpi=500)
plt.show()