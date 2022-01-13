import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats


df_corr = pd.read_csv("Data/Corr_variables.txt", sep='\t', index_col=0)
print(df_corr)

print(df_corr["Rsquared_lex"].median())
print(df_corr["Kendall_lex"].median())


plt.figure(figsize=(10, 7))
sns.regplot(df_corr["Rsquared_lex"], df_corr["DV200"], label="HTSeq,           {:.2f}".format(df_corr["Rsquared_lex"].median()))
sns.regplot(df_corr["Rsquared_sal"], df_corr["DV200"], label="Salmon,          {:.2f}".format(df_corr["Rsquared_sal"].median()))
plt.xlabel("R-squared - NanoString versus RNAseq")
plt.legend(title="Quantification Algorithm, Median R-sq.")
plt.savefig("/Users/nicholas/Desktop/Corrlation_Analysis_Rsq.png", dpi=500)
#plt.show()

plt.figure(figsize=(10, 7))
sns.regplot(df_corr["Kendall_lex"], df_corr["DV200"], label="HTSeq,           {:.2f}".format(df_corr["Kendall_lex"].median()))
sns.regplot(df_corr["Kendall_sal"], df_corr["DV200"], label="Salmon,          {:.2f}".format(df_corr["Kendall_sal"].median()))
plt.xlabel("Kendall's Tau - NanoString versus RNAseq")
plt.legend(title="Quantification Algorithm, Median Kendall")
plt.savefig("/Users/nicholas/Desktop/Corrlation_Analysis_kendall.png", dpi=500)
plt.show()