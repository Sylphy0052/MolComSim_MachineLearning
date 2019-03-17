import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class Analyzer:
    def __init__(self, column):
        df = pd.read_csv('output.csv', header=None, sep=',')
        df.columns = column
        # print(df.head())

        sns.pairplot(df[column])
        plt.tight_layout()
        # plt.savefig('img/img01.png', dpi=300)
        plt.close('all')

        cm = np.corrcoef(df[column].values.T)
        hm = sns.heatmap(cm,
                 cbar=True,
                 annot=True,
                 square=True,
                 fmt='.2f',
                 annot_kws={'size': 15},
                 yticklabels=column,
                 xticklabels=column)

        plt.tight_layout()
        plt.savefig('img/img02.png', dpi=300)
        plt.close('all')
