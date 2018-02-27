import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class Analyzer:
    def __init__(self, column):
        df = pd.read_csv('output.csv', header=None, sep=',')
        df.columns = column
        print(df.head())

        sns.pairplot(df[column])
        plt.tight_layout()
        plt.savefig('img01.png', dpi=300)
