from turtle import color
from matplotlib import markers, style
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv("table_for_plots.tsv", delimiter="\t")

plt.figure(figsize=(20,13))
sns.stripplot(y="pan_transcriptome_size", x = "number_genotypes", data = data, jitter=True, marker="o", alpha=0.3, palette="tab10")
plt.savefig("raw_saturation.png")



