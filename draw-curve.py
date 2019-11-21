import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import seaborn as sns


def main():
    figname = sys.argv[-1]
    title   = sys.argv[2]
    datas   = pd.read_csv(sys.argv[1])
    
    plot  = sns.lineplot(x=datas["size"], y=datas["perf"],
                         hue=datas['function'],
                         err_style='band', legend='full')
    plt.xlabel("Size (n)")
    plt.ylabel("Performance (Flops/s)")
    plt.title(title)
    plt.savefig(figname)
    """with open(sys.argv[1], "r") as curvefile:
        curveslines = curvefile.readlines()
    for line in curveslines:
        if not line:continue
        if not line[0].isnumeric(): continue
        size, value, duration = map(float, line.split(","))
        size = int(size)"""
        
main()
