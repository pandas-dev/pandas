import pandas as pd
def plotManhattan(df,chrom_column, pos_column, pval_column, colors=['black','gray']):
    def findOutliers(x):
        m=x.mean()
        std=x.std()
        return  x[x>(m+3*std)]
    from itertools import  cycle
    import pylab as plt
    def plotOne(a,c,name,chroms):
        plt.scatter(a.index, a, s=2, c=c, alpha=0.8, edgecolors='none')
        outliers=findOutliers(a)
        if len(outliers):
            plt.scatter(outliers.index, outliers, s=2, c='r', alpha=0.8, edgecolors='none')
        plt.axis('tight');plt.xlim(0, a.index[-1]);
        plt.xticks([x for x in chroms.mid], [str(x) for x in chroms.index], rotation=-90, size=20);plt.ylim(ymin=0);plt.ylabel(name)
    df=df.sort_values([df.columns[chrom_column],df.columns[pos_column]])
    chroms=pd.DataFrame(df.groupby(df.columns[chrom_column])[df.columns[pos_column]].max()) +1000;chroms.columns=['len']
    print df
    chroms['offset']=chroms.len.cumsum();chroms['offset'].iloc[1:]=chroms['offset'].iloc[:-1].values;chroms['offset'].iloc[0]=0

    colors=cycle(colors)
    chroms['color']=[colors.next() for i in range(chroms.shape[0])]
    chroms['mid']=[x+y/2  for x,y in zip(chroms.offset,chroms.len)]
    df['gpos']=df.iloc[:,pos_column]+ chroms.offset.loc[df.iloc[:,chrom_column]].values
    df['color']=chroms.color.loc[df.iloc[:,chrom_column]].values
    df.set_index('gpos',inplace=True);df.sort_index(inplace=True)
    plt.figure(figsize=(20,6));
    plotOne(df.iloc[:,pval_column], df.color, 'pval',chroms)
    plt.xlabel('Chromosome')

#Download sample data file from  'https://www.dropbox.com/s/n6lqtx4kny8nd0b/AFR.win500K.df?dl=0'
df=pd.read_pickle('/home/arya/Dropbox/genome_scan_results/results/pandas/AFR.win500K.df').reset_index()
plotManhattan(df,chrom_column=0,pos_column=1,pval_column=3)
