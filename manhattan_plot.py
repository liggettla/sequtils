bucket=''

sumstats_list = [
    'mysumstats'
]

for sumstat in sumstats_list:

    sumstats = pd.read_csv(
        'gs://{}/sumstats/{}'.format(bucket, sumstat),
        sep='\t',
        compression='gzip',
    )

    # remove 'chr'
    sumstats.CHR = sumstats.CHR.str.slice(start=3,stop=5,step=1)

    # calculate -logp
    sumstats['-logp'] = -np.log10(sumstats['p.value'])

    # relabel X to treat as int
    sumstats.loc[sumstats.CHR=='X', 'CHR'] = '23'
    sumstats.CHR = sumstats.CHR.astype(int)

    # order the data
    # this doesn't sort in proper order when CHR is str
    running_pos = 0
    cumulative_pos = []

    for chrom, group_df in sumstats.groupby('CHR'):
        cumulative_pos.append(group_df['POS'] + running_pos)
        running_pos += group_df['POS'].max()

    sumstats['cumulative_pos'] = pd.concat(cumulative_pos)

    # reindex
    sumstats['SNP number'] = sumstats.index

    # group and recolor hits
    evens = list(range(2,23,2))

    sumstats['color_group'] = np.where(sumstats['CHR'].isin(evens), 'B', 'A')

    # label significant hits with different color
    # -log10(5x10^-8)
    sumstats['pos_abbrev'] = sumstats.POS.astype(str).str.slice(start=0,stop=2,step=1)
    sig = sumstats[sumstats['-logp'] > 7.301029995663981].drop_duplicates(subset = ['pos_abbrev', 'CHR']).sort_values('POS')

    # recolor significant hits
    for index, row in sig.iterrows():
        sumstats.loc[(sumstats['CHR'] == row['CHR']) &
                     (sumstats['POS'] >= row['POS'] - 300000) &
                     (sumstats['POS'] <= row['POS'] + 300000), 
                     'color_group'] = 'C'

    # downsample 1% -log10(p) < 4
    temp = sumstats[(sumstats.color_group == 'C') | (sumstats['-logp'] > 4)]
    temp = temp.append(sumstats[
        ~(sumstats.color_group == 'C') |
        ~(sumstats['-logp'] > 4)
    ].sample(
        frac=0.01,
        replace=False,
        random_state=1))
    
    g = sns.relplot(
        data = temp.sort_values('color_group'),
        x = 'cumulative_pos',
        y = '-logp',
        aspect = 4, # width and length
        hue = 'color_group',
        # palette = ["#276FBF", "#183059", '#ff0df7'],
        # https://coolors.co/palettes/popular
        # https://i.stack.imgur.com/D5TdK.jpg
        palette = ["#81A1C1", "#5E81AC", '#99d98c'],
        # palette = ["#81A1C1", "#5E81AC", '#94ffda'],
        linewidth=0,
        s=30, # dot size
        legend=None
    )

    # set xtick positions to be at center of chromosomes
    g.ax.set_xticks(sumstats.groupby('CHR')['cumulative_pos'].median())
    
    # convert CHR numbers back to str and relabel 23 as X
    labels = sumstats['CHR'].unique().astype(str)
    labels = np.where(labels == '23', 'X', labels)   
    g.ax.set_xticklabels(labels)

    plt.tight_layout()
    plt.ylim(0,)
    sns.despine(offset=10, trim=True, bottom=False)
    plt.xticks(rotation=90)
    plt.xlabel('')
    plt.ylabel('$-log_{10}$(P)')
    plt.title('{}'.format(sumstat))
    
    # -log10(5x10^-8) significance line
    plt.hlines(y=7.301029995663981, xmin=0, xmax=sumstats.cumulative_pos.max(), colors='#a8dadc', linestyles='dashed', label='',)
    pretty('talk')
    
    # plt.savefig('{}.pdf'.format(sumstat), format="pdf", bbox_inches="tight")

    plt.show()
