import pandas as pd # 1.5.3

# Loads ped file into a pandas dataframe
indiv_ped = pd.read_table("Data/Example1/indiv.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
indiv_map = pd.read_table("Data/Example1/indiv.map", sep="\t", header=None)


# this function sorts the snp column content

def order_snp (df):
    #
    df[df.columns[6:]] = df.iloc[:,6:].applymap(lambda x: int(''.join(sorted(x)).strip()) if isinstance(x, str) else x)
    return df

sorted = order_snp(indiv_ped)
# now apply them only to the snp columns
