import pandas as pd # 1.5.3

# Loads ped file into a pandas dataframe
indiv_ped = pd.read_table("Data/Example1/indiv.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
indiv_map = pd.read_table("Data/Example1/indiv.map", sep="\t", header=None)


# sort cell content
sort_indiv_ped = indiv_ped.applymap(lambda x: ''.join(sorted(x)) if isinstance(x, str) else x)
sort_indiv_peds = sort_indiv_ped
sort_indiv_peds.iloc[:,6:]= indiv_ped.iloc[:,6:].applymap(lambda x: ''.join(sorted(x)) if isinstance(x, str) else x)

#crate a proper function


def order_snp (df):

    df.iloc[:,6:] = df.iloc[:,6:].applymap(lambda x: ''.join(sorted(x)).strip() if isinstance(x, str) else x)
    return df

sorted = order_snp(indiv_ped)
a= sorted.dtypes
# now apply them only to the snp columns
