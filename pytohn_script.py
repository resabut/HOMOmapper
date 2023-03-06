import pandas as pd # 1.5.3

# Loads ped file into a pandas dataframe
indiv_ped = pd.read_table("Data/Example1/indiv.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
indiv_map = pd.read_table("Data/Example1/indiv.map", sep="\t", header=None)


# sort cell content
sort_indiv_ped = indiv_ped.applymap(lambda x: ''.join(sorted(x)) if isinstance(x, str) else x)

#crate a proper function
# now apply them only to the snp columns
