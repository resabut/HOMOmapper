import pandas as pd # 1.5.3

# Loads ped file into a pandas dataframe
indiv_ped = pd.read_table("Data/Example1/indiv.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
indiv_map = pd.read_table("Data/Example1/indiv.map", sep="\t", header=None)

# Loads ped file into a pandas dataframe
ref_ped = pd.read_table("Data/Example1/ref.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
red_map = pd.read_table("Data/Example1/ref.map", sep="\t", header=None)


# this function sorts the snp column content

def order_snp (df):
    # creates genotype cells
    # ignores the first 6 columns
    # the rest of the columns are sorted to  eliminate genotype 13 and 31 ambiguity
    df[df.columns[6:]] = df.iloc[:,6:].applymap(lambda x: int(''.join(sorted(x)).strip()) if isinstance(x, str) else x)
    # returns the dataframe with the genotype data where each column is a loci
    return df


# sort the ped files
sorted_indiv_ped = order_snp(indiv_ped)
sorted_ref_ped = order_snp(ref_ped)



# now find loci with same genotype in both individuals and reference individuals

for i in range(0, len(sorted_indiv_ped.columns)):
    # this is the loci name
    loci = sorted_indiv_ped.iloc[:,1][i-6]
    # this is the genotype of the individual
    indiv_genotype = sorted_indiv_ped.iloc[:,i][0]
    # this is the genotype of the reference individual
    ref_genotype = sorted_ref_ped.iloc[:,i][0]
    # if the genotypes are the same
    if indiv_genotype == ref_genotype:
        # print the loci name
        print(loci)