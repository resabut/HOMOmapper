"""
info
"""

import pandas as pd  # 1.5.3
# remove false positve warning related to chained assignment
pd.options.mode.chained_assignment = None  # default='warn'
import time

# config
# minimum length of ROH - it used for the window scanning
min_roh = 15
# max number of mismatches
max_mismatch = 2
# for testing, shorten ped files to 50 cols
short = True
short_len = 50

print("Loading data...")

# Loads ped file into a pandas dataframe
indiv_ped = pd.read_table("Data/Example1/indiv.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
indiv_map = pd.read_table("Data/Example1/indiv.map", sep="\t", header=None)

# Loads ped file into a pandas dataframe
ref_ped = pd.read_table("Data/Example1/ref.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
ref_map = pd.read_table("Data/Example1/ref.map", sep="\t", header=None)

# check if map files are the same
print("Adding chromosome information...")

# include chromosome information in the ped file
empty_cols = dict(zip(range(6), [None] * 6))
chr_row = pd.concat([pd.DataFrame(empty_cols, index=[0]), pd.DataFrame(indiv_map[0]).T], axis=1, ignore_index=True)
# rename rows
chr_row.rename(index={0: "Chromosome"}, inplace=True)
indiv_ped.rename(index={0: "Ind"}, inplace=True)
# concate the chromosome row and the two ped files
all_data_df = pd.concat([chr_row, indiv_ped, ref_ped], axis=0)
# make temporary shorter ones, if preferred
if short:
    all_data_df = all_data_df.iloc[:, 0:short_len]

# this function sorts the snp column content
def order_snp(df):
    print("Sorting data...")
    start_time = time.time()
    # creates genotype cells
    # ignores the first 6 columns
    # the rest of the columns are sorted to  eliminate genotype 13 and 31 ambiguity
    df[df.columns[6:]] = df[df.columns[6:]].applymap(lambda x: int(''.join(sorted(x)).strip()) if isinstance(x, str) else x)
    # returns the dataframe with the genotype data where each column is a loci
    end_time = time.time()
    print(f"Time to sort data: {end_time - start_time:.4f}")
    return df



# sort the ped files

all_data_df.iloc[1:, :] = order_snp(all_data_df.iloc[1:, :])




# find matches
def find_matches(df):
    matches = df.copy()
    # turn all loci to 1 ~ they'll become 0 if they match
    matches.iloc[1:, 6:] = 1
    # find matches
    print("Finding matches...")
    start_time = time.time()
    for pair in range(2, len(df.index)): # number of pairwise comparisons
        # print(f"Pair {pair}")
        for snp in range(6, len(df.columns)): # number of snp
            # print(f"SNP number {snp}")
            if df.iloc[1, snp] == df.iloc[pair, snp]:
                # match is 0
                matches.iloc[pair, snp] = 0
    end_time = time.time()
    print(f"Time to find matches: {end_time - start_time:.4f}")
    return matches


matches_df = find_matches(all_data_df)

# find homologous regions
def find_roh(df):
