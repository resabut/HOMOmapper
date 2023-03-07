import pandas as pd # 1.5.3
import time
# config
# minimum length of ROH - it used for the window scanning
min_roh = 15
# max number of mismatches
max_mismatch = 2
# for testing, shorten ped files to 50 cols
short = True
short_len = 10000


print("Loading data...")

# Loads ped file into a pandas dataframe
indiv_ped = pd.read_table("Data/Example1/indiv.ped", sep="\t", header=None)


# Loads map file into a pandas dataframe
indiv_map = pd.read_table("Data/Example1/indiv.map", sep="\t", header=None)

# Loads ped file into a pandas dataframe
ref_ped = pd.read_table("Data/Example1/ref.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
ref_map = pd.read_table("Data/Example1/ref.map", sep="\t", header=None)

# include chromosome information in the ped file

# remove first 6 columns
indiv_ped = indiv_ped.iloc[:,6:]
new_col_names = dict(zip(indiv_ped.columns, range(0,len(indiv_ped.columns))))
indiv_ped = indiv_ped.rename(columns=new_col_names)
ref_ped = ref_ped.iloc[:,6:].reset_index(drop=True)
ref_ped = ref_ped.rename(columns=new_col_names)
# include chromosome information in the ped file
indiv_ped= pd.concat([indiv_ped, pd.DataFrame(indiv_map[0]).T]) # transpose the chromosome column


# this function sorts the snp column content
def order_snp (df):
    # creates genotype cells
    # ignores the first 6 columns
    # the rest of the columns are sorted to  eliminate genotype 13 and 31 ambiguity
    df = df.applymap(lambda x: int(''.join(sorted(x)).strip()) if isinstance(x, str) else x)
    # returns the dataframe with the genotype data where each column is a loci
    return df
print("Sorting data...")
#make temporary shorter ones
if short:
    indiv_ped = indiv_ped.iloc[:,0:short_len]
    ref_ped = ref_ped.iloc[:,0:short_len]
# sort the ped files
start_time = time.time()
sorted_indiv_ped = order_snp(indiv_ped)
sorted_ref_ped = order_snp(ref_ped)
end_time = time.time()

execution_time = end_time - start_time
print(f"Sorting: Execution time: {execution_time:.4f} seconds")



# now find loci with same genotype in both individuals and reference individuals
# craete empty dataframe with the matches
matches = sorted_ref_ped.copy()
# turn all loci to 1 ~ they'll become 0 if they match
matches.fillna(value=1, inplace=True)
# loops through every column
print("Finding matches...")
for row in range(0, len(sorted_ref_ped.index)):
    pers = ref_ped.iloc[:, 6][row] # recovers metadata from the original ped file
    for loci in range(0, len(sorted_indiv_ped.columns)):
        # this is the genotype of the individual
        indiv_genotype = sorted_indiv_ped.iloc[0,loci] # always the first row, 2nd is the chr info
        # this is the genotype of the reference individual
        ref_genotype = sorted_ref_ped.iloc[:,loci][row]
        # print(pers, indiv_genotype, ref_genotype)
        # if the genotypes are the same
        if indiv_genotype == ref_genotype:
            # match is 0 ~ 0 distance
            matches.iloc[row, loci] = 0
# print(matches)

print("Finding ROH...")
# scan for ROH with moving window

for row in range(0, len(sorted_ref_ped.index)): # for every individual
    print("Individual", row)
    for chr_nr in sorted_indiv_ped.iloc[1,:].unique().tolist(): # for every chromosome
        print("Chromosome", chr_nr)
        # get the data for the chromosome
        chr_data_ref = matches[sorted_indiv_ped.iloc[1,:] == chr_nr]
        num_loci = len(chr_data_ref)  # number of loci
        print("In this chromosome:", num_loci, "loci")
        for start_pos in range(0, num_loci):
            # if the window is not longer than the minimum ROH
            if start_pos + min_roh <= num_loci:
                # get the window sum
                window = chr_data_ref.iloc[row, start_pos:(start_pos + min_roh)]
                # if the window has no mismatches
                # print(row, window, window.values.tolist(), "window sum", window.sum())
                if window.sum() <= max_mismatch: # mismatch limit
                    # print the start and end position of the ROH
                    print("ROH found")
                    print(row, start_pos, start_pos + min_roh, "Mismatch:", window.sum())