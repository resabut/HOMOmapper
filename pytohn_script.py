import pandas as pd # 1.5.3
import time
# config
# minimum length of ROH - it used for the window scanning
min_roh = 100
# max number of mismatches
max_mismatch = 1
# for testing, shorten ped files to 50 cols
short_50 = False


print("Loading data...")

# Loads ped file into a pandas dataframe
indiv_ped = pd.read_table("Data/Example1/indiv.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
indiv_map = pd.read_table("Data/Example1/indiv.map", sep="\t", header=None)

# Loads ped file into a pandas dataframe
ref_ped = pd.read_table("Data/Example1/ref.ped", sep="\t", header=None)

# Loads map file into a pandas dataframe
ref_map = pd.read_table("Data/Example1/ref.map", sep="\t", header=None)


# this function sorts the snp column content
def order_snp (df):
    # creates genotype cells
    # ignores the first 6 columns
    # the rest of the columns are sorted to  eliminate genotype 13 and 31 ambiguity
    df[df.columns[6:]] = df.iloc[:,6:].applymap(lambda x: int(''.join(sorted(x)).strip()) if isinstance(x, str) else x)
    # returns the dataframe with the genotype data where each column is a loci
    return df
print("Sorting data...")
#make temporary shorter ones
if short_50:
    indiv_ped = indiv_ped.iloc[:,0:50]
    ref_ped = ref_ped.iloc[:,0:50]
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
matches.iloc[:,6:] = 1
# loops through every column
print("Finding matches...")
for row in range(0, len(sorted_ref_ped.index)):
    pers = sorted_ref_ped.iloc[:, 6][row]
    for loci in range(6, len(sorted_indiv_ped.columns)):

        # this is the genotype of the individual
        indiv_genotype = sorted_indiv_ped.iloc[:,loci][0]
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
num_loci = len(matches.columns) - 6 # number of loci
print(num_loci)
for row in range(0, len(sorted_ref_ped.index)): # for every individual
    print("Individual", row)
    for start_pos in range(6, num_loci):
        # if the window is not longer than the minimum ROH
        if start_pos + min_roh <= num_loci:
            # get the window sum
            window = matches.iloc[row, start_pos:(start_pos + min_roh)]
            # if the window has no mismatches
            # print(row, window, window.values.tolist(), "window sum", window.sum())
            if window.sum() <= max_mismatch: # mismatch limit
                # print the start and end position of the ROH
                print("ROH found")
                print(row, start_pos, start_pos + min_roh)