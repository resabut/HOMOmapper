#! /usr/bin/env python3
"""
usage: HOMOmapper.py [-h] -i indiv_ped -r ref_ped [-o output_txt] [-s] [-sl short_len]
Description: HOMOmapper is a tool for finding homologous regions between two
individuals. It uses a reference individual to find homologous regions in the
other individual. The reference individual can be a close relative of the
individual of interest, or a population representative. The tool uses a
sliding window approach to find homologous regions. The window size is set by

Description: HOMOmapper is a tool for finding homologous regions between two
individuals. It uses a reference individual to find homologous regions in the
other individual. The reference individual can be a close relative of the
individual of interest, or a population representative. The tool uses a
"""
import random
import time
import pandas as pd  # 1.5.3
import argparse
from pathlib import Path
import sys

# remove false positve warning related to chained assignment
pd.options.mode.chained_assignment = None  # default='warn'

# Get argparse arguments
parser = argparse.ArgumentParser(description='HOMOmapper')
parser.add_argument('-i', '--indiv', help='Individual ped file', metavar="indiv_ped",
                    type=Path, required=True)
parser.add_argument('-r', '--ref', help='Reference ped file', metavar="ref_ped",
                    type=Path, required=True)
parser.add_argument('-o', '--out', help='Output file', metavar="output_txt", type=Path,
                    default="out_roh_table.txt", required=False)
parser.add_argument('-s', '--short', action="store_true", help='Shorten ped files. Useful for testing', required=False)
parser.add_argument('-sl', '--short_len',
                    help='Number of columns (Phenotype + genotype data) to be kept in the shortened ped files. \
                    By default, 500 columns. Useful for testing',
                    type=int, default=500, required=False)
parser.add_argument('-l', '--min-roh-length',
                    help='Minimum ROH length. It is used for the length of the sliding window',
                    type=int, required=True)
parser.add_argument('-m', '--max-mismatch', help='Maximum number of mismatches within an ROH.',
                    type=int, required=True)
args = parser.parse_args()

# config
# minimum length of ROH
min_roh = args.min_roh_length
# max number of mismatches
max_mismatch = args.max_mismatch

print("Loading data...")


# find the corresponding map file
def find_map_file(ped_file):
    map_file = ped_file.with_suffix(".map")  # check fi it works with ped file
    if Path(map_file).is_file():
        print(f"Map file {map_file.name} found. Continuing...")
    else:
        print(f"Map file {map_file.name} not found.\nPlease, make sure that the .map file is in the same path folder \
        as the .ped file and has the same name. Then, run again.\nExiting...")
        sys.exit()
    return map_file


ind_map_file = find_map_file(args.indiv)
ref_map_file = find_map_file(args.ref)

# Loads ped file into a pandas dataframe
indiv_ped = pd.read_table(args.indiv, sep="\t", header=None)

# Loads map file into a pandas dataframe
indiv_map = pd.read_table(ind_map_file, sep="\t", header=None)

# Loads ped file into a pandas dataframe
ref_ped = pd.read_table(args.ref, sep="\t", header=None)

# Loads map file into a pandas dataframe
ref_map = pd.read_table(ref_map_file, sep="\t", header=None)

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
if args.short:
    print("Shortening data...")
    all_data_df = all_data_df.iloc[:, 0:args.short_len]


# this function sorts the snp column content
def ordersnp(df):
    print("------------------------------------")
    print("Sorting data...")
    start_time = time.time()
    # creates genotype cells
    # ignores the first 6 columns
    # the rest of the columns are sorted to  eliminate genotype 13 and 31 ambiguity
    df[df.columns[6:]] = df[df.columns[6:]].applymap(
        lambda x: int(''.join(sorted(x)).strip()) if isinstance(x, str) else x)
    # returns the dataframe with the genotype data where each column is a loci
    end_time = time.time()
    print(f"Time to sort data: {end_time - start_time:.4f}")
    return df


# sort the ped files

all_data_df.iloc[1:, :] = ordersnp(all_data_df.iloc[1:, :])


# find matches
def find_matches(df):
    print("------------------------------------")
    print("Finding matches...")
    matches = df.copy()
    # turn all loci to 1 ~ they'll become 0 if they match
    matches.iloc[1:, 6:] = 1
    # find matches

    start_time = time.time()
    for pair in range(2, len(df.index)):  # number of pairwise comparisons
        # print(f"Pair {pair}")
        for snp in range(6, len(df.columns)):  # number of snp
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
    print("------------------------------------")
    print("Finding ROH...")
    roh_df_temp = pd.DataFrame(
        {"Pair": [], "Chromosome": [], "First_SNP": [], "Last_SNP": [], "Mismatch": [], "Original_row": []})
    indiv_str = "-".join([str(x) for x in df.iloc[1, 0:2].tolist()])
    # print(indiv_str)
    start_time = time.time()
    for pair in range(2, len(df.index)):  # number of pairwise comparisons
        ref_str = "-".join([str(x) for x in df.iloc[pair, 0:2].tolist()])
        pair_str = ";".join([indiv_str, ref_str])
        print("Comparing ", pair_str)
        roh_found = 0
        for chr_nr in df.iloc[0, 6:].unique().tolist():
            # print("Chromosome", chr_nr)
            trans_df = df.T  # make a transposed copy to make the chromosome selection easier
            chr_df = trans_df[trans_df["Chromosome"] == chr_nr].T  # transpose again
            # print(df.iloc[df.iloc[0, :] == chr_nr])
            num_loci = len(chr_df.columns[6:])
            # print("Number of loci", num_loci)
            for start_pos in range(0, num_loci - min_roh):
                # print("Start position", start_pos)
                window = chr_df.iloc[pair, start_pos:start_pos + min_roh]
                # print(window)
                if window.sum() <= max_mismatch:  # mismatch limit
                    # print("Found ROH")
                    roh_found += 1
                    # print(f"Start position in {chr_nr} in position {start_pos}")
                    new_row = pd.DataFrame(
                        {"Pair": pair_str, "Chromosome": chr_nr, "First_SNP": start_pos,  # it starts at 0
                         "Last_SNP": start_pos + min_roh, "Mismatch": window.sum(),
                         "Original_row": pair},
                        index=[0])
                    roh_df_temp = pd.concat([roh_df_temp, new_row],
                                            ignore_index=True)
        print(f"Found {roh_found} ROHs in pair {pair_str}.")

    end_time = time.time()
    print(f"Time to find ROH: {end_time - start_time:.4f}")
    roh_df_temp.iloc[:, 1:] = roh_df_temp.iloc[:, 1:].astype(int)
    return roh_df_temp


preroh_df = find_roh(matches_df)
print("------------------------------------")
print(preroh_df)


# join overlapping roh
# find overlapping roh
def find_overlap_roh(df):
    df['overlap'] = df['Last_SNP'] - df['First_SNP'].shift(-1)
    # if overlap is positive, then the roh are overlapping. If it is 0, they are contiguous
    roh_candidates = df[df['overlap'] >= 0]

    print(roh_candidates)

    # join roh and evaluate number of mismatches
    # udpate roh table


# find_overlap_roh(roh_df)
# print(len(matches_df.T[matches_df.T['Chromosome'] <=8].T))
# extend the existing roh until threshold
def extend_roh(roh_df, match_df, mismatch_threshold):
    # find roh that are not overlapping
    # find roh that are overlapping
    # extend roh
    # loop through every roh, every row
    for index, roh in roh_df.iterrows():
        # print("ROH")
        # print(roh)

        while roh['Mismatch'] < mismatch_threshold:
            # print("Extending ROH", index)
            extend = False  # reset the extend flag

            # get the chr subset
            trans_df = match_df.T  # make a transposed copy to make the chromosome selection easier
            chr_df = trans_df[trans_df["Chromosome"] == roh["Chromosome"]].T  # transpose again
            num_loci = len(chr_df.columns[6:])

            # extend roh
            # new limits
            new_first_snp = roh['First_SNP'] - 1
            new_last_snp = roh['Last_SNP'] + 1
            # check that the snp are within the chromosome
            if new_first_snp < 0:
                new_first_snp = 0
            elif new_last_snp > num_loci:
                new_last_snp = num_loci
            # calc mismatch vals for new positions
            up_mismatch = chr_df.iloc[roh['Original_row'], new_first_snp]  # upstream
            down_mismatch = chr_df.iloc[roh['Original_row'], new_last_snp]  # downstream
            # evaluate the mismatch values
            if up_mismatch == 0:  # the snp upstream is a match
                # print("good up")
                # it doesn't matter if the down is good or not,
                # the while forces it to be at least one mismatch under the threshold,
                # so the roh is extended
                extend = True
            elif down_mismatch == 0:  # the snp downstream is a match
                # print("good down")
                # idem
                extend = True
            else:  # both ends are mismatches
                # check if the new limits are within the threshold
                if roh["Mismatch"] + 2 <= mismatch_threshold:
                    # Can extend both ends under the threshold
                    extend = True
                else:  # both ends are mismatches, so extending both ends would exceed the threshold
                    # extend only one end, choose randomly
                    if random.randint(0, 1):  # up
                        roh['First_SNP'] = new_first_snp
                        roh['Mismatch'] = roh["Mismatch"] + up_mismatch
                    else:  # down
                        roh['Last_SNP'] = new_last_snp
                        roh['Mismatch'] = roh["Mismatch"] + down_mismatch

                    # update roh table
                    roh_df.iloc[index, :] = roh
            if extend:
                roh['First_SNP'] = new_first_snp
                roh['Last_SNP'] = new_last_snp
                roh['Mismatch'] = roh["Mismatch"] + up_mismatch + down_mismatch
                # update roh table
                roh_df.iloc[index, :] = roh
                # print(roh)

    return roh_df


extended_roh_df = extend_roh(preroh_df, matches_df, max_mismatch)
print(extended_roh_df)

print("Done.")
