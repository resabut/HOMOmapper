#! /usr/bin/env python3
"""
usage: HOMOmapper.py [-h] -i indiv_ped -r ref_ped [-o output_txt] [-s] [-sl short_len]
Description: HOMOmapper is a tool for finding homologous regions between two
individuals. It uses a reference individual to find homologous regions in the
other individual. The reference individual can be a close relative of the
individual of interest, or a population representative. The tool uses a
sliding window approach to find homologous regions. The window size is set by
the minimum ROH length. The maximum number of mismatches within the window is
set by the user. The tool outputs a table with the start and end position of
the ROHs, the chromosome, the number of mismatches and the length of the ROH.
positional arguments:
    indiv_ped             Individual ped file
    ref_ped               Reference ped file
optional arguments:
    -h, --help            show this help message and exit
    -o output_txt, --out output_txt
                            Output file (default: out_roh_table.txt)
    -s, --short           Shorten ped files. Useful for testing
    -sl short_len, --short_len short_len
                            Number of columns (Phenotype + genotype data) to be kept in the shortened ped files. By default, 500 columns. Useful for testing (default: 500)
    -l min_roh_length, --min-roh-length min_roh_length
                            Minimum ROH length. It is used for the length of the sliding window (default: 100)
    -m max_mismatch, --max-mismatch max_mismatch
                            Maximum number of mismatches within the sliding window (default: 5)

"""
import random
import time
import pandas as pd  # 1.5.3
import argparse
from pathlib import Path
import sys
from tqdm import tqdm
import subprocess

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
parser.add_argument('-b', '--bash', action="store_true", help='Use bash to sort ped files. '
                                                              'Faster but requires UNIX', required=False)

args = parser.parse_args()

# define basic variables
start_total_time = time.time()
# config
# minimum length of ROH
min_roh = args.min_roh_length
# max number of mismatches
max_mismatch = args.max_mismatch

# define functions
def find_map_file(ped_file):
    map_file = ped_file.with_suffix(".map")  # check fi it works with ped file
    if Path(map_file).is_file():
        print(f"Map file {map_file.name} found. Continuing...")
    else:
        print(f"Map file {map_file.name} not found.\nPlease, make sure that the .map file is in the same path folder \
        as the .ped file and has the same name. Then, run again.\nExiting...")
        sys.exit()
    return map_file

def sort_bash():
    print("Sorting ped files with bash... (Used -b/--bash option)")
    # Define the Bash command to run, including the user input
    bash_command_indiv = "gawk -F '\\t' '{for(i=7;i<=NF;i++) {split($i,a,\"\"); asort(a); $i=\"\"; for(j=1;j<=length(a);j++) { $i=$i\"\"a[j]; } }}1' " \
                         + str(args.indiv) + " | tr -s ' ' '\t' | sed 's/\t$//' > temp_ind.ped"  # it is a  Path object
    bash_command_ref = "gawk -F '\\t' '{for(i=7;i<=NF;i++) {split($i,a,\"\"); asort(a); $i=\"\"; for(j=1;j<=length(a);j++) { $i=$i\"\"a[j]; } }}1' " \
                       + str(args.ref) + " | tr -s ' ' '\t' | sed 's/\t$//' > temp_ref.ped"  # it is a  Path object

    # Run the Bash command and capture its output
    a = subprocess.check_output(bash_command_indiv, shell=True)
    b = subprocess.check_output(bash_command_ref, shell=True)
    # print(sorted)
    args.indiv = Path("temp_ind.ped")
    args.ref = Path("temp_ref.ped")
    # exit()

def create_all_data_df(indiv_map, indiv_ped, ref_ped):
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
    elif args.short_len != 500:  # the user has specified a length but not the short option
        print("\033[1mIgnoring short length option. Please, use the -s option to shorten the data.\033[0m")
        print("Continuing...")
    # print(all_data_df)
    return all_data_df

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

def find_matches(df):
    print("------------------------------------")
    print("Finding matches...")
    matches = df.copy()
    # turn all loci to 1 ~ they'll become 0 if they match
    matches.iloc[1:, 6:] = 1
    # find matches

    start_time = time.time()
    for pair in tqdm(range(2, len(df.index)), ascii="░▒█", leave=False):  # number of pairwise comparisons
        # print(f"Pair {pair}")
        for snp in range(6, len(df.columns)):  # number of snp
            # print(f"SNP number {snp}")
            if df.iloc[1, snp] == df.iloc[pair, snp]:
                # match is 0
                matches.iloc[pair, snp] = 0
    end_time = time.time()
    print(f"Time to find matches: {end_time - start_time:.4f}")
    return matches


def find_roh(df):
    print("------------------------------------")
    print("Finding ROH...")
    roh_df_temp = pd.DataFrame(
        {"Pair": [], "Chromosome": [], "First_SNP": [], "Last_SNP": [], "Mismatch": [], "Original_row": []})
    indiv_str = "-".join([str(x) for x in df.iloc[1, 0:2].tolist()])
    # print(indiv_str)
    start_time = time.time()
    for pair in tqdm(range(2, len(df.index)), ascii="░▒█", leave=False):  # number of pairwise comparisons
        ref_str = "-".join([str(x) for x in df.iloc[pair, 0:2].tolist()])
        pair_str = ";".join([indiv_str, ref_str])
        # print("Comparing ", pair_str)
        roh_found = 0
        for chr_nr in df.iloc[0, 6:].unique().tolist():
            # print("Chromosome", chr_nr)
            trans_df = df.T  # make a transposed copy to make the chromosome selection easier
            chr_df = trans_df[trans_df["Chromosome"] == chr_nr].T  # transpose again
            # print(df.iloc[df.iloc[0, :] == chr_nr])
            num_loci = len(chr_df.columns)
            # print("Number of loci", num_loci)
            for start_pos in range(0, num_loci - min_roh):
                # print("Start position", start_pos)
                window = chr_df.iloc[pair, start_pos:start_pos + min_roh]
                # print(window)
                # print the window values

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
        # print(f"Found {roh_found} ROHs in pair {pair_str}.")

    end_time = time.time()
    print(f"Time to find ROH: {end_time - start_time:.4f}")
    roh_df_temp[roh_df_temp.columns[1:]] = roh_df_temp[roh_df_temp.columns[1:]].astype(int)
    return roh_df_temp

def find_overlap_roh(df):
    df['overlap'] = df['Last_SNP'] - df['First_SNP'].shift(-1)
    # if overlap is positive, then the roh are overlapping. If it is 0, they are contiguous
    roh_candidates = df[df['overlap'] >= 0]

    print(roh_candidates)

    # join roh and evaluate number of mismatches
    # udpate roh table

def extend_roh(roh_df, match_df, mismatch_threshold):
    print("------------------------------------")
    print("Extending ROH...")
    start_time = time.time()
    # find roh that are not overlapping
    # find roh that are overlapping
    # extend roh
    # loop through every roh, every row
    for index, roh in tqdm(roh_df.iterrows(), ascii="░▒█", leave=False):

        # print("ROH")
        # print(roh)
        # print("Extending ROH", index)
        up_limit = False
        down_limit = False
        while roh['Mismatch'] <= mismatch_threshold:
            # print("Extending ROH", index)

            # get the chr subset
            trans_df = match_df.T  # make a transposed copy to make the chromosome selection easier
            chr_df = trans_df[trans_df["Chromosome"] == roh["Chromosome"]].T  # transpose again
            num_loci = len(chr_df.columns)

            # extend roh

            # new limits
            new_first_snp = roh['First_SNP'] - 1
            new_last_snp = roh['Last_SNP'] + 1
            # check that the snp are within the chromosome
            if new_first_snp < 0:
                # print("First SNP is out of range")
                up_limit = True
                new_first_snp = 0  # reset position
            if new_last_snp > num_loci:
                # print("Last SNP is out of range")
                down_limit = True
                new_last_snp = num_loci
            # calc mismatch vals for new positions
            #check for NA values
            # print("Checking for NA values")
            # print(chr_df.isna(), chr_df.isnull())
            # print(num_loci)

            up_mismatch = chr_df.iloc[roh['Original_row'], new_first_snp]  # upstream
            down_mismatch = chr_df.iloc[roh['Original_row'], new_last_snp - 1]  # downstream
            # print("Got downstream and upstream mismatch values")
            # print()
            # evaluate the mismatch values
            if up_mismatch == 0:  # the snp upstream is a match
                # print("good up")
                roh['First_SNP'] = new_first_snp
                # no need to update mismatch, since up mismatch is 0
                if down_mismatch == 0 or roh["Mismatch"] + 1 <= mismatch_threshold:
                    # the snp downstream is a match too or the threshold is not exceeded (if it is mismatch)
                    # print("good down")
                    roh['Last_SNP'] = new_last_snp
                    roh['Mismatch'] = roh["Mismatch"] + down_mismatch  # update mismatch in case it is 1
                elif up_limit:  # the snp downstream is a mismatch and the threshold is exceeded and the upstream limit is reached
                    # print("bad down")
                    break
            elif down_mismatch == 0:  # the snp downstream is a match, but the upstream is a mismatch
                # print("only good down")
                # extend only the downstream
                roh['Last_SNP'] = new_last_snp
                # no need to update mismatch, since down mismatch is 0
                if roh["Mismatch"] + 1 <= mismatch_threshold:
                    # the snp upstream is a match too or the threshold is not exceeded (if it is mismatch)
                    # print("good up")
                    roh['First_SNP'] = new_first_snp
                    roh['Mismatch'] = roh["Mismatch"] + up_mismatch
                elif down_limit:  # the snp upstream is a mismatch and the threshold is exceeded
                    # and the downstream limit is reached
                    # no more extension possible
                    break

            else:  # both ends are mismatches
                # check if the new limits are within the threshold
                if roh["Mismatch"] + 2 <= mismatch_threshold:
                    # Can extend both ends under the threshold
                    roh['First_SNP'] = new_first_snp
                    roh['Last_SNP'] = new_last_snp
                    roh['Mismatch'] = roh["Mismatch"] + up_mismatch + down_mismatch
                elif roh["Mismatch"] + 1 <= mismatch_threshold:  # extending both ends exceeds the threshold, but one doesn't
                    # extend only one end, choose randomly
                    if random.randint(0, 1) and not up_limit:  # up, and upstream limit not reached
                        roh['First_SNP'] = new_first_snp
                        roh['Mismatch'] = roh["Mismatch"] + up_mismatch
                    elif not down_limit:  # down and downstream limit not reached
                        roh['Last_SNP'] = new_last_snp
                        roh['Mismatch'] = roh["Mismatch"] + down_mismatch
                    else:  # both limits reached, no extension is possible
                        # break the loop
                        break
                else:  # extending both ends exceeds the threshold, no extension is possible
                    # break the loop
                    break
            # update roh table
            roh_df.iloc[index, :] = roh
            # check if both ends of the chromosome have been reached
            if up_limit and down_limit:
                break
            # print(roh)
        # print("Done extending ROH", index)

    # also remove the original row column
    roh_df = roh_df.drop(columns=['Original_row'])
    end_time = time.time()
    print(f"Time to extend ROH: {end_time - start_time:.4f}")
    return roh_df

def eliminate_duplicate_roh(df):
    print("------------------------------------")
    print("Eliminating duplicate ROH...")
    start_time = time.time()
    df = df.drop_duplicates(subset=['Pair', 'Chromosome', 'First_SNP', 'Last_SNP'], keep='first')
    end_time = time.time()
    print(f"Time to eliminate duplicate ROH: {end_time - start_time:.4f}")
    return df

def find_genome_position(df, map_df):
    print("------------------------------------")
    print("Finding genome position...")
    start_time = time.time()
    for index, roh in tqdm(df.iterrows(), ascii="░▒█", leave=False):
        df.loc[index, 'First_Genome_Pos'] = map_df.iloc[roh['First_SNP'], 3].astype(int)
        df.loc[index, 'Last_Genome_Pos'] = map_df.iloc[roh['Last_SNP'], 3].astype(int)
        # print()
    end_time = time.time()
    print(f"Time to find genome position: {end_time - start_time:.4f}")
    return df

def save_to_file(df, filename):
    print("------------------------------------")
    print("Saving to file...")
    start_time = time.time()
    df.to_csv(filename, sep='\t', index=False)
    end_time = time.time()
    print(f"Time to save to file: {end_time - start_time:.4f}")

### MAIN ###
print("Loading data...")

# find the corresponding map file
ind_map_file = find_map_file(args.indiv)
ref_map_file = find_map_file(args.ref)

# Faster way with bash
if args.bash:
    sort_bash()


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
all_data_df=create_all_data_df(indiv_map, indiv_ped, ref_ped)


if not args.bash: # if not using bash, sort the ped files with pandas
    # sort the ped files
    all_data_df.iloc[1:, :] = ordersnp(all_data_df.iloc[1:, :])


# find matches
matches_df = find_matches(all_data_df)
print(matches_df)

# find homologous regions
preroh_df = find_roh(matches_df)
print(preroh_df)

# join overlapping roh
# find overlapping roh
# find_overlap_roh(roh_df)
# deprecated

# extend the existing roh until threshold is reached
extended_roh_df = extend_roh(preroh_df, matches_df, max_mismatch)
print(extended_roh_df)

# eliminate duplicate roh
clean_roh_df = eliminate_duplicate_roh(extended_roh_df)
print(clean_roh_df)

# find genome position
final_df = find_genome_position(clean_roh_df, indiv_map)
print(final_df)

# save to file
save_to_file(final_df, args.out)

print("------------------------------------")
end_total_time = time.time()

print(f"Total time: {end_total_time - start_total_time:.4f}")
print("Done.")
