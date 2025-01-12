import numpy as np
import pandas as pd

# Defining a function that filters the specified chromosome from the user's exome data
def filter_via_chr_func(chr, user_file_path):

    # Reading the csv file from user_file_processor() function
    df_user = pd.read_csv(user_file_path)

    # Step 2: Determining the start and end indices for each chromosome considered in user's data
    index_start_list = [0]
    index_end_list = []

    # Step 3: Determining the chromosomes that are considered in the user's exome data
    chr_list = df_user["CHROM"].unique()

    for i in range(len(list(chr_list))):
        chr_count = df_user["CHROM"].value_counts()[chr_list].iloc[i]
        index_end = index_start_list[i] + chr_count - 1
        index_start = index_end + 1
        index_end_list.append(index_end)

        if i != len(list(chr_list)) - 1:
            index_start_list.append(index_start)

    # Step 4: Determining the user's data indices to be considered based on chr number
    if type(chr) == str:
        JH_index_start = index_start_list[22]
        JH_index_end = index_end_list[22]
    else:
        JH_index_start = index_start_list[chr - 1]
        JH_index_end = index_end_list[chr - 1]

    # Step 5: Filtering the DataFrame based on chr number
    df2 = df_user.iloc[JH_index_start: JH_index_end + 1, :]

    # Returning the filtered DataFrame
    return df2