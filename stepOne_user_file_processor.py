def user_file_processor(user_vcf_file_path, first_letter_unwanted_sample):

  """
  This function performs the following processing tasks onto the user's exome data file:
  (1) Removes sample IDs that start with a particular letter (e.g., "T")
  (2) Removes SNPs that are INDELS
  (3) Converts genotypes to 0, 1, 2 and 3 values
  (4) Removes SNPs for chromosome Y

  """
  print("=====================================================================================================================")

  # Step 1: Importing useful libraries
  import numpy as np
  import pandas as pd
  import matplotlib.pyplot as plt
  import seaborn as sns
  import io
  import os

  # Step 2: Creating a function that stores data from .vcf file to Pandas DataFrame
  def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

  # Step 3: Running the read_vcf() function and storing the vcf data in a Pandas DataFrame
  df = read_vcf(user_vcf_file_path)
  print("Before removing unnecessary samples: ", df.shape[0], "SNPs and", df.shape[1], "samples.")

  # Step 4: Removing columns that hold Temiar data as well as 'ID', 'QUAL', 'FILTER', 'INFO' and 'FORMAT' columns
  # Initialising a new DataFrame
  df2 = df

  for i in range(df.shape[1]):

    # Extracting each column title
    col_title = df.columns[i]

    # Extracting the first character of each column title
    first_char = col_title[0]

    # Removing the columns of not our interest
    if first_char == first_letter_unwanted_sample:
      df2 = df2.drop(columns = col_title, axis=1)
    elif col_title == "ID":
      df2 = df2.drop(columns = col_title, axis=1)
    elif col_title == "QUAL":
      df2 = df2.drop(columns = col_title, axis=1)
    elif col_title == "FILTER":
      df2 = df2.drop(columns = col_title, axis=1)
    elif col_title == "INFO":
      df2 = df2.drop(columns = col_title, axis=1)
    elif col_title == "FORMAT":
      df2 = df2.drop(columns = col_title, axis=1)

  # Step 5: Saving the new DataFrame in a csv file
  print("After removing unnecessary samples: ", df2.shape[0], "SNPs and", df2.shape[1], "samples.")

  # Step 6: Removing SNPs of indel mutation
  # Initialising a new DataFrame
  #df3 = df2
  index_INDEL_list = []

  for i in range(df2.shape[0]):

    # Extracting the length of REF data for each row
    len_ref = len(df2.iloc[i,2])

    # Extracting the length of ALT data for each row
    len_alt = len(df2.iloc[i,3])

    if len_ref > 1 or len_alt > 1:
      index_INDEL_list.append(i)

  df2 = df2.drop(index=index_INDEL_list, axis=0).reset_index(drop=True)

  print("After removing INDELS: ", df2.shape[0], "SNPs.")

  # Step 7: Encoding data entries to just genotypes in the form of 0, 1, 2 and 3

  # Converting the DataFrame into a numpy array
  df2_np = np.array(df2.iloc[:, 4:])

  for i in range(df2_np.shape[0]):
    for j in range(df2_np.shape[1]):
      if df2_np[i,j][0:3] == "0/0" or df2_np[i,j][0:3] == "0|0":
        df2_np[i,j] = 2

      elif df2_np[i,j][0:3] == "0/1" or df2_np[i,j][0:3] == "0|1" or df2_np[i,j][0:3] == "1/0" or df2_np[i,j][0:3] == "1|0":
        df2_np[i,j] = 1

      elif df2_np[i,j][0:3] == "1/1" or df2_np[i,j][0:3] == "1|1":
        df2_np[i,j] = 0

      else:
        df2_np[i,j] = 3

  df2.iloc[:, 4:] = pd.DataFrame(df2_np)

  print("- Genotypes have been encoded as numerical values successfully!")

  # Step 8: Removing rows that do not contain any mutations (i.e., only 0/0 genotypes)
  # Initialising an empty list to store row indices of data entries that do not contain any mutations
  index_no_mutation_list = []

  for i in range(df2_np.shape[0]):

    genotype_count_df = pd.DataFrame(df2.iloc[i,4:].value_counts())

    # if there is only one type of genotype, and the type of genotype is of value 2, and there are (df2.shape[1] - 4) counts of it
    if genotype_count_df.shape[0] == 1 and genotype_count_df.index[0] == 2 and genotype_count_df.iloc[0,0] == df2.shape[1] - 4:
      index_no_mutation_list.append(i)

  df2 = df2.drop(index=index_no_mutation_list, axis=0).reset_index(drop=True)

  print("After removing rows without mutations: ", df2.shape[0], "SNPs.")

  # Step 9: Removing chrY rows
  # Getting indices for rows that pertain to chrY
  index_Y_list = [index for index, item in enumerate(list(df2["CHROM"])) if item == "chrY"]

  # Dropping rows of chrY
  df2 = df2.drop(index_Y_list)

  # Saving dataframe as a csv file
  df2.to_csv("User_Exome_Analysis.csv", index=False)
  print("After removing chromosome Y observations: ", df2.shape[0], "SNPs.")
  print("- User file has been processed successfully!")

  print("=====================================================================================================================")


  #return df2