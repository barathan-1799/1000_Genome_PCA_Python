def Thousand_Genome_SNP_extractor(chr, user_file_path, batch_size):

  print("==============================================================================================================")
  print("Extracting SNPs from 1000 Genome data for chromosome", chr, "...")

  # Step 1: Importing useful libraries
  import numpy as np
  import pandas as pd
  import io
  import os
  import math
  import gzip
  import requests
  import gc
  from io import BytesIO

  # Step 3: Download the file
  url = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr" + str(chr) + ".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

  # Step 4: Defining a function that stores data from vcf file into a string
  response = requests.get(url, stream=True)
  response.raise_for_status()  # Ensure the request was successful

  # Step 5: Parse the VCF data
  line_list = []
  count = 0
  df = 0

  with gzip.open(BytesIO(response.content), 'rt') as f:
    for l in f:
      if not l.startswith('##'):

        count += 1

        if count == 1:

          # Step 6: Specifying the column titles of the DataFrame
          df = pd.read_csv(io.StringIO(''.join(l)),
                 dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                        'QUAL': str, 'FILTER': str, 'INFO': str},
                 sep='\t').rename(columns={'#CHROM': 'CHROM'})

          print("Extracted column titles successfully!")

        line_list.append(l[0:15])
        del l

  print("Extracted data as string type successfully!")
  gc.collect()

  # Step 7: Filtering the specified chromosome from the user's exome data
  from filter_via_chr_func import filter_via_chr_func

  # Step 8: Storing POS data from vcf file into a string
  from read_POS_vcf import read_POS_vcf

  # Step 9: Executing read_POS_vcf function
  #POS_list = read_POS_vcf([line[0:14] for line in line_list], line_list)
  POS_list = read_POS_vcf([line[0:14] for line in line_list], line_list)
  gc.collect()

  print("POS have been extracted successfully!")

  # Step 10: Replacing first item of numpy array to "POS"
  POS_list[0] = "1kG POS"

  # Step 11: Creating index for POS_list array
  # Initializing an empty list to store the indices
  index_POS_list = []

  for i in range(len(POS_list)):
    index_POS_list.append(i)

  # Converting the list to a numpy array
  index_POS_list = np.array(index_POS_list)

  # Step 12: Concatenate the POS_list and index_POS_list arrays rowwise and transpose the resultant array
  POS_list = np.concatenate([[index_POS_list], [POS_list]]).T

  # Step 13: Filtering the DataFrame for data concerning chromosome 22
  df2 = filter_via_chr_func(chr, user_file_path)
  gc.collect()
  print("Data from user's file has been extracted successfully!")

  # Step 14: Extracting POS from user's exome data and converting DataFrame to numpy array
  JH_POS_list = np.unique(df2.iloc[:, 1].to_numpy())

  # Step 15: Converting the numpy array to a set
  JH_POS_list = JH_POS_list.tolist()
  JH_POS_list = list(map(str, JH_POS_list))

  # Step 16: Determining indices that denote data to be extracted from 1000 Genome vcf file
  print("Comparing 1000 Genome data with user's data...")

  # Initialising indices that denote data to be extracted from 1000 Genome vcf file
  index_list = [list(map(int, np.where(POS_list == POS)[0])) for POS in JH_POS_list]

  # Step 17: Removing empty indices
  index_list = [index for index in index_list if index != []]
  gc.collect()

  index_to_keep_list = [index[0] for index in index_list]

  print("Indices of matching POS from 1000 Genome data have been extracted successfully!")

  del POS_list
  del JH_POS_list
  del line_list
  gc.collect()


  # Step 18: Creating a new list that stores data corresponding to indices in index_list_new
  updated_line_list = []
  count = 0
  index_Thousand_Genome_list = []
  common_index_list = []
  gc.collect()

  with gzip.open(BytesIO(response.content), 'rt') as f:
      for index_Thousand_Genome, line_Thousand_Genome in enumerate(f):
        if not line_Thousand_Genome.startswith('##'):
          count += 1
          del line_Thousand_Genome
          index_Thousand_Genome_list.append(index_Thousand_Genome)

      if count == batch_size:
        common_index_list += list(set(index_Thousand_Genome_list).intersection(index_to_keep_list))
    
      else:
        if count == batch_size + 1:
          count = 0
          index_Thousand_Genome_list = []
          gc.collect()
          
          count += 1
          index_Thousand_Genome_list.append(index_Thousand_Genome)
    
        else:
          common_index_list += list(set(index_Thousand_Genome_list).intersection(index_to_keep_list))

  with gzip.open(BytesIO(response.content), 'rt') as f:
        updated_line_list = [line_Thousand_Genome for i, line_Thousand_Genome in enumerate(f) if not line_Thousand_Genome.startswith('##') and (i in common_index_list) == True]


  print(len(updated_line_list), "out of", len(index_list) ,"data added successfully!")

  # Step 19: Converting the list of strings into a DataFrame
  df2 = pd.read_csv(io.StringIO(''.join(updated_line_list)), sep="\t", header=None)
  df2.columns = df.columns

  # Step 20: Saving the DataFrame to a csv file with title that is dynamically based on chr number
  csv_title = "Thousand_Genome_chr" + str(chr) + ".csv"
  exec(f'df2.to_csv("Thousand_Genome_chr{chr}.csv", index=False)')

  print("Code has run successfully!")

  # Deleting updated_line_list to save memory
  del updated_line_list

  print("==============================================================================================================")

  return csv_title