def Thousand_Genome_data_processor(SNP_extracted_list, ID_to_population_file_list, LD_threshold):
  print("=====================================================================================================================")
  print("Processing 1000 Genome data...")

  # Step 1: Importing dependencies
  import numpy as np
  import pandas as pd

  # Step 2: Storing data from csv file for chromosome 1 in a DataFrame
  df = pd.read_csv(SNP_extracted_list[0])

  # Step 3: Expanding the DataFrame to contain data from all chromosomes
  for csv_title in SNP_extracted_list[1:]:
    df = pd.concat([df, pd.read_csv(csv_title)])

  # Step 4: Printing total number of SNPs extracted for all chromosomes and total number of samples
  print("After compiling data, there are", df.shape[0], "SNPs involving", df.shape[1], "individuals (samples).")

  # Step 5: Replacing the matching sample IDs with their corresponding population
  df2 = df
  sample_ID_vcf_list = df.columns[9:]
  count = 0

  for path in ID_to_population_file_list:

    # Storing the data from the csv file in a Pandas DataFrame
    df_ID = pd.read_csv(path)

    # Extracting the population acronym from the csv file path
    population = path[9:12]

    # Removing the gender identifier within the sample IDs
    sample_ID_list = []

    for index in range(df_ID.shape[0]):
      sample_ID = df_ID["Sample (Male/Female/Unknown)"][index][:-4]
      sample_ID_list.append(sample_ID)

    # Replacing the matching sample IDs with their corresponding population
    for ID in sample_ID_list:
      for item in sample_ID_vcf_list:
        replace = lambda item: population if item == ID else item
        sample_ID_vcf_list = list(map(replace, sample_ID_vcf_list))

    df2.columns = list(df2.columns[0:9]) + sample_ID_vcf_list
    sample_ID_vcf_list = df2.columns[9:]
    count += 1
    print("- Number of populations that have undergone sample ID encoding: ", count, "/", len(ID_to_population_file_list))

  # Step 6: Saving the DataFrame as a csv file
  df2.to_csv("1000_Genome_Data_Processed_Step1.csv")

  # Step 7: Dropping unnecessary columns from the DataFrame
  df3 = df2
  df3 = df3.drop(columns = "ID")
  df3 = df3.drop(columns = "QUAL")
  df3 = df3.drop(columns = "INFO")
  df3 = df3.drop(columns = "FILTER")
  df3 = df3.drop(columns = "FORMAT")

  # (COMMENT OUT) Step 9: Importing csv file data and storing it in DataFrame
  #df3 = pd.read_csv("1000_Genome_Data_Processed_Step2.csv")

  # (COMMENT OUT) Dropping the first column
  #df3.drop(columns=df3.columns[0], axis=1, inplace=True)

  # Removing the extra decimals in header of DataFrame for the population acronyms
  df4 = df3

  column_header_list = list(df3.columns)[0:4]

  for column_header in list(df3.columns)[4:]:
    if column_header[0:3] == "AFR" or column_header[0:3] == "AMR" or column_header[0:3] == "EUR" or column_header[0:3] == "EAS" or column_header[0:3] == "SAS":
      column_header = column_header[0:3]

    column_header_list.append(column_header)

  df4.columns = column_header_list

  # Step 10: Dropping columns that are of unknown population
  for column_header in list(df4.columns)[4:]:
    if len(column_header) > 3:
      df4 = df4.drop(columns = column_header, axis=1)

  # Step 11: Removing rows of data that signify INDEL
  indel_list_one = [item for item in df4["REF"].unique() if len(item) > 1]
  indel_list_two = [item for item in df4["ALT"].unique() if len(item) > 1]

  # Obtaining indices of rows that need to be removed
  remove_index_list = []
  for item in indel_list_one:
    remove_index_list += df4.index[df4['REF'] == item].tolist()

  for item in indel_list_two:
    remove_index_list += df4.index[df4['ALT'] == item].tolist()

  # Obtaining unique indices of rows that need to be removed
  from functools import reduce

  def unique(list1):

    # Print directly by using * symbol
    ans = reduce(lambda re, x: re+[x] if x not in re else re, list1, [])

    return ans

  remove_index_list_unique = unique(remove_index_list)

  # Removing the rows of DataFrame based on the indices
  df5 = df4
  df5 = df5.drop(remove_index_list_unique, axis=0)

  print("After removing INDELs, there are", df5.shape[0], "number of SNPs for", df5.shape[1], "individuals (samples).")

  # Removing the extra decimals in header of DataFrame for the population acronyms
  df6 = df5

  column_header_list = list(df5.columns)[0:4]

  for column_header in list(df5.columns)[4:]:
    if column_header[0:3] == "AFR" or column_header[0:3] == "AMR" or column_header[0:3] == "EUR" or column_header[0:3] == "EAS" or column_header[0:3] == "SAS":
      column_header = column_header[0:3]

    column_header_list.append(column_header)

  df6.columns = column_header_list

  df7 = df6

  # Dropping the first column
  df7 = df7.drop(columns=[df7.columns[0], df7.columns[1], df7.columns[2], df7.columns[3]], axis=1)
  df7 = df7.sort_index(axis=1)

  df8 = df6.iloc[:, 0:4].join(df7)

  # Removing the extra decimals in header of DataFrame for the population acronyms
  df9 = df8

  column_header_list = list(df8.columns)[0:4]

  for column_header in list(df8.columns)[4:]:
    if column_header[0:3] == "AFR" or column_header[0:3] == "AMR" or column_header[0:3] == "EUR" or column_header[0:3] == "EAS" or column_header[0:3] == "SAS":
      column_header = column_header[0:3]

    column_header_list.append(column_header)

  df9.columns = column_header_list

  # Step 12: Executing the helper function to filter data based on LD
  from LD_calculator import LD_calculator
  df_10_list = LD_calculator(df9, LD_threshold)

  population_list = ["AFR", "AMR", "EUR", "EAS", "SAS"]

  for index, population in enumerate(population_list):
    csv_title = '1000_Genome_Data_Processed_Step5_' + population + '.csv'
    df_10_list[index].to_csv(csv_title)
    print("Done saving processed data for", population, "population!")
  print("Code has run successfully!")

  print("=====================================================================================================================")

  return df_10_list