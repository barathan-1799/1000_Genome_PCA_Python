def Thousand_Genome_data_processor(SNP_extracted_list, ID_to_population_file_list, LD_threshold):
  print("=====================================================================================================================")
  print("Processing 1000 Genome data...")

  # Step 1: Importing dependencies
  import numpy as np
  import pandas as pd
  import gc

  # Step 2: Storing data from csv file for chromosome 1 in a DataFrame
  df = pd.read_csv(SNP_extracted_list[0])

  # Step 3: Expanding the DataFrame to contain data from all chromosomes
  for csv_title in SNP_extracted_list[1:]:
    df = pd.concat([df, pd.read_csv(csv_title)])

  # Step 4: Printing total number of SNPs extracted for all chromosomes and total number of samples
  print("After compiling data, there are", df.shape[0], "SNPs involving", df.shape[1] - 9, "individuals (samples).")

  # Step 5: Replacing the matching sample IDs with their corresponding population abbreviations
  df2 = df
  sample_ID_vcf_list = df.columns[9:]
  count = 0

  del df
  gc.collect()

  for path in ID_to_population_file_list:

    # Storing the data from the csv file in a Pandas DataFrame
    df_ID = pd.read_csv(path)

    # Extracting the population acronym from the csv file path
    population = path[9:12]

    # Removing the gender identifier within the sample IDs
    sample_ID_list = []

    for index in range(pd.read_csv(path).shape[0]):
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

    del df_ID
    del sample_ID
    del sample_ID_list
    gc.collect()

    print("- Number of populations that have undergone sample ID encoding: ", count, "/", len(ID_to_population_file_list))

  print("After sample ID encoding, there are", df2.shape[0], "SNPs involving", df2.shape[1] - 9, "individuals (samples).")

  # Step 6: Saving the DataFrame as a csv file
  df2.to_csv("1000_Genome_Data_Processed_Step1.csv")
  gc.collect()


  df3 = pd.read_csv("1000_Genome_Data_Processed_Step1.csv")
  df3.drop(columns=df3.columns[0], axis=1, inplace=True)

  # Step 7: Dropping unnecessary columns from the DataFrame
  df3 = df3.drop(columns = "ID")
  df3 = df3.drop(columns = "QUAL")
  df3 = df3.drop(columns = "INFO")
  df3 = df3.drop(columns = "FILTER")
  df3 = df3.drop(columns = "FORMAT")

  print("After dropping ID, QUAL, INFO, FILTER, and FORMAT columns, there are", df3.shape[0], "SNPs involving", df3.shape[1] - 4, "individuals (samples).")
  

  # Step 8: Dropping columns that are of unknown population
  del df2
  gc.collect()

  df4 = df3
  del df3
  gc.collect()

  # Step 9: Removing rows of data that signify INDEL
  remove_index_list_REF = [(list(df4["REF"]).index(item)) for item in list(df4["REF"]) if len(item) != 1]
  remove_index_list_ALT = [(list(df4["ALT"]).index(item)) for item in list(df4["ALT"]) if len(item) != 1]

  remove_index_list_unique = list(pd.Series(remove_index_list_REF + remove_index_list_ALT).unique())

  print("Number of SNPs of INDEL type:", len(remove_index_list_unique))

  # Removing the rows of DataFrame based on the indices of rows containing INDELs
  df5 = df4
  df5 = df5.drop(remove_index_list_unique, axis=0)
  print("After removing INDELs, there are", df5.shape[0], "number of SNPs involving", df5.shape[1] - 4, "individuals (samples).")

  del df4
  gc.collect()

  # Removing the extra decimals in header of DataFrame for the population acronyms
  column_header_list = list(df5.columns)[0:4]

  for column_header in list(df5.columns)[4:]:
    if column_header[0:3] == "AFR" or column_header[0:3] == "AMR" or column_header[0:3] == "EUR" or column_header[0:3] == "EAS" or column_header[0:3] == "SAS":
      column_header = column_header[0:3]

    column_header_list.append(column_header)

  df5.columns = column_header_list
  print("Removed the extra decimals from the population acronyms!")
  
  # Step 10: Removing unknown population IDs
  for column_header in list(df5.columns)[4:]:
    if len(column_header) > 3:
      df5 = df5.drop(columns = column_header, axis=1)

  print("After removing samples of unknown populations, there are", df5.shape[0], "SNPs involving", df5.shape[1] - 4, "individuals (samples).")

  # Step 10: Rearranging the population alphabetically
  df5_a = pd.concat([df5["CHROM"], df5["POS"], df5["REF"], df5["ALT"]], axis=1)
  df5_b = df5.drop(columns = df5_a.columns)
  df5_b = df5_b.sort_index(axis=1)
  #df6 = pd.concat([df5_a, df5.drop(columns = df5_a.columns).sort_index(axis=1)], axis=1).reset_index(drop=True)
  df6 = pd.concat([df5_a, df5_b], axis=1).reset_index(drop=True)

  
  del df5
  del df5_a
  del df5_b
  gc.collect()
  print("Before performing LD filtering, there are", df6.shape[0], "number of SNPs involving", df6.shape[1] - 4, "individuals (samples).")

  # Step 11: Filtering data based on LD value
  from LD_calculator_v2 import LD_calculator

  print("Performing LD filtering...")
  df_10_list = LD_calculator(df6, LD_threshold)

  population_list = ["AFR", "AMR", "EUR", "EAS", "SAS"]

  for index, population in enumerate(population_list):
    csv_title = '1000_Genome_Data_Processed_Step5_' + population + '.csv'
    df_10_list[index].to_csv(csv_title)
    print("Done saving processed data for", population, "population!")
  print("Code has run successfully!")

  print("=====================================================================================================================")

  return df_10_list