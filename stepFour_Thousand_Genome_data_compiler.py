def Thousand_Genome_data_compiler(user_data_df, df_10_list, user_data_population):

  print("=====================================================================================================================")
  print("Compiling 1000 Genome data...")

  """
  This function compiles the processed data from:
  (1) 1000 Genome, and
  (2) User's exome data.

  """

  import numpy as np
  import pandas as pd
  import gc
  gc.collect()

  # Step 1: Creating a chromosome number list
  chr_num_list = []

  for chr_num in range(1, 23):
    chr_num_list.append(chr_num)

    if chr_num == 22:
      chr_num_list = list(map(str, chr_num_list))
      chr_num_list.append("X")

  # Step 3: Removing "chr" from chromosome number
  for chr_num in chr_num_list:
    chrom = "chr"+str(chr_num)

  user_data_df = user_data_df.replace(to_replace=chrom, value=chr_num)

  # Step 4: Replacing Jehai sample IDs with "JH" acronym

  user_data_df2 = user_data_df
  column_header_list = list(user_data_df.columns)[0:4]

  for column_header in list(user_data_df.columns)[4:]:
    if column_header[0:2] == user_data_population:
      column_header = column_header[0:2]

    column_header_list.append(column_header)

  user_data_df2.columns = column_header_list

  # Step 5: Establishing data for each population
  df_AFR = df_10_list[0]
  df_AMR = df_10_list[1]
  df_EUR = df_10_list[2]
  df_EAS = df_10_list[3]
  df_SAS = df_10_list[4]

  df_thousandG = [df_AFR, df_AMR, df_EAS, df_EUR, df_SAS]

  # Step 6: Defining a function that finds matching POS between SNPs from compiled Thousand Genome data and user's exome data
  def matching_POS_extractor(df, user_data_df2):

    import gc
    gc.collect()

    POS_ThousandG = np.array(df.iloc[:,1])
    POS_JH = np.array(user_data_df2.iloc[:,1])

    matching_POS = set(POS_ThousandG).intersection(set(POS_JH))

    index_matching_POS_list_JH = [list(POS_JH).index(POS) for POS in matching_POS]
    index_matching_POS_list_df = [list(POS_ThousandG).index(POS) for POS in matching_POS]

    np_POS = np.zeros([user_data_df2.shape[0], df.shape[1] - 2]).astype(int)
    np_POS = np.where(np_POS == 0, "8", 0)

    for indexOne, POS in enumerate(index_matching_POS_list_JH):
      np_POS[index_matching_POS_list_JH[indexOne]] = np.array([df.iloc[index_matching_POS_list_df[indexOne], 2:]])[0]

    df_POS = pd.DataFrame(np_POS, columns = df.columns[2:])

    return df_POS
  
  df_POS_list =[]
  count = 0

  # Executing the helper function
  for df in df_thousandG:
    df_POS = matching_POS_extractor(df, user_data_df2)
    df_POS_list.append(df_POS)
    count += 1
    print("Matching POS with SNPs from user's exome data have been extracted for", count, "/ 5 populations.")
  
  gc.collect()

  # Concatenating SNP data of matching POS with user's exome data from all population columnwise
  df_final = pd.concat(df_POS_list, axis=1)

  # Deleting df_thousandG to clear up some space
  del df_thousandG
  del df_POS_list
  del df_AFR
  del df_AMR
  del df_EUR
  del df_EAS
  del df_SAS
  gc.collect()

  # Encoding genotypes
  df_final = df_final.replace("0|0", 2)
  df_final = df_final.replace("0/0", 2)
  df_final = df_final.replace("0|1", 1)
  df_final = df_final.replace("0/1", 1)
  df_final = df_final.replace("1|0", 1)
  df_final = df_final.replace("1/0", 1)
  df_final = df_final.replace("1|1", 0)
  df_final = df_final.replace("1/1", 0)
  gc.collect()

  # Saving the dataframe to a csv file
  df_final.to_csv('1000_Genome_Final.csv', index=False)

  # Concatenating user's exome data with Thousand Genome data
  df_final = pd.concat([user_data_df2, df_final], axis=1)
  gc.collect()

  # Deleting SNPs that do not have matching POS between user's exome data and Thousand Genome data
  index_delete = [index for index in range(df_final.shape[0]) if ((df_final.iloc[index, :].isin(["8"])).any()) == True]
  df_final = df_final.drop(index_delete)
  gc.collect()

  print(len(index_delete), "SNPs are removed as their POS don't match between user's exome data and Thousand Genome data.")

  #df_final.to_csv('1000_Genome_Final_v2.csv', index=False)
  df_final_transposed = df_final.T
  df_final_transposed = pd.concat([df_final_transposed.reset_index(drop=True), pd.DataFrame(np.array(df_final.columns), columns = ["Population"])], axis=1)
  df_final_transposed.to_csv('1000_Genome_Final_v3.csv', index=False)

  print("1000 Genome data compiled successfully!")
  print("=====================================================================================================================")

  return df_final_transposed