import numpy as np
import pandas as pd

# Creating a function that filters data based on LD value
def LD_calculator(df9, LD_threshold):


    # Obtaining DataFrames for CHROM and POS
    df_CHROM_POS =  pd.merge(df9["CHROM"], df9["POS"], right_index = True, left_index = True)

    # Obtaining DataFrames for each population
    df_AFR = df9["AFR"]
    df_AMR = df9["AMR"]
    df_EUR = df9["EUR"]
    df_EAS = df9["EAS"]
    df_SAS = df9["SAS"]

    # Storing the DataFrames for each population in a list
    df_population_list = [df_AFR, df_AMR, df_EUR, df_EAS, df_SAS]

    # Initialising an empty list
    df10_list = []

    for df in df_population_list:

      # Initialising and resetting to empty lists
      n_zero_zero = []
      n_zero_one = []
      n_one_zero = []
      n_one_one = []
      n_total = []

      # Calculating number of genotypes for each genotype type
      n_zero_zero = [df.iloc[index,:].tolist().count("0|0") for index in range(df9.shape[0])]
      n_zero_one = [df.iloc[index,:].tolist().count("0|1") for index in range(df9.shape[0])]
      n_one_zero = [df.iloc[index,:].tolist().count("1|0") for index in range(df9.shape[0])]
      n_one_one = [df.iloc[index,:].tolist().count("1|1") for index in range(df9.shape[0])]

      # Calculating total number of genotypes
      n_total = list(map(sum, zip(n_zero_zero, n_zero_one, n_one_zero, n_one_one)))

      # Converting lists to DataFrame
      df_LD = pd.DataFrame(n_zero_zero, columns = ["Number of 0|0"])
      df_LD = df_LD.join(pd.DataFrame(n_zero_one, columns = ["Number of 0|1"]))
      df_LD = df_LD.join(pd.DataFrame(n_one_zero, columns = ["Number of 1|0"]))
      df_LD = df_LD.join(pd.DataFrame(n_one_one, columns = ["Number of 1|1"]))
      df_LD = df_LD.join(pd.DataFrame(n_total, columns = ["Total number of genotype"]))

      # Calculating probabilities
      Pa = list((df_LD.iloc[:,1] + df_LD.iloc[:,3])/df_LD.iloc[:,-1].unique()[0])
      Pb = list((df_LD.iloc[:,2] + df_LD.iloc[:,3])/df_LD.iloc[:,-1].unique()[0])
      Pab = list(df_LD.iloc[:,3]/df_LD.iloc[:,-1].unique()[0])

      # Calculating LD
      numerator = (np.array(Pab) - np.multiply(np.array(Pa), np.array(Pb)))**2
      denominator_partOne = np.multiply(np.array(Pa), 1 - np.array(Pa))
      denominator_partTwo = np.multiply(np.array(Pb), 1 - np.array(Pb))
      denominator = np.multiply(denominator_partOne, denominator_partTwo)
      LD = np.divide(numerator, denominator)

      # Converting lists to DataFrame and appending to earlier DataFrame
      df_Pa = pd.DataFrame(Pa, columns = ["P_A"])
      df_Pb = pd.DataFrame(Pb, columns = ["P_B"])
      df_Pab = pd.DataFrame(Pab, columns = ["P_AB"])
      df_rsquared = pd.DataFrame(LD, columns = ["LD"])

      df_P_list = [df_Pa, df_Pb, df_Pab, df_rsquared]

      for df_P in df_P_list:
        df_LD = df_LD.join(df_P)

      # Extracting row indices that denote LD > 0.2
      index_remove_LD_one = df_LD.index[(df_LD['LD'] > LD_threshold) == True].tolist()
      index_remove_LD_two = df_LD.index[(df_LD['LD'] == "NaN") == True].tolist()
      index_remove_LD = index_remove_LD_one + index_remove_LD_two

      # Incorporating CHROM and POS data to genotype DataFrame
      df2 = pd.merge(df_CHROM_POS, df, right_index = True, left_index = True)

      # Remove the dataFrame based on index_remove_LD
      df2 = df2.drop(index_remove_LD, axis=0)

      # Add the filtered dataFrame to a list
      df10_list.append(df2)

    return df10_list