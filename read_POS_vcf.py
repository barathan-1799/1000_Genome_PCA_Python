import numpy as np
import pandas as pd

# Defining a function that stores POS data from vcf file into a string
def read_POS_vcf(shortened_line_list, line_list):

    index_space = [[index for index, character in enumerate(item) if character == "\t"][0:2] for item in [line[0:14] for line in line_list]]

    print("Indices of spaces have been extracted successfully!")

    index_first_space = [item[0] for item in index_space]
    index_second_space = [item[1] for item in index_space]

    print("Extracting POS...")

    # Extracting the POS from each line in vcf file
    # Converting the resultant list to a numpy array
    POS_list = []

    for i in range(len(shortened_line_list)):
      POS = shortened_line_list[i][index_first_space[i] + 1:index_second_space[i]]
      POS_list.append(POS)

    POS_list = np.array(POS_list)

    return POS_list