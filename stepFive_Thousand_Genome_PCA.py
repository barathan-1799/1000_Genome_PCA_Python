def Thousand_Genome_PCA(df, user_data_population):

  # Step 1: Importing useful python packages
  import numpy as np
  import pandas as pd
  import matplotlib.pyplot as plt
  import seaborn as sns
  from sklearn.preprocessing import StandardScaler
  from sklearn.decomposition import PCA
  import gc
  gc.collect()

  # Step 2: Extract the values of all features' data
  start_x_index_row = 5
  x = np.array(df.iloc[start_x_index_row:df.shape[0] + 1,:-1]).astype(int)
  x = pd.DataFrame(x).values

  # Step 3: Find standard deviation and mean of the existing data of features
  current_mean = np.mean(x)
  current_std = np.std(x)

  # Step 4: Defining targets of the data
  target_list = pd.DataFrame(np.array(df.iloc[start_x_index_row:df.shape[0] + 1,-1:]).astype(str))

  # Step 5: Scale the data of features such that new mean = 0 and new std = 1
  if current_mean > 1e-10 or current_std > 1:
    x2 = StandardScaler().fit_transform(x)
    new_mean = np.mean(x2)
    new_std = np.std(x2)

  df3 = pd.DataFrame(x2)

  # Adding back model column at the last column
  df4 = pd.concat([df3, target_list], axis=1)
  target = target_list.values

  # Step 6: Conducting PCA for standardised data of feature
  print("=====================================================================================================================")
  print("Conducting 2D PCA...")
  pca = PCA(n_components=2)
  principalComponents = pca.fit_transform(x2)

  principal_df = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])

  print('Explained variability per principal component: {}'.format(pca.explained_variance_ratio_))

  PC1_var_capture = pca.explained_variance_ratio_[0]*100
  PC2_var_capture = pca.explained_variance_ratio_[1]*100
  print("- Thus, PC1 captures " + str(round(PC1_var_capture)) + "%")
  print("- Thus, PC2 captures " + str(round(PC2_var_capture)) + "%")

  # Step 7: Plotting a 2D PCA plot
  plt.figure()
  plt.figure(figsize=(10,10))
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=14)
  PC1_title = "PC1 (" + str(round(PC1_var_capture)) + "%)"
  PC2_title = "PC2 (" + str(round(PC2_var_capture)) + "%)"
  plt.xlabel(PC1_title,fontsize=18, labelpad=20)
  plt.ylabel(PC2_title,fontsize=18, labelpad=20)
  target_legend_list = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS', user_data_population]
  colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

  for target_legend, color in zip(target_legend_list,colors):
    indicesToKeep = target == target_legend
    plt.scatter(principal_df.loc[indicesToKeep, 'principal component 1'],
                principal_df.loc[indicesToKeep, 'principal component 2'], c = color, s = 50)

  plt.legend(target_legend_list,prop={'size': 15})
  plt.savefig("PCA_2D.png")
  print("2D PCA plot has been generated and saved successfully!")

  # Step 8: Conducting PCA for standardised data of feature (3 PCs)
  print("Conducting 3D PCA...")
  pca = PCA(n_components=3)
  principalComponents = pca.fit_transform(x2)

  principal_df = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2', 'principal component 3'])

  print('Explained variability per principal component: {}'.format(pca.explained_variance_ratio_))

  PC1_var_capture = pca.explained_variance_ratio_[0]*100
  PC2_var_capture = pca.explained_variance_ratio_[1]*100
  PC3_var_capture = pca.explained_variance_ratio_[2]*100
  print("- Thus, PC1 captures " + str(round(PC1_var_capture)) + "%")
  print("- Thus, PC2 captures " + str(round(PC2_var_capture)) + "%")
  print("- Thus, PC3 captures " + str(round(PC3_var_capture)) + "%")

  # Step 9: Plotting 3D PCA 
  from matplotlib.colors import ListedColormap
  sns.set_style("whitegrid", {'axes.grid' : False})

  # get colormap from seaborn
  cmap = ListedColormap(sns.color_palette("husl", 256).as_hex())

  plt.figure(figsize=(14, 14))
  ax = plt.axes(projection='3d')
  plt.xticks(fontsize=24)
  plt.yticks(fontsize=24)

  PC1_title = "PC1 (" + str(round(PC1_var_capture)) + "%)"
  PC2_title = "PC2 (" + str(round(PC2_var_capture)) + "%)"
  PC3_title = "PC3 (" + str(round(PC3_var_capture)) + "%)"

  ax.set_xlabel(PC1_title,fontsize=24, labelpad=40)
  ax.set_ylabel(PC2_title,fontsize=24, labelpad=40)
  ax.set_zlabel(PC3_title,fontsize=24, labelpad=40)

  font = {'size': 24}
  ax.tick_params('z', labelsize=font['size'])

  ax.tick_params(axis='both', which='major', pad=20)

  targets = ['AFR', 'AMR', 'EUR', 'EAS', 'SAS', user_data_population]

  for target_legend, color in zip(target_legend_list,colors):
    indicesToKeep = target == target_legend
    ax.scatter3D(principal_df.loc[indicesToKeep, 'principal component 1'],
                principal_df.loc[indicesToKeep, 'principal component 2'],
                principal_df.loc[indicesToKeep, 'principal component 3'], cmap=cmap, s = 80)

  ax.legend(targets,prop={'size': 20}, ncol=6, labelspacing=0.4, columnspacing=1.4, loc='upper center')
  ax.view_init(20, 50)

  import matplotlib.ticker as ticker
  tick_spacing = 40
  ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
  ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
  ax.zaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

  plt.savefig("PCA_3D.png")
  print("3D PCA plot has been generated and saved successfully!")
  print("=====================================================================================================================")