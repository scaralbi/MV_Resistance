import pandas as pd

# Load the csv files into pandas dataframes
df1 = pd.read_csv('Data/variant_analysis/NC000911_allmutations.csv')
df2 = pd.read_csv('Data/genome_annotations/NC_000911_Annotations_Locus_Tags.csv')

# Create a dictionary from file2 where locus_tag is the key and tag is the value
locus_to_tag = df2.set_index('tag')['locus_tag'].to_dict()

# Replace the locus_tag values in file1 based on the dictionary
df1['locus_tag'] = df1['locus_tag'].map(locus_to_tag).fillna(df1['locus_tag'])

# Save the updated dataframe to a new csv file
df1.to_csv('Data/variant_analysis/NC000911_allmutations_locus_tag.csv', index=False)
