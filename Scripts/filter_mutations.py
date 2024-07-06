import pandas as pd
import os


def filter_mutations(strains, ref):

    # Load the csv file
    df = pd.read_csv(f'../Data/variant_analysis/{ref}_allmutations.csv')

    # Print the column names
    print("Variant analysis data imported succesfully")

    # # Create a tuple of the relevant columns for comparison
    df['comparison_key'] = list(zip(df['Minimum'], df['Maximum'], df['product'], df['Change'], df['Amino Acid Change'], df['CDS Position']))


    # Filter the data based on input strains
    df = df[df['Strain'].isin(strains)]

    # Create a dictionary to store wild-type data
    wt_data = {}

    # Iterate through the DataFrame
    for index, row in df.iterrows():
        strain = row['Strain']
        # Check if strain is wild-type
        if strain.startswith('wt'):
            # Get the strain's background
            background = strain.split('_')[1]
            # If the background is not in wt_data, add it
            if background not in wt_data:
                wt_data[background] = []
                print("Found unique mutation")
            # Add the row data to the background list in wt_data
            wt_data[background].append(row['comparison_key'])

    # Create a new DataFrame to store the final data
    new_df = pd.DataFrame()

    # Iterate through the DataFrame again
    for index, row in df.iterrows():
        strain = row['Strain']
        # Check if strain is mutant
        if strain.startswith('mvR'):
            # Get the strain's background
            background = strain.split('_')[1]
            # If the mutation is not in the respective background in wt_data, add it to new_df
            if row['comparison_key'] not in wt_data[background]:
                new_df = pd.concat([new_df, pd.DataFrame(row).T])

    # Drop the 'comparison_key' column as it is not needed in the final csv
    new_df = new_df.drop(columns='comparison_key')

    # Save new_df as a csv file
    new_df.to_csv(f'../Data/variant_analysis/{ref}_filtered_mutations.csv', index=False)

    return 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Filter mutations and plot a circos diagram.')
    parser.add_argument('-strains', nargs='+', help='The list of strains to include in the analysis.')
    parser.add_argument('-ref', help='The reference genome.')
    args = parser.parse_args()

    filter_mutations(args.strains, args.ref)