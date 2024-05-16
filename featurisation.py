import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

exec(open("Data_ext.py").read())
df1 = globals()["df"]
df2=globals()["result"]
df2.columns=["Filename","natom", "nelectrons", "denergy", "smiley code", 'C',"H","S","N","O","Misc","numberofrings"]
df1['Modified_File_name'] = df1['Filename'].str.split('.').str[0]
df2['Modified_File_name'] = df2['Filename'].str.split('.').str[0]

# Join the data frames based on the modified File_name
merged_df = pd.merge(df1, df2, on='Modified_File_name', how='inner')

# Drop the temporary 'Modified_File_name' column
merged_df.drop(columns=['Modified_File_name',"smiley code",'numberofrings'], inplace=True)
merged_df.dropna(subset=['denergy'],inplace=True)
merged_df.reset_index()

# Function to generate fingerprint from SMILES code
def generate_fingerprint(smiles_code):
    mol = Chem.MolFromSmiles(smiles_code)
    if mol is not None:
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        # Convert fingerprint to numpy array
        features = np.array(fingerprint)
        return features
    else:
        return None

# Example DataFrame 'merged_df' with 'SMILES' column
# merged_df = pd.read_csv("your_file.csv")  # Load your DataFrame from CSV or any other source

# Apply generate_fingerprint function to 'SMILES' column
merged_df['Fingerprint'] = merged_df['SMILES'].apply(generate_fingerprint)
num_bits = merged_df['Fingerprint'].iloc[0].shape[0]

# Create a DataFrame to store individual features
features_df = pd.DataFrame(merged_df['Fingerprint'].tolist(), columns=[f'Feature_{i+1}' for i in range(num_bits)])
merged_df = pd.concat([merged_df, features_df], axis=1)