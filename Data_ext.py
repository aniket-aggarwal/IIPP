import os
import numpy as np
import pandas as pd
from rdkit import Chem
import xyz2mol as x2m

class DataExtractor:
    def __init__(self, folder_path, repo_path):
        self.folder_path = folder_path
        self.repo_path = repo_path

    def get_files_in_folder(self):
        files = os.listdir(self.folder_path)
        return files

    def read_first_line(self, file_path):
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
        return first_line

    def extract_data(self, folder_path, file_name):
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, 'r',encoding='utf-8') as file:
            lines = file.readlines()
            natom = None
            denergy = None
            nelectrons= None
            for line in lines:
                if "No. of atoms     " in line:
                    natom = line.split()[-1]
                elif "No. of electrons " in line:
                    nelectrons = line.split()[-1] # Extracting 3 characters after "No of atoms"
                elif "Total DFT energy" in line:
                    denergy = line.split()[-1]  # Extracting data after "Total DFT energy"
            data = []
            with open(file_path, 'r', encoding='utf-8') as file:
                lines = file.readlines()
                start_reading = False
                for line in lines:
                    if "No.       Tag          Charge          X              Y              Z" in line:
                        start_reading = True
                        continue
                    if start_reading:
                        parts = line.split()
                        if len(parts) == 6:  # Assuming each line has exactly 6 columns
                            data.append(parts)
                        if len(data) == float(natom)+1:  # Read two lines after encountering the header
                            break
            df = pd.DataFrame(data, columns=["No.", "Tag", "Charge", "X", "Y", "Z"])
            #df.drop(0,inplace=True)

        return df,natom,nelectrons, denergy
    
    def count_number_pairs(self,string):
        stack = []  # Initialize stack
        pair_count = 0  # Initialize count for number pairs

        # Iterate through the characters in the string
        for char in string:
            if char.isdigit():
                if stack and stack[-1] == char:
                    stack.pop()  # Remove top digit from stack
                    pair_count += 1  # Increment count for number pairs
                else:
                    stack.append(char)  # Append current digit to stack

        return pair_count
    def extracter_smiles(self,folder_path_smiles):
            data=[]
            for filename in os.listdir(folder_path_smiles):
                if filename.endswith(".xyz"):  # Assuming files have .txt extension
                    file_path = os.path.join(folder_path_smiles, filename)
                    with open(file_path, 'r') as file:
                        lines = file.readlines()
                        smiles_code = lines[1].strip().split(":")[1]
                        digits_list = []  # Initialize list to store digits
                        pair_count=0
                    # Iterate through the characters in the string
                        for char in str(smiles_code):
                            if char.isdigit():
                                if char in digits_list:
                                    digits_list.remove(char)  # Remove the digit from the list
                                    pair_count += 1  # Increment count for number pairs
                                else:
                                    digits_list.append(char)  # Add the digit to the list

                        num_atoms = int(lines[0].strip())  # Extract number of atoms from the first line
                        data.append({"Filename": filename, "SMILES": smiles_code, "Num_Atoms": num_atoms,"nrings":pair_count})
            return data
    def process_files(self):
        data = np.empty((0, 12))
        file_names = self.get_files_in_folder()

        for file_name in file_names:
            df, natom, nelectrons, denergy = self.extract_data(self.folder_path, file_name)
            da=pd.Series(df.drop('No.',axis=1).value_counts('Tag')).reindex()
            element_order = ['C', 'H', 'S', 'N', 'O', 'Misc']
            for i, element in enumerate(element_order):
                if element is not None:
                    if element in da.index:
                        globals()[element]= da[element]
                    else:
                        globals()[element]=0
            
            #number_pairs_count = self.count_number_pairs(c)
            row = np.array([[file_name,natom, nelectrons, denergy, "smiley code", C,H,S,N,O,Misc,"numberofrings"]])
            data = np.append(data, row, axis=0)
            array = np.zeros((1, len(element_order)), dtype=int)

        return data

    def run(self):
        return self.process_files()

folder_path = r'xyz files'
repo_path = r'git_repo'
extractor = DataExtractor(folder_path, repo_path)
result= pd.DataFrame(extractor.run())
folder_path_smiles = r"C:\Users\Admin\nsysu\IIPP\smiles2xyz\\"
df = pd.DataFrame(DataExtractor.extracter_smiles(DataExtractor,folder_path_smiles))
x=2
