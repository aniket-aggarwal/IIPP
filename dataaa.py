### Code to extract  .nwo files from folders and subfolder 
import os
import shutil
def find_files(directory_path):
    for root, dirs, files in os.walk(directory_path):
        for dir in dirs:
            subfolder_path = root + dir
            for sub_root, sub_dirs, sub_files in os.walk(subfolder_path):
                for file in sub_files:
                    file_path = os.path.join(sub_root, file)
                    yield file_path

# Define the directory path
directory_path = r"data\\"

def generate_custom_file_name(file_path):
    # Extract relevant components from the file path
    components = file_path.split("\\")
    mod=components[-5][0:4]
    monomer_type = components[-4].split("_")[-3]  # Extract monomer type from the third last component
    method = components[-3][0:4]  # Extract method from the second last component
    molecule_name = components[-2].split("_")[-1]  # Extract molecule name from the second last component
    product_name = components[-1].split("_")[-1]  # Extract product name from the last component
    
    # Construct the custom file name
    custom_file_name = f"{mod}_{monomer_type}_{method[:3]}_{molecule_name}_{product_name.split('.')[0]}"
    return custom_file_name
destination_directory="xyz files\\"
destination_directory2="smiles2xyz\\"
# Loop through files in subfolders and find those with the specified extension
l=[]
for file_path in find_files(directory_path, file_extension):
    # Check if the file type is "nwo"
    if file_path.endswith(".nwo"):
        # Extract the file name
        file_name = os.path.basename(file_path)
        custom_file_name = generate_custom_file_name(file_path)
        print(custom_file_name)
        print(file_path)
        custom_file_name = custom_file_name+ ".nwo"  # Add the file extension
        destination_path = os.path.join(destination_directory, custom_file_name)
        l.append([custom_file_name,file_path])
        shutil.copy(file_path, destination_path)
    elif file_path.endswith(".xyz"):
        file_name = os.path.basename(file_path)
        custom_file_name = generate_custom_file_name(file_path)
        print(custom_file_name)
        print(file_path)
        custom_file_name = custom_file_name+ ".xyz"  # Add the file extension
        destination_path = os.path.join(destination_directory2, custom_file_name)
        l.append([custom_file_name,file_path])
        shutil.copy(file_path, destination_path)
