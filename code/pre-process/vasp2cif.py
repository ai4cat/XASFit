import os
from pymatgen.core import Structure

input_dir = "INPUT_PATH"      
output_dir = "OUTPUT_PATH"     

os.makedirs(output_dir, exist_ok=True)

for file in os.listdir(input_dir):
    if file.endswith(".vasp"):
        vasp_path = os.path.join(input_dir, file)
        cif_name = os.path.splitext(file)[0] + ".cif"
        cif_path = os.path.join(output_dir, cif_name)

        try:
            structure = Structure.from_file(vasp_path)
            structure.to(fmt="cif", filename=cif_path)
            print(f"Done : {file} → {cif_name}")
        except Exception as e:
            print(f"Failed : {file}, {e}")

print("\nFinished!")
