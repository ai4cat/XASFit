from ase.io import read, write
import os

# Convert all CONTCAR files in current folder (and subfolders) to CIF
for root, dirs, files in os.walk("PATH"):
    for file in files:
        if file.endswith("CONTCAR") or file.endswith("POSCAR"):
            filepath = os.path.join(root, file)
            try:
                atoms = read(filepath, format="vasp")
                cif_name = os.path.splitext(file)[0] + ".cif"
                output_path = os.path.join("OUT_PATH", cif_name)
                write(output_path, atoms, format="cif")
                print(f"[OK] Converted: {filepath} → {output_path}")
            except Exception as e:
                print(f"[ERROR] {filepath}: {e}")

