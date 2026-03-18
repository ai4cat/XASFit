import os
import shutil

class XmuProcessor:
    """
    1. Collect all xmu.dat files from source folders into a destination folder.
    2. Extract lines after a specific marker in each .dat file.
    """
    def __init__(self, src_path: str, dst_path: str, marker: str):
        self.src_path = src_path
        self.dst_path = dst_path
        self.marker = marker

    def collect(self):
        """Locate all xum.dat files in every folder under the source directory, 
           and copy them into a single destination folder. 
           Each file is renamed to reflect its original source path, 
           ensuring uniqueness while preserving location information."""
        for root, dirs, files in os.walk(self.src_path):
            if 'xmu.dat' in files:
                # calculate relative path to preserve structure
                rel_path = os.path.relpath(root, self.src_path)
                target_dir = os.path.join(self.dst_path, rel_path)
                os.makedirs(target_dir, exist_ok=True)

                src_file = os.path.join(root, 'xmu.dat')
                dst_file = os.path.join(target_dir, 'xmu.dat')

                shutil.copy2(src_file, dst_file)
                #print(f'copy done: {src_file} -> {dst_file}')

    def extract_data(self):
        """Extract lines after marker for each .dat file in the destination folder,
           keep col. 0 and 3."""
        for dirpath, dirnames, filenames in os.walk(self.dst_path):
            for file in filenames:
                if file.endswith('xmu.dat'):
                    file_path = os.path.join(dirpath, file)
                    with open(file_path, 'r') as f:
                        lines = f.readlines()
                    try:
                        idx = lines.index(self.marker + '\n')
                        data_lines = lines[idx+1:]
                    except ValueError:
                        print(f"marker not found in {file_path}")
                        continue

                    # find the x and y
                    processed_lines = []
                    for line in data_lines:
                        parts = line.strip().split()  
                        if len(parts) >= 4:
                            processed_lines.append(f"{parts[0]} {parts[3]}\n")

                    base_name, ext = os.path.splitext(file)
                    save_name = f"{base_name}_extracted{ext}"
                    save_path = os.path.join(dirpath, save_name)

                    with open(save_path, 'w') as f_out:
                        f_out.writelines(processed_lines)

                    #print(f"extracted from {file_path} -> {save_path}")