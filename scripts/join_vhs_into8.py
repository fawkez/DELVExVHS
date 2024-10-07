import os
import glob
from astropy.table import Table, vstack
from concurrent.futures import ThreadPoolExecutor

def process_files(file_chunk, output_file):
    tables = []
    for filename in file_chunk:
        try:
            table = Table.read(filename)
            tables.append(table)
        except Exception as e:
            print(f"Error reading {filename}: {e}")
    if tables:
        concatenated_table = vstack(tables)
        concatenated_table.write(output_file, overwrite=True)
        print(f"Written {output_file} successfully.")

def main(input_dir, num_chunks=8):
    files = glob.glob(os.path.join(input_dir, '*.fits'))
    chunk_size = len(files) // num_chunks
    file_chunks = [files[i:i + chunk_size] for i in range(0, len(files), chunk_size)]

    with ThreadPoolExecutor(max_workers=num_chunks) as executor:
        for i, file_chunk in enumerate(file_chunks):
            output_file = f"/fast_scratch3/mncavieres/VHS_merged/vhs_{i+1}.fits"
            executor.submit(process_files, file_chunk, output_file)

if __name__ == "__main__":
    input_dir = '/fast_scratch3/mncavieres/VHS'
    main(input_dir)
