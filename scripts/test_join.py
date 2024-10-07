import os
import shutil
from astropy.table import Table
import numpy as np
import tempfile
import pytest
import join_vhs_into8  # Ensure your script is importable

def create_dummy_fits(file_path, num_rows=10):
    """Create a dummy FITS file with random data."""
    data = np.random.rand(num_rows, 3)  # Change dimensions as needed
    table = Table(data, names=['a', 'b', 'c'])
    table.write(file_path, overwrite=True)

def verify_concatenated_files(temp_dir, num_expected_files=8):
    """Verify the concatenated files contain all the original data."""
    concatenated_data = []
    for i in range(1, num_expected_files + 1):
        file_path = os.path.join(temp_dir, f"concatenated_{i}.fits")
        assert os.path.exists(file_path), f"{file_path} does not exist"
        concatenated_data.append(Table.read(file_path))
    concatenated_table = Table(np.concatenate([t.as_array() for t in concatenated_data]))
    return concatenated_table

def test_concatenation_of_fits():
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create dummy FITS files
        num_dummy_files = 100  # Adjust based on your needs
        for i in range(num_dummy_files):
            create_dummy_fits(os.path.join(temp_dir, f"dummy_{i}.fits"))

        # Run your concatenation script
        join_vhs_into8.main(temp_dir)

        # Verify the results
        concatenated_table = verify_concatenated_files(temp_dir)
        # Perform additional checks as necessary, e.g., on the size of the concatenated_table
        assert len(concatenated_table) == num_dummy_files * 10  # Adjust based on your dummy data creation

# Run the test
if __name__ == "__main__":
    pytest.main()
