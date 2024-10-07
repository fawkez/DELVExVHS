# test_mpi.py
from mpi4py import MPI

def main():
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    print(f"Hello from process {rank} out of {size}")

if __name__ == "__main__":
    main()
