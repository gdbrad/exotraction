import h5py 
import glob 
import argparse

def merge_hdf5_files(output_file, input_files):
    with h5py.File(output_file, "w") as h5f_out:
        # h5_group_out = h5f_out.create_group('pion_2pt')

        for input_file in input_files:
            print(f"Merging {input_file}...")
            with h5py.File(input_file, "r") as h5f_in:
                h5_group_in = h5f_in['pion_000']

                for tsrc in h5_group_in:
                    tsrc_in = h5_group_in[tsrc]
                    if tsrc not in h5f_out:
                        tsrc_out = h5f_out.create_group(tsrc)
                    else:
                        tsrc_out = h5f_out[tsrc]
                # Copy all datasets and groups from input file to output file
                    for cfg in tsrc_in:
                        if cfg not in tsrc_out:
                            data = tsrc_in[cfg][:]
                            tsrc_out.create_dataset(cfg,data=data)
                        # h5f_in.copy(h5_group_in[cfg], h5_group_out)

    print(f"All files merged into {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="merge batched contractions into a single h5 file.")
    parser.add_argument('--nv',type=int,required=True)
    args = parser.parse_args()

    # List all task-generated HDF5 files (adjust the pattern as necessary)
    input_files = glob.glob(f'pion_2pt_nvec_{args.nv}_tsrc_24_task*.h5')
    
    # Specify the output file
    output_file = f'pion_nv-{args.nv}_tsrc_24_master.h5'

    # Merge the input files into the final HDF5 file
    merge_hdf5_files(output_file, input_files)



