The MATLAB algorithm consists of the following steps:

Read the data file (read.m): Read the input data file and store the relevant information, including atom, bond, and angle data.

Reorder the data (ordering.m): Sort the atom, bond, and angle information into the required order for subsequent processing.

Sort molecules (Sort_molecules.m): Execute this script to organize the molecular information based on the reordered data.

Random allocation (rand_allo.m): Define the simulation box dimensions and the number of polymer chains to be allocated.
For example: rand_allo([200 200 200], 60).
Here, [200 200 200] specifies the simulation box dimensions along the (x), (y), and (z) directions, respectively, and 60 specifies
the number of polymer chains to be allocated randomly within the box.

Generate the final data file (print_new_data_file.m): Execute this script to generate the final data file containing the processed molecular information.
