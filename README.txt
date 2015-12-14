Hello! I'm surprised you're reading this!

This code was written by Cem Bagdatlioglu for computational methods class in Fall 15

Once you have copied the git repository to your local linux system, change directory 
to the main folder (where this file is) and run:
$ ./install.py
to install the software. The install script will automatically run the input: Input/input.txt

See the input file template to write your own input file.

To run the code after installation:
$ ./transport

In order to run your input file, use the "ipath" (or simply "i") tag:
$ ./transport ipath Input/your_input_file.txt
-or-
$ ./transport i Input/your_input_file.txt

The command line tags are (input file overrides if there's a conflict):
- groups: energy groups to use
- dpath: database path
- ipath (or i): input file path and name
- f: ffactor order
- s: legendre order
- threads (or t): number of software threads to use (you probably read this far for the threading)

The folders are:
- benchmarks: some saved input files for benchmarking
- Data: this is where the isotope data goes, unless specified otherwise
- Input: this is a place to store input files
- Saved: saved outputs to compare benchmark results
- src: the source code for the software

Go transport some neutrons!