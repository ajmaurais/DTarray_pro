# DTarray_pro
## Combine DTASelect-filter files.
DTarray_pro extracts Uniprot ID numbers, molecular weights, and spectral counts from .dtafilter files  stored  in  the  working directory.  Protein  data  is  combined  into  one  dataset and written to the working directory as a tab delimitated text file (.tsv). At runtime, a file list named dtarray_pro_flist.txt, containing all the valid DTAfilter  files  found  in  the  working directory,  is generated. If a file list already exists, the existing file list is used. Various options are available to, modify how files are treated in the working directory, extract additional data from the DTAfilter files, and to change the  output file format.

## Installation

DTarray_pro is written in c++ 11 and can be compiled on linux or OS-X with gcc.  Command line developer tools are required to build DTarray_pro on OS-X.  

DTarray_pro expects to be installed in `~/local`.  The program needs data stored in text files in `~local/<program_dir>/db` for some features to work.  First make `~/local` if it doesn't exist.  
```bash
mkdir -p ~/local
```
Next, download tarball of the latest release to `~/local`.  If you do not wish to build the program in `~/local`, you can run the configure script with the `--progdir <path_to_dir>` option. 

To build the executable, run:
```bash
cd <program_dirrectory>
./configure
make
```
Optionally you can install the executable to `/usr/local/bin`  with.
```bash
make install
```
Or, can create an alias for `DTarray` to `~/local/<program_directory>/bin/DTarray`.  

## Uninstallation
```bash
make uninstall
```

## Usage
From folder containing DTAFilter files with the format `<sample_name>.dtafilter`.  
```bash
DTarray
```
See `~/local/<program_directory>/helpFile.pdf` for detailed explanation of optional arguments or use `DTarray -h` to see help file from the terminal.  
