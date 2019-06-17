# DTarray_pro
## Combine DTASelect-filter files.
`DTarray_pro` extracts Uniprot ID numbers, molecular weights, and spectral counts from .dtafilter files  stored  in  the  working directory.  Protein  data  is  combined  into  one  dataset and written to the working directory as a tab delimitated text file (.tsv). At runtime, a file list named `dtarray_pro_flist.txt`, containing all the valid DTASelect-filter files  found  in  the  working directory,  is generated. If a file list already exists, the existing file list is used. Various options are available to, modify how files are treated in the working directory, extract additional data from the DTASelect-filter files, and to change the  output file format.

## Installation

`DTarray_pro` is written in c++ 11 and can be compiled on linux or OS-X with gcc. 

`DTarray_pro` expects to be installed in `~/local`.  The program needs data stored in text files in `~local/<program_dir>/db` for some features to work.  First make `~/local` if it doesn't exist.  
```bash
mkdir -p ~/local
```
If you do not wish to build the program in `~/local`, you can run the configure script from the desired location with the `--progdir ./` option. 

`DTarray_pro` depends on the [utils](https://github.com/ajmaurais/utils) c++ library which is included as a submodule in this repository. To clone this repository and the `utils` submodule, run:

```bash
git clone --recurse-submodules --single-branch --branch stable https://github.com/ajmaurais/DTarray_pro
```

To build the executable, run:
```bash
cd ~/local/DTarray_pro
./configure
make
```

After compilation, the `bin` directory should contain two files:

* `DTarray` - the `DTarray_pro` executable.
* `DTSetup` - a bash script which optionally be used to automatically setup the input directory for `DTarray`

## Usage
From a folder containing DTASelect-filter files with the format `<sample_name>.dtafilter`. 
```bash
DTarray
```

See `~/local/<program_directory>/helpFile.pdf` for a detailed explanation of optional arguments or use `DTarray -h` to see the help file from the terminal.  
