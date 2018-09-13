# DTarray_pro
## Combine DTASelect-filter files.
DTarray_pro extracts Uniprot ID numbers, molecular weights, and spectral counts from .dtafilter files  stored  in  the  working directory.  Protein  data  is  combined  into  one  dataset and written to the working directory as a tab delimitated text file (.tsv). At runtime, a file list named dtarray_pro_flist.txt, containing all the valid DTAfilter  files  found  in  the  working directory,  is generated. If a file list already exists, the existing file list is used. Various options are available to, modify how files are treated in the working directory, extract additional data from the DTAfilter files, and to change the  output file format.

## Installation

DTarray_pro is written in c++ 98 and can be compiled on linux or OS-X with gcc.  Command line developer tools are required to build DTarray_pro on OS-X.  

DTarray_pro expects to be installed in `~/local`.  The program needs data stored in text files in `~local/<program_dir>/db` for some features to work.  First make `~/local` if it doesn't exist.  
```bash
mkdir -p ~/local
```
Download tarball of the latest release to `~/local`.  The following commands should be executed from a `bash` shell.  To check your shell enviroment use `echo $SHELL`.  If the output is not something like `/bin/bash`, enter `bash` to start a `bash` shell.   The default shell on `pleiades.bc.edu` is `tcsh`.  
```bash
cd ~/local
curl -L $(curl -s https://api.github.com/repos/ajmaurais/DTarray_pro/releases/latest|grep tarball_url|sed s/'\s*"tarball_url": "'//|sed s/'",$'//) -o DTarray.tar.gz
```
When you have downloaded the tarball, you can use `exit` to return to your default shell environment.  

Unpack tarball in current directory.
```bash
tar -zxf DTarray.tar.gz
```
Build executable.
```bash
cd <program_dirrectory>
make
```
Optionally you can install the executable to `/usr/local/bin`  with.
```bash
make install
```
If you don't have admin privileges you can create an alias for `DTarray` to `~/local/<program_directory>/bin/DTarray`.  To create an shortcut for `DTarray_pro`, you will need to add to your shell config file.  For `tcsh`, the shell config file is `~/.tcshrc`.  
```bash
nano ~/.tcshrc
```


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
