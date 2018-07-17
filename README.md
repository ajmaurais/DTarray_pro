# DTarray_pro
## Combine DTASelect-filter files.
DTarray_pro extracts Uniprot ID numbers, molecular weights, and spectral counts from .dtafilter files  stored  in  the  working
directory.  Protein  data  is  combined  into  one  dataset and written to the working directory as a tab delimitated text file
(.tsv). At runtime, a file list named dtarray_pro_flist.txt, containing all the valid DTAfilter  files  found  in  the  working
directory,  is generated. If a file list already exists, the existing file list is used. Various options are available to, mod-
ify how files are treated in the working directory, extract additional data from the DTAfilter files, and to change the  output
file format.

## Installation

DTarray_pro expects to be installed in `~/local`.  First make dirrectory if it dosen't exist.  
```bash
mkdir -p ~/local
```
Download tarball from latest release to `~/local`.
```bash
cd ~/local
bash
curl -L $(curl -s https://api.github.com/repos/ajmaurais/DTarray_pro/releases|grep tarball_url|sed s/'\s*"tarball_url": "'//|sed s/'",$'//) -o DTarray.tar.gz
exit
```

Unpack tarball in current dirrectory.
```bash
tar -zxf DTarray.tar.gz
```
Build executable.
```bash
make
```
Optionally you can install the executable to `/usr/local/bin`  with.
```bash
make install
```
Or if you don't have admin privalages you can create an alias for `DTarray` to `~/local/<program_directory>/bin/DTarray`.


## Uninstallation
```bash
make uninstall
```

## Usage
From folder containnig DTAFilter files with the format `<sample_name>.dtafilter`.  
```bash
DTarray
```
See `~/local/<program_directory>/helpFile.pdf` for detailed explanation of optional arguments or use `DTarray -h` to see help file from the terminal.  
