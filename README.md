# DTsetup
## Automatically setup DTASelect-filter to be read into DTarray_pro
`DTsetup` searches all directories one level below where it is run for a file named `DTASelect-filter.txt`. If a matching file is found, it is copied into a folder named `<date>_<time>_dtarraySetup` with the name `<sample_name>.dtafilter`.  Once `DTsetup` is run, [DTarray_pro](https://github.com/ajmaurais/DTarray_pro) can be run from the newly created directory to comapre spectral counts for proteins found in multiple runs.

## Installation
The script is a stand alone `bash` executable no installation is required.

## Usage
From a directory containing subdirectories with DTASelect-filter files
```bash
DTsetup
```
