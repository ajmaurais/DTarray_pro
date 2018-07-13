groff -mandoc > temp.ps `man -w ./db/helpFile.man`
ps2pdf temp.ps helpFile.pdf
rm temp.ps
