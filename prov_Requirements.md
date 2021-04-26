# Requirements 

To install Svhip you need to install:


## Git 
To clone the repository. 

https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

## Newest versions of Pip and Python3 

Python
https://docs.python-guide.org/starting/install3/linux/

Pip
https://www.tecmint.com/install-pip-in-linux/


## GCC C++ compiler 

https://linuxize.com/post/how-to-install-gcc-compiler-on-debian-10/

Or on Mageia 
```bash
sudo uprmi gcc-c++
```


## The Vienna RNAz package

https://github.com/ViennaRNA/RNAz

If it doesn't work like instructed in the repository, you can try using GNU-89 as compiler or try work-around with Conda. 

https://conda.io/projects/conda/en/latest/user-guide/install/index.html 
https://anaconda.org/bioconda/viennarna

## The ViennaRNA package

https://www.tbi.univie.ac.at/RNA/tutorial/

If you find yourself an error like "EXTERN.h not found" or "xlocale.h not found" try 
```bash
ln -s /usr/include/locale.h /usr/include/xlocale.h
```
and redo the installation.

## Biopython
```bash
pip3 install biopython
```

## clustalw2

https://bioinformaticsreview.com/20210126/how-to-install-clustalw2-on-ubuntu/

Note you might need to use sudo make (install). 
