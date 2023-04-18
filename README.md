# ECM-SemesterProject
This repository contains the implementation of two methods for RSA modulus factorisation, namely the williams p+1 and ECM methods. These implementations are available both in Sage and in C and can be found in the respective directories.   

## Structure
The code is organized with the following structure:
``` 
c/
  ---> Makefile
  ---> ecm.c
  ---> parameters.h

sage/
  ---> williams.sage
  ---> ecm.sage
  ---> parameters.py
```

## Run
The `sage` files can be just executed using `<filename>.sage` command (SageMath is required). 
The c files can be compiled using the `make` command and run with `./<filename>`.
