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

The `.sage` files can be executed inside the [sage](sage/) directory using `sage <filename>.sage` ([SageMath](https://www.sagemath.org) is required). 
The `.c` files can be compiled inside the [c](c/) directory using the `make` command and executed with `./<filename>` (the [GMP](https://gmplib.org) multiprecision library has to be already installed).
