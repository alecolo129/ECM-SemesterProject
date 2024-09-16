# ECM-SemesterProject
This folder contains all the code developed for Re-Discovering ECM Semester project. 

## Files in this folder


* `sage/`
    
    -  `ecm.sage`—Implementation of the ECM algorithm.
    -  `williams.sage`—Implementation of the Williams p+1 algorithm.

    -  `ecm_anomalous.sage`—Implementation of ECM on anomalous curves with symbolic j-invariants.
    
    - `generate.sage`—Script to generate RSA moduli of the form 4p=D(2m+1)^2 + 1.

    - `parameters_supersingular.sage`—Candidate supersingular curves and vulnerable RSA moduli. 

    - `parameters_anomalous.sage`—Candidate anomalous curves and vulnerable RSA moduli. 

    - `test_anomalous.sage`—Tests ECM on anomalous curves.
    - `test_supersingular.sage`—Tests ECM on supersingular curves comparing it to Williams p+1.

* `C/`
    
    -  `curves.h`—Contains curves and points representations.

    -  `ecm.h`—Implementation of the ECM algorithm and ECM with anomalous curves

    -  `test.c`—Test for the implementation of ECM with anomalous curves
    
    - `Makefile`—To compile everithing.


## Setting up the environment

To execute the C implementation you need to install the GMP multiprecision library. All the instructions can be found [here](https://gmplib.org). Note that the library is not fully supported by Apple M1 chips.


## Sage code

We have two test cases for RSA moduli with different vulnerabilities:

1. Each RSA modulus `n=pq` is such that `p+1` is B-powersmooth for `B=2^10`.
2. Each RSA modulus `n=pq` is such that `4p=1+D(2m+1)^2`, for `D<1000`.

The methods that exploit the vulnerability are ECM with supersingular and anomalous curves respectively.

### ECM with Anomalous Curves
Here I compared the implementation of ECM with anomalous curves and Williams p+1. As I reported in the report, running ECM over supersingular curves transforms it into a p+1 method.

Inside `sage` directory, you can run the tests with:
```
sage test_supersingular.sage
```
We clearly see that the Williams p+1 method is way faster than ECM with supersingular curves.

### ECM with Anomalous Curves
Here we tested that the method proposed in the report successfully factors RSA moduli satisfying `4p=1+D(2m+1)^2`. You can generate and test other vulnerable moduli of specified length by using the `generate.sage` file (inside the `sage` directory):
```
sage generate.sage <bits>
```
It is also possible to specify the discriminant `-D` satisfying the condition:
```
sage generate.sage <bits> <-D>
```
If `-D` is not specified, it will be included in `-3, ..., -999`.

We can run the test for ECM on curves with known and unkwnown j-invariants using the `test_anomalous.sage` file:
```
sage test_anomalous.sage
```



# C code

Inside the `C` directory, to compile and test the anomalous implementation, execute the following (after having installed the GMP library):
```
make
./test
```
As this is just a proof of concept, here we just test vulnerable moduli such that  `-D=11`. However, by adding the curves in `parameters_anomalous.sage`,the code can be executed on every `-D` belonging to `{-3, -11, -19, -27, -43, -67, -163}`.
