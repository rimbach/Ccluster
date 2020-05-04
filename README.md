# Ccluster

Ccluster is a C library implementing an algorithm for local clustering of the complex roots
of a univariate polynomial whose
coefficients are any complex numbers.

It has a wrapper for Julia (https://github.com/rimbach/Ccluster.jl)
and is available in Singular (https://github.com/Singular).

The inputs of the clustering algorithm are
a polynomial *P*, 
a square complex box *B* and
a rational number *eps*.

It outputs a set of *eps-natural clusters* of roots together with the sum of multiplicities
of the roots in each cluster.
An *eps-cluster* is a complex disc *D* of radius at most *eps* containing at least one root, 
and it is *natural when* *3D* contains the same roots
than *D*.
Each root of *P* in *B* is in exactly one cluster of the output, and clusters may contain
roots of *P* in *2B*.

The implemented algorithm is described here:
https://dl.acm.org/citation.cfm?id=2930939

Please cite:
https://link.springer.com/chapter/10.1007/978-3-319-96418-8_28
if you use it in your research.

## Dependencies

Ccluster depends on flint2 (>= 2.5.2, https://github.com/wbhart/flint2)
and arb (>= 2.17.0, https://github.com/fredrik-johansson/arb).

## Installation

1. clone Ccluster

2. In the Ccluster directory,

```
./configure --with-flint="path to flint" --with-arb="path to arb"
make
```
Specifying the location of flint and arb is optional if these libraries are in /usr/local. 
If you want to install Ccluster in /usr/local and have permissions:
```
make install
```
## Binary ccluster: compute *eps-natural clusters* of complex roots of a polynomial with rational coefficients

We provide a binary file in the directory bin that takes in input a polynomial with rational coefficients, given in an input file.

### Input files

A suitable input for ccluster for the polynomial *1z^3 + (1/2)z^2 + (1/3)z + (1/4)* is 
```
4 1/4 1/3 1/2 1
```
i.e the first string is the length, then the coefficients are given according to increesing degree.

We provide a programm called genPolFile that generates input files for state of the art polynomials. Try command
```
bin/genPolFile
```
to see the possible options. In particular, genPolFile can also generate input files for anewdsc or mpsolve.

### Basic options
Try command
```
bin/ccluster
```
to see the possible options and default values for ccluster.

Let's say you want to compute a set of clusters of size less *2^-53* for all roots of Mignotte polynomial of degree 64 and bitsize 14.
Do:
```
bin/genPolFile Mignotte 1 64 14 bin/Mignotte_64_14.ccl
```
to generate the input file, then 
```
bin/ccluster bin/Mignotte_64_14.ccl -e -53 -o 5
```
means that 5-bit approximations of the computed clusters are printed.
For a graphical output, if gnuplot is install on your system, do
```
bin/ccluster bin/Mignotte_64_14.ccl -o -2 | gnuplot
```
and use option -o -3 to draw the subdivision tree.

## Use with Singular
Requires the ccluster library to be installed on your system.

1. checkout the latest version of Singular

2. compile it 

```
./autogen.sh
./configure --with-flint="path to flint" --with-ccluster="path to ccluster"
make 
make install
```

Specifying the location of flint and ccluster is optional if these libraries are in /usr/local. 

3. Launch Singular

Clustering the roots of a polynomial with real rational coefficients:

```
ring r=QQ,x,dp;
poly f=x^2 + x +1;
list l=system("ccluster",f,0/1,0/1,10/1,1/2^53);
l;
```

Clustering the roots of a polynomial with complex rational coefficients can be done
by passing to Ccluster 2 polynomials: one whose coefficients are the real parts,
one whose coefficients are the imaginary parts:

```
ring r=QQ,x,dp;
poly fre=x^2 + x +1;
poly fim=3*x^2 + 2*x +2;
list l=system("ccluster",fre,fim,0/1,0/1,10/1,1/2^53);
l;
```

Arguments 0/1,0/1,10/1 define the initial box (real part of the center, imaginary part of the center
and width) and 1/2^53 is the value for *eps*.

## Use with Julia

See (https://github.com/rimbach/Ccluster.jl).
