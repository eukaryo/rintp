# rintp
rintp

# usage

$ make

to compile the source and obtain the executable "rintp". I succeed to compile them in the following environment:

* gcc 8.2.0 on CentOS 7
* gcc 5.4.0 on Ubuntu 16.04

then

$ ./rintp "CCCCAAAAGGGG" "((((....))))" 12 RintPwithDFT

# output

If the input is valid, output is given from Standard Output in the following format:

> R  
> 0 z[0]  
> 1 z[1]  
> :  
> (R-1) z[R-1]  
> 0 0 bulge P[0][0]["bulge"]  
> 0 0 exterior P[0][0]["exterior"]  
> 0 0 hairpin P[0][0]["hairpin"]  
> 0 0 internal P[0][0]["internal"]  
> 0 0 multi P[0][0]["multi"]  
> 0 0 stem P[0][0]["stem"]  
> 0 (N-1) bulge P[0][N-1]["bulge"]  
> :  
> 0 (N-1) stem P[0][N-1]["stem"]  
> :  
> (R-1) (N-1) bulge P[R-1][N-1]["bulge"]  
> :  
> (R-1) (N-1) stem P[R-1][N-1]["stem"]  



