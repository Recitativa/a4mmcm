Simulation
simqq:	  simulations for the discretization error of quantiles. 
	  Relevent source: 
	  theoretic.hpp	   theoretic.cpp	simquantile.cpp
	  
	  Usage:
	  simqq Terms Rb Rm Re sigma mu OutputIndex
	  Terms:      number of simulation
	  Rb:	      smallest number of points is:2^Rb, e.g. Rb=4 in my paper
	  Rm:	      lagest number of points for computing quantiles is 2^Rm,
	  	      Rm=16 in my paper.  
	  Re:	      Finest discretization is 2^Re points, e.g. Re=25 
	  sigma:      volatility
	  mu:	      drift
	  OutputIndex:	number will be append at the end of outputfile, e.g 10
	  
	  Remark: Quantile will be computed at time T=1. 
	  
	  Result:
	  Messages will be print to the conslon, to indicate how many times
	  simulation has be down. 
	  An binary file nout_Rb_Rm_Re_sigma_mu_OutputIndex.bin, 
	  e.g. nout_4_16_25_1_0_10.bin
	  
treat_m.r: Purpose:
	   To treat the binary file to get a plot, use this R program.
	   
	   Usage:	   
	   Suppose we have several binary files with index from 100:114,
	   >>> ns <- c(paste("nout_4_16_25_1_3_",100:114,".bin",sep=""))
	   >>> source("treat_m.r"); rrr <- run(ns);
	   the return "rrr" is a list containing the log of errors, 
	   theoretical and simulation results of the intersections.  
	   relevent PDF files will be created.	   
	   
	   Purpose:
	   To print the Figure of coefficients for different mu.
	   use:
	   >>> source("treat_m.r");pdf("mu.pdf");plotI();dev.off()

Quantile options:
bfsg:	 compute the price of quantile call option by FSG.
	 bfsg S K alpha r sigma T Bn En
	 S:   initial stock price
	 K:   strike price
	 Bn:  smallest number of the steps e.g.  10,
	 Be:  largest number of the steps e.g.  100.
	 A text file will be created to store the result. 
biam:	 compute the price of American quantile call option by binomial tree.  
	 same parameter and output pattern as above. 
biit:	 compute the price of European quantile call option
	 same parameter and output pattern as above. 
	 
simEro.r: simulation methord to compute the quantile option.
	  Useage:
	  >>> S=100;K=100;r=0.05;sigma=0.2;T=1;alpha=0.5;
	  >>> source("brod.r");
	  >>> rr<-price(S,K, r, sigma, T, alpha, n=10000000);

treat_bitree.r: extrapolation relevant programs.

Compile the c++ programs:

All programs are compiled in CYGWIN, using OMP and GSL libaries. 
You can compile all the programs simply by using makefile.
e.g.:
c:\prog> make

Or compile them on any linux system, for example, the HPC in NUS. 
