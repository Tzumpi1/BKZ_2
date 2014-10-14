BKZ_2
=====

This repo adds on top of NTL (http://www.shoup.net/ntl/doc/tour-win.html) by updating the lattice reduction technique known as Base Korkin Zolotarev (BKZ), originated by C. Schnorr in the article - "
Lattice Basis Reduction: Improved Practical Algorithms and Solving Subset Sum Problems. (1993)".

This implementation (BKZ v2) follows closely the article "BKZ 2.0: Better Lattice Security Estimates" (2012) by  Yuanmi Chen and Phong Q. Nguyen.

We therefore updated the code of NTL (namely one single file: LLL_FP) with a new method: 

long BKZ_FP_v2(mat_ZZ& BB, long alfa,  long beta, double gamma, double delta, 
	       double prunningProbability, long vNumRetry, BKZCheckFct check);

where:

- BB is the lattice to reduce
- alfa is the preprocessing block size (see Phong's article)
- beta is the block size (same as the original BKZ)
- gamma is the Gaussian Heuristic enumeration radius (see Phong's article)
- delta is for LLL (see LLL reduced basis)
- prunning probability is between 0.05 and 0.95 and establishes a curve that we computed according to Phong's article (for now - all curves are hard-coded and the probability only makes the selection) - see Phong's article
- vNumRetry is by default 1 and is not usable at this stage. Its purpose would be to give a number of retries when reducing one block (given that we randomize the GS norms...). At this stage we did not manage to randomize the GS norms of a local block !!!
- check is a callback function that can be used to inspect the current progress of the reduction.



Everyone is encouraged to freely use and test this... 

To use this, take the LLL_FP file, overwrite it on top of your NTL's sources, then recompile !
No changes have been done in the original NTL functions !!!


For any questions, pls send an email to dan underscore durbaca at yahoo dot com.


