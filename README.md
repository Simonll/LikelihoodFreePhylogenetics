# LikelihoodFreePhylogenetics


./HS -m CodonMutSelFinite ANKRD28-M0GTRW-A-M0GTRWCpG.conf

\###### Configuration file ANKRD28-M0GTRW-A-M0GTRWCpG.conf ######

\#SUMMARIES  A C D E F G H I K L M N P Q R S T V W Y <br>
\#PARAM  chainID root    lambda_CpG  lambda_TBL  lambda_omega    nucsA   nucsC   nucsG   nucsT   nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT <br>
\#MAP    gtnrAC  gtnrAG  gtnrAT  gtnrCA  gtnrCG  gtnrCT  gtnrGA  gtnrGC  gtnrGT  gtnrTA  gtnrTC  gtnrTG  dinuc12CGTG dinuc12CGCA dinuc23CGTG dinuc23CGCA dinuc31CGTG dinuc31CGCA Nsub  Nsynsub <br>
\#GENES ANKRD28.puz <br>
\#CHAINS ANKRD28-M0GTRW-A <br>
\#NRUN 10 1000 <br>
\#NTHREADS 1 <br>
\#OUTPUT /ppred/ANKRD28-M0GTRW-A-LFP-M0GTRWCpG <br>
\#OLDPARAMS -d ANKRD28.puz -chain ANKRD28-M0GTRW-A -iscodon -code Universal -start 600 -every 4 -until 1000 -fixroot Echinops Procavia 0.9 -lambdaCpG 1.0 -lambdaTBL 1.0 -lambdaOmega 1.0 -tophylip<br>

Here 10 replicate per MCMC pt will be generated for 100 MCMC pt, leading to 1000 simulations <br>
log2 amino acid usage will be save in a ANKRD28-M0GTRW-A-LFP-M0GTRWCpG.post file as well as squared discrepancies (e.g., D_A)  and the sum of squared discrepancies (i.e., D_sum) <br>
Statistics on sampled evolutionary history will be save in that same file, where gtnrAC is the number of substitutions from A->C occuring along the tree for a given simulation, and were dinuc12CGTG is the number of substitutions C->T occuring in CpG context within codon position 12. Finaly, the number of synonymous substitution is recorded (i.e., Nsynsub) as well as the total number of substitutions (i.e., Nsub). <br> 


Where lambdaCpG is the mutation rate of transitions related to CpG context (i.e., C->T, G->A) <br>
Where lambdaTBL modulate treelength or global mutation rate (i.e., TotalTreeLength X lambdaTBL or mutation rate X lambdaTBL) <br>
Where lambdaOmega modulate Omega and Omega* (i.e., Omega X lambdaOmega or Omega* X lambdaOmega) <br>






