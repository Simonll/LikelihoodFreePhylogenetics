# LikelihoodFreePhylogenetics


./LFP -m CodonMutSelFinite ANKRD28-M0GTRW-A-M0GTRWCpG.conf

\###### Configuration file ANKRD28-M0GTRW-A-M0GTRWCpG.conf ######

\#SUMMARIES  A C D E F G H I K L M N P Q R S T V W Y <br>
\#PARAM  chainID root    lambda_CpG  lambda_TBL  lambda_omega    nucsA   nucsC   nucsG   nucsT   nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT <br>
\#MAP    gtnrAC  gtnrAG  gtnrAT  gtnrCA  gtnrCG  gtnrCT  gtnrGA  gtnrGC  gtnrGT  gtnrTA  gtnrTC  gtnrTG  dinuc12CGTG dinuc12CGCA dinuc23CGTG dinuc23CGCA dinuc31CGTG dinuc31CGCA Nsub  Nsynsub <br>
\#NRUN 10 1000 <br>
\#NTHREADS 1 <br>
\#OUTPUT /ppred/ANKRD28-M0GTRW-A-LFP-M0GTRWCpG <br>
\#LOCALPARAM -d ANKRD28.puz -chain ANKRD28-M0GTRW-A -iscodon -code Universal -start 600 -every 4 -until 1000 -fixroot Echinops Procavia 0.9 -lambdaCpG 1.0 -lambdaTBL 1.0 -lambdaOmega 1.0 -tophylip -rootlength 100 <br>

Where lambdaCpG is a multiplicative parameter that modulate the mutation rate of CpG context (i.e., C->T, G->A) <br>
Where lambdaTBL is a multiplicative parameter that modulate the branches lengths <br>
Where lambdaOmega is a multiplicative parameter that modulate Omega or Omega* (i.e., Omega X lambdaOmega or Omega* X lambdaOmega) <br>






