# LikelihoodFreePhylogenetics<br>
## Options <br>
./LFP <br>
version 1.0 <br>
\########################### <br>
-m < stats | show | CodonMutSelFiniteABC | CodonMutSelSBDPABC | CodonDegMutSelFiniteABC | CodonDegMutSelSBDPABC | CodonMutSelFinite | CodonMutSelSBDP | CodonMutSelFinitePPred | CodonMutSelSBDPPPred > <controlfile> <br>
\###########################<br>
\#SUMMARIES<br>
\#ANCSUMMARIES<br>
\#ACCSUMMARIES<br>
\#PARAM<br>
\#SSMAP<br>
\#MAP<br>
\#ANCESTRALMAP<br>
\#DIST<br>
\#TRANS<br>
\#SPEUDODATA<br>
\#NRUN<br>
\#SAMPLING<br>
\#LOCALPARAM<br>
\###########################<br>


## Running rejection sampling<br>
./LFP -m CodonMutSelSBDPABC ACO1-2673-MSAAGTRW-A-GTRWCpG.conf<br>

### Configuration file:<br>
\#SUMMARIES	pwAC	pwAG	pwAT	pwCG	pwCT	pwGT	dinuc31CG	dinuc31TG	dinuc31CA	nuc3A	nuc3C	nuc3G	nuc3T	pwAA<br>
\#PARAM	chainID	root	lambda_CpG	lambda_TBL	lambda_omega	nucsA	nucsC	nucsG	nucsT	nucrrAC	nucrrAG	nucrrAT	nucrrCG	nucrrCT	nucrrGT<br>
\#MAP	gtnrAC	gtnrAG	gtnrAT	gtnrCA	gtnrCG	gtnrCT	gtnrGA	gtnrGC	gtnrGT	gtnrTA	gtnrTC	gtnrTG	dinuc31CGTG	dinuc31CGCA	Nsub	Nsynsub<br>
\#SAMPLING 600 4 999<br>
\#NRUN 1000000 10000<br>
\#NTHREADS 6<br>
\#OUTPUT /ABC/ACO1-2673-MSAAGTRW-A-GTRWCpG<br>
\#LOCALPARAM -d ACO1.puz -chain /step_1/ACO1-2673-MSAAGTRW-A -code Universal -iscodon -freeroot Echinops Procavia -freelambdaTBL -freelambdaCpG -freegtr -freelambdaomega -rootlength 10 -priorlambdaomega log2Unif -priorlambdaCpG log10Unif -priorlambdaTBL log2Unif<br>

## Simulating from M[GTR+λCpG]-[1CatAA] or M[GTR+λCpG]-[1CatCodon]<br>

./LFP -m CodonMutSelFinite ANKRD28-M0GTRW-A-M0GTRWCpG.conf<br>

### Configuration file: ANKRD28-M0GTRW-A-M0GTRWCpG.conf<br>

\#SUMMARIES  A C D E F G H I K L M N P Q R S T V W Y <br>
\#PARAM  chainID root    lambda_CpG  lambda_TBL  lambda_omega    nucsA   nucsC   nucsG   nucsT   nucrrAC nucrrAG nucrrAT nucrrCG nucrrCT nucrrGT <br>
\#MAP    gtnrAC  gtnrAG  gtnrAT  gtnrCA  gtnrCG  gtnrCT  gtnrGA  gtnrGC  gtnrGT  gtnrTA  gtnrTC  gtnrTG  dinuc12CGTG dinuc12CGCA dinuc23CGTG dinuc23CGCA dinuc31CGTG dinuc31CGCA Nsub  Nsynsub <br>
\#NRUN 10 1000 <br>
\#NTHREADS 1 <br>
\#OUTPUT /ppred/ANKRD28-M0GTRW-A-LFP-M0GTRWCpG <br>
\#LOCALPARAM -d ANKRD28.puz -chain ANKRD28-M0GTRW-A -iscodon -code Universal -start 600 -every 4 -until 1000 -fixroot Echinops Procavia 0.9 -lambdaCpG 1.0 -lambdaTBL 1.0 -lambdaOmega 1.0 -tophylip -rootlength 10 <br>

Where  λCpG is a multiplicative parameter that modulate the transition rate of CpG context<br>
Where  λTBL is a multiplicative parameter that modulate the all the branch lengths<br>
Where  λomega \  λomega* is a multiplicative parameter that modulate Omega or Omega*<br>
