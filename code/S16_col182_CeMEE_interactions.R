#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = 1)
args  = commandArgs(trailingOnly = T)
NP    = as.numeric(args[1]) # parallel threads
job   = as.numeric(args[2]) # 1-3: curvature/circularity, rol/sqt SDA, length/width
nperm = 50                  # dump a file every nperm permutations

require(data.table())
require(dplyr)
require(plyr)
require(parallel)
require(car)

load('S7_CeMEE_phenoGeno.rda', verbose=T)
# gtf  : SNP genotypes
# cold : col-182 indel genotypes in CeMEE RILs
# loco : CeMEE Multi-Worm Tracker-derived least squares plate mean phenotypes
# exm  : CeMEE exploratory behaviour line means (using log circularity only)
# sdas : CeMEE rotated traits for rol-1/sqt-3 discriminant functions
# utility functions `getBLUPs`, `filterMAF`

# do the tests, save results and permutations
runmvt <- function(opref, traits, blups, nperm, NP=1, ix=1){
  
  MVGInt <- function(gtf, colhaps, blups, traits, minclass=10, perm=F, np=4){
    
    popix = grep('popr', names(blups))
    if(perm) blups = ldply(split(blups, blups[,popix]), function(x) {x[,traits] = x[sample(nrow(x)),traits]; x})
    
    iline = names(gtf)[names(gtf) %in% blups$line]
    gtf <- filterMAF(gtf, MAFgt = 0.05, lines = iline)
    X = as.matrix(gtf[,names(gtf) %in% iline])
    colhaps = subset(colhaps, line %in% blups$line)
    snps = gtf[,c('chrom', 'pos')]
    ix = 1:nrow(snps)
    yst = paste(traits, collapse=',')
    
    colint = bind_rows(mclapply(ix, mc.cores = np, function(i) {
      dx = cbind(colhaps, x = X[i,])
      dx = subset(dx, x%%1==0)
      if(nrow(dx)>0){
        ht = table(dx$x, dx$col182)
        if(min(ht)>=minclass){
          dx = merge(dx, blups)
          m00 = lm(as.formula(sprintf('cbind(%s)~1', yst)), dx)
          m0 = lm(as.formula(sprintf('cbind(%s)~x+col182', yst)), dx)
          m1 = lm(as.formula(sprintf('cbind(%s)~x*col182', yst)), dx)
          cbind(snps[i,c('chrom', 'pos')], 
                full = -log10(anova(m1, m00)[2,8]), 
                int = -log10(anova(m0, m1)[2,8])
          )
        }  
      }
    }))
    colint[order(colint$chrom, colint$pos),]
  }
  
  blups = subset(blups, line %in% cold$line)
  
  handle = sprintf('%s_N2_MVindelGenInt', opref)
  resf = sprintf('%s_res.rda', handle)
  if(!file.exists(resf)) {
    colint = MVGInt(gtf, cold, blups, traits, np=NP)
    save(colint, file = resf)
  } else {
    load(resf, verbose=T)
  }
  geno = subset(gtf, paste(gtf$chrom, gtf$pos) %in% paste(colint$chrom, colint$pos))
  
  perms = bind_rows(lapply(1:nperm, function(i) {
    o = MVGInt(geno, cold, blups, traits, perm=T, np=NP)
    o = melt(o, 1:2)
    cbind(merge(o, aggregate(data=o, value ~ variable, max)), perm=i)
  }))
  
  permf = sprintf('%s_%s.rda', handle, ix)
  while(file.exists(permf)){
    ix = ix+1
    permf = sprintf('%s_%s.rda', handle, ix)
  }
  save(perms, blups, traits, file = permf)
}

# prepare the data for a single trait pair
if(job==1){
  
  traits = c('curvature', 'ln.circ')
  loco$curvature = log(loco$curvature)
  blups1 = lapply(traits[1], function(x) getBLUPs(loco, x, popres=T))
  blups = blups1 %>% Reduce(function(d1,d2) full_join(d1,d2, by=c('line', 'popr')), .)
  blups = merge(blups, getBLUPs(subset(exm, env=='NGM' & line %in% names(gtf)), traits[2], popres=T))
  maha <- mahalanobis(blups[,traits], center = apply(blups[,traits], 2, mean), cov = cov(blups[,traits]))
  outl = blups$line[which(maha>quantile(maha, 0.995, na.rm=T))]
  blups <- subset(blups, line %in% names(gtf) & !(line %in% outl))
  while(T) runmvt('cuci', traits, blups, nperm, NP=NP)
}

if(job==2){
  traits = c('rolsd', 'sqtsd')
  blups = sdas
  # mean center pops
  blups = bind_rows(lapply(split(blups, blups$popr), function(x) {x[,traits] = scale(x[,traits], scale = F); x}))
  maha <- mahalanobis(blups[,traits], center = apply(blups[,traits], 2, mean), cov = cov(blups[,traits]))
  outl = blups$line[which(maha>quantile(maha, 0.995, na.rm=T))]
  blups <- subset(blups, line %in% names(gtf) & !(line %in% outl))
  while(T) runmvt('sda', traits, blups, nperm, NP=NP)
}

if(job==3){
  traits = c('length.F', 'width.F')
  blups1 = lapply(traits, function(x) getBLUPs(subset(loco, env=='NGM' & line %in% names(gtf)), x, popres=T))
  blups = blups1 %>% Reduce(function(d1,d2) full_join(d1,d2, by=c('line', 'popr')), .)
  maha <- mahalanobis(blups[,traits], center = apply(blups[,traits], 2, mean), cov = cov(blups[,traits]))
  outl = blups$line[which(maha>quantile(maha, 0.995, na.rm=T))]
  blups <- subset(blups, line %in% names(gtf) & !(line %in% outl))
  while(T) runmvt('lw', traits, blups, nperm, NP=NP)
}

# test interactions against an additive model by parametric bootstrap
# (for genome-wide significant QTL)
hier2dMV <- function(dfi, traits, np=4, nboot=1000, returnBoots=F){
  
  # multivariate bootstrap for pairwise interaction against an additive model
  # see https://www.ncbi.nlm.nih.gov/pubmed/20384625
  # data.frame dfi has `x1`, `x2` numeric predictors to test for interaction
  # and responses in `traits` vector
  # if returnBoots, return [[summary stats]], [[null stats]]
  stopifnot(sum(names(dfi) %in% c('x1', 'x2', traits))==(length(traits)+2))
  require(parallel)
  
  boots <- function(dfj, np, n) {
    unlist(mclapply(1:n, mc.cores=np, function(x){
      # joint resampling with replacement
      ys = scale(m1$fitted.values)
      dfj <- cbind(dfj[,c('x1', 'x2')], ys[sample(nrow(ys), replace=T),])
      m01 = lm(sprintf('cbind(%s)~x1+x2', yst), dfj)
      m02 = lm(sprintf('cbind(%s)~x1*x2', yst), dfj)
      anova(m01, m02)[2,8]
    }))
  }
  
  dfi$x1 = scale(as.numeric(dfi$x1))
  dfi$x2 = scale(as.numeric(dfi$x2))
  yst = paste(traits, collapse=',')
  m0 = lm(sprintf('cbind(%s)~1', yst), dfi)
  m1 = lm(sprintf('cbind(%s)~x1+x2', yst), dfi)
  m2 = lm(sprintf('cbind(%s)~x1*x2', yst), dfi)
  
  # Pillai trace pvals
  t1p = anova(m0, m2)[2,8]
  t2p = anova(m1, m2)[2,8]
  bs = boots(dfi, np, nboot)
  o = data.frame(t1p, t2p, bs = sum(bs < t2p, na.rm=T)/nboot)
  if(returnBoots) return(list(o, bs))
  return(o)
}



