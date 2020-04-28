require(plyr)
require(dplyr)
require(parallel)
require(sp)
require(rptR)
require(heatmap3)
require(ggplot2)
require(sparseLDA)

np=4 # parallel threads

# raw data (merged from Choreography)
load('S12_merged_roller_MWT.rda')
# per track summary stats
load('S13_merged_roller_MWT_trackStat.rda', verbose=T)
colo = merge(colo, muts)

# add mutant genotypes, excluding NA
cresv = na.exclude(cres)

# convex hull in 30s bins per track
hulld <- function(res, nlong=100, minPersistence=30, discard=4, bin=NULL) {
  
  # for longest nlong tracks (in time), 
  # get length traveled, area of path (convex hull), adjusted track length
  # after discarding first discard minutes
  # if bin, split tracks by bin seconds
  
  res$id = as.character(res$id)
  print(res$id[1])
  res <- res[,c('id', 'block', 'Time', 'IndWorm', 'XPosition', 'Yposition', 'Bias', 'Persistence')]
  names(res)[5:6] = c('x', 'y')
  # discard early data
  df <- subset(res, Time > (discard*60) & Persistence > minPersistence)
  # recalc track lengths
  trackl = aggregate(data = df, Time ~ IndWorm+block, function(x) diff(range(x)))
  trackl = subset(trackl, Time > minPersistence)
  df = subset(df, IndWorm %in% trackl$IndWorm)
  df <- merge(df, trackl[,c('IndWorm', 'block')])
  if(nrow(trackl) > nlong){
    # take the longest (time)
    trackl <- trackl[order(trackl$Time, decreasing = T),]  
    df <- merge(df, trackl[1:nlong,c('IndWorm', 'block')])
  }
  df <- df[order(df$IndWorm, df$Time),]
  df$bin = 1
  if(!is.null(bin)) df$bin = round(df$Time/bin)
  
  df <- ldply(split(df, list(df$IndWorm, df$bin)), .id=NULL, function(x) if(nrow(x)>3) cbind(x, disp = c(0, sapply(2:nrow(x), function(i) dist(x[(i-1):i,c('x', 'y')])))))
  df$block = as.character(df$block)
  maxd = ddply(df, .(block, bin), function(b) {xd = aggregate(data=b, disp~IndWorm+block, sum); xd[order(xd$disp, decreasing = T),][1:nlong,]})
  ldf <- merge(df, maxd[,c('IndWorm', 'block', 'bin')])
  
  # ratio of distance travelled to area of convex hull, matched for time
  yars = ldply(split(ldf, list(ldf$IndWorm, ldf$block, ldf$bin)), .id=NULL, function(y) {
    if(nrow(y)>4) cbind(y[1,c('block', 'IndWorm', 'bin')], peri = perimeter(y$x, y$y), hulla = Polygon(as.matrix(y[chull(y[,c('x','y')]),c('x','y')]))@area)
  })
  
  binl = aggregate(data = df, Time ~ IndWorm+block+bin, function(x) diff(range(x)))
  
  yars = yars[order(yars$IndWorm, yars$bin, yars$block),]
  yars = merge(yars, maxd, sort=F)
  yars = merge(yars, binl)
  rownames(yars) <- NULL
  cbind(id=res$id[1], yars)
}

hullas = bind_rows(mclapply(split(cresv, list(cresv$id, cresv$block)), mc.cores=np, function(x) hulld(x,bin=30)))
hullas = merge(muts, hullas)
hullas$circ = 4*pi*hullas$hulla/hullas$peri^2
# per worm means from 30s bins
hullm = hullas %>% group_by(line, block, col.geno, int, IndWorm) %>% summarise_at(.vars = c('hulla', 'disp', 'Time', 'circ'), mean)
# plate circularity means
circm = hullm %>% group_by(line, col.geno, int, block) %>% summarise_at(.vars = 'circ', function(x) mean(log(x), na.rm=T))

tracks = subset(colo, trackStartTime > 240 & persistence>30); dim(tracks); table(tracks$id)
tracks <- colo[,c('line', 'id', 'col.geno', 'int', 'block', 'date', 'xtime', 'worm', ctraits)]
tracks = na.exclude(tracks)

# log TF where appropriate (> minImp in -log10 units)
minimp=6
for(i in ctraits) {
  x = tracks[,i]
  if(min(x, na.rm=T) == 0) {
    x <- x + min(abs(x)[abs(x)>0], na.rm=T)
  } else {
    if(min(x, na.rm=T) < 0) x <- (x-min(x, na.rm=T)) + abs(min(x, na.rm=T))
  }
  # sample a few times since n is capped
  a = median(unlist(lapply(1:10, function(j) -log10(shapiro.test(sample(x, 5000))$p))))
  al = median(unlist(lapply(1:10, function(j) -log10(shapiro.test(log(sample(x, 5000)))$p))))
  if((a-al)>minimp){
    tracks[,i] = log(x)
    names(tracks)[names(tracks)==i] <- paste0('ln.', i)
  }
}
ctraits = names(tracks[-(1:8)])

# take out block effects for all MWT traits
tracks$xtime = as.numeric(as.character(tracks$xtime))
tres <- data.frame(do.call(cbind, lapply(ctraits, function(i) scale(resid(lm(as.formula(sprintf('%s ~ block', i)), tracks))))))
names(tres) = ctraits
tres <- cbind(tracks[,1:8], tres)

# get repeatabilities
reps = unlist(mclapply(ctraits, mc.cores = np, function(y) {
  df = cbind(data.frame(y = tres[,y]), tres[,c('col.geno', 'int', 'block', 'id')])
  vv <- aggregate(data = df, y ~ id+block, mean)
  rpt(y~(1|id), grname = 'id', data = vv, datatype = 'Gaussian', nboot=0)$R
}))
names(reps) = ctraits
# add circularity
circr = rpt(circ~(1|line), grname = 'line', data = circm, datatype = 'Gaussian', nboot=0)$R; names(circr) = 'circ'
reps = c(reps, unlist(circr))
(reps = reps[order(reps, decreasing = T)])

selectTraitsBySDA <- function(reps, K=5, maxr = 0.5, minrep = 0.5){

  # prune on repeatability, trait correlations  
  tix = names(reps[reps>minrep])
  tix = tix[tix!='circ']
  ccor = cor(tres[,tix])^2
  print(heatmap3(ccor, symm = T))
  diag(ccor) = 0
  dropt = NULL; for(i in 1:length(tix)) if(!tix[i] %in% dropt) dropt = unique(c(dropt, tix[which(ccor[,i]>maxr)]))
  (tix = tix[!tix %in% dropt])
  ccor = cor(tres[,tix])^2
  print(heatmap3(ccor, symm = T))
  
  # select based on sparse LDA, differentiating rol-1~ and sqt-3~ separately.
  # using plate means, get top K-1(=5) predictors.
  # add mean circularity
  platem <- tres %>% group_by(line, col.geno, int, block) %>% summarise_at(.vars = tix, mean)
  platem = merge(platem, circm)
  platem$dummy = platem$col.geno=='anc.'; platem$dummy[platem$line %in% c('N2', 'COP1834')] = T
  table(platem$line, platem$dummy)
  inix = c(tix, 'circ')
  # rol-1
  dd = platem[platem$int!='sqt-3(sc8)',c('dummy', 'line', inix)]; table(dd$line, dd$dummy)
  (dafit = sda(x=as.matrix(dd[,inix]), y=model.matrix(~-1+dummy, dd), stop=-K))
  (rot = dafit$varNames[order(abs(dafit$beta[,1]), decreasing = T)][1:K])
  rolbetas = matrix(dafit$beta[match(rot, dafit$varNames)]); rownames(rolbetas) = rot
  # plot the functions
  roro <- as.matrix(platem[,rot]) %*% rolbetas
  print(ggplot(cbind(platem, rot=roro), aes(line, rot)) + stat_summary() + facet_grid(.~int+col.geno, scales='free') + ggtitle('rol-1 D'))
  
  # sqt-3
  dd = platem[platem$int!='rol-1(sc22)',c('dummy', 'line', inix)]; table(dd$line, dd$dummy)
  (dafit = sda(x=as.matrix(dd[,inix]), y=model.matrix(~-1+dummy, dd), stop=-K))
  (sqt = dafit$varNames[order(abs(dafit$beta[,1]), decreasing = T)][1:K])
  sqtbetas = matrix(dafit$beta[match(sqt, dafit$varNames)]); rownames(sqtbetas) = sqt
  # plot the functions
  sqro <- as.matrix(platem[,sqt]) %*% sqtbetas
  print(ggplot(cbind(platem, rot=sqro), aes(line, rot)) + stat_summary() + facet_grid(.~int+col.geno, scales='free') + ggtitle('sqt-3 D'))
  
  # return loadings
  list(rolbetas, sqtbetas)
}

# use all (thinned) traits for Figure 3
o = selectTraitsBySDA(reps)
(rol=o[[1]])
(sqt=o[[2]])

# -kink for CeMEE mapping functions (wasn't measured for RIL locomotion data)
o = selectTraitsBySDA(reps[grep('kink', names(reps), invert=T)])
(rolk=o[[1]])
(sqtk=o[[2]])

(tix = unique(c(rownames(rol), rownames(sqt))))
save(rol, sqt, rolk, sqtk, tix, file = 'RolSqt_SDA_betas.rda')

univariate_plots <- function(traitx, tres, circm){
  
  pth = theme_classic() + theme(strip.text = element_text(face='italic', size=6), strip.background = element_rect(size=0), axis.ticks.x = element_blank(), panel.spacing = unit(1, 'mm'), axis.text = element_text(size=8), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  
  uplot <- function(trait, pout, splot=T){
    
    # plate means from all tracks (block residuals)
    df = cbind(data.frame(y = tres[,trait]), tres[,c('col.geno', 'int', 'block', 'line', 'id', 'worm')])
    dfm = aggregate(data = tres, as.formula(sprintf('%s~line+block+col.geno+int+id', trait)), mean); names(dfm)[6] = 'y'
    dfm$line = factor(dfm$line)
    dfm$line = factor(dfm$line, levels = c('N2','COP1834','BE22','QG2957','BE8','QG3070'))
    
    p <- ggplot(df, aes(line, y)) + 
      geom_boxplot(data = dfm, size=1, outlier.shape = NA) + 
      stat_summary(aes(group = block), fun.y = mean, geom='point', alpha=0.5) + 
      labs(x='', y = pout) + pth +
      geom_text(data = data.frame(x = 1:6, y = rep(max(dfm$y) + diff(range(dfm$y))*.2, 6), l=rep(c('-', '+'), 3)), aes(x, y, label=l), size=5) +
      geom_text(data = data.frame(x = c(3.5, 5.5), y = rep(max(dfm$y) + diff(range(dfm$y))*.3, 2), l=c('rol-1', 'sqt-3')), aes(x, y, label=l), size=3, fontface='italic') + 
      coord_cartesian(ylim = c(min(dfm$y), max(dfm$y) + diff(range(dfm$y))*.3)) +
      geom_tile(data = data.frame(x = 3.5, y = 0), aes(x, y, width=2, height=diff(range(dfm$y)*3)), alpha=0.1)
    
    if(splot) ggsave(sprintf('mwt_%s.png', trait), h = 4, w = 1.5)
    p
  }
  
  ptraits = traitx[traitx!='circ']
  plabs = capwords(gsub('.F', ' (F)', gsub('.S', ' (S)', gsub('.B', ' (B)', gsub('.var', ' variance', gsub('ln.', '', ptraits))))))
  lapply(1:len(ptraits), function(i) uplot(ptraits[i], plabs[i]))
  
  # and circularity
  circm$line = factor(circm$line)
  circm$line = factor(circm$line, levels=levels(circm$line)[c(5,4,1,6,2,7)])
  p <- ggplot(circm, aes(line, circ)) + 
    geom_boxplot(data = platem, size=1, outlier.shape = NA) + 
    stat_summary(aes(group = block), fun.y = mean, geom='point', alpha=0.5) + 
    labs(x='', y = 'Track circularity') + pth + 
    geom_text(data = data.frame(x = 1:6, y = rep(max(platem$circ) + diff(range(platem$circ))*.2, 6), l=rep(c('-', '+'), 3)), aes(x, y, label=l), size=5) +
    geom_text(data = data.frame(x = c(3.5, 5.5), y = rep(max(platem$circ) + diff(range(platem$circ))*.3, 2), l=c('rol-1', 'sqt-3')), aes(x, y, label=l), size=3, fontface='italic') +
    coord_cartesian(ylim = c(min(platem$circ), max(platem$circ) + diff(range(platem$circ))*.3)) +
    geom_tile(data = data.frame(x = 3.5, y = 0), aes(x, y, width=2), alpha=0.1)
  ggsave(sprintf('mwt_%s.png', 'circ'), h = 4, w = 1.5)
  
}

# plot plate means for traitx
univariate_plots(tix, tres, circm)

multivariate_plots <- function(traitx, tres, circm, opref){
  
  traitx = traitx[traitx!='circ']
  platem <- tres %>% group_by(line, col.geno, int, block) %>% summarise_at(.vars = traitx, function(x) mean(x, na.rm=T))
  if('circ' %in% traitx){
    platem = merge(platem, circm)
    traitx = c(traitx, 'circ')  
  }
  
  md = cmdscale(dist(scale(platem[,traitx], center = T, scale = F)))
  md = cbind(data.frame(md), data.frame(platem[,1:4]))
  md$col.geno = factor(md$col.geno, labels = c('ancestral', 'N2'))
  md$int = as.factor(as.character(md$int))
  levels(md$int)[1] = 'None'
  md$int = factor(md$int, labels = sprintf('italic("%s")', levels(md$int)))
  p <- ggplot(md, aes(X1, X2, shape=col.geno, col=int)) + geom_point(size=3, stroke=0, alpha=0.75) + theme_classic() + 
    stat_ellipse(aes(linetype=col.geno), level=0.8, size=1.2) + scale_color_brewer('Interaction', type='qual', palette = 2, labels = scales::parse_format()) +
    labs(x = 'MDS coordinate 1', y = 'MDS coordinate 2', shape = expression(italic('col-182')), linetype = expression(italic('col-182'))) +
    guides(col = guide_legend(override.aes = list(shape=15, linetype=NA, size=4)), shape = guide_legend(override.aes = list(alpha=0.5))) +
    theme(legend.position = 'top', legend.margin = margin(0,0,0,0), plot.margin = margin(c(0,1,0,1)))
  ggsave(sprintf('mwt_mds_%s.png', opref), h = 2.5, w = 6)
  
}

multivariate_plots(c(ctraits, 'circ'), tres, circm, opref='rep')
multivariate_plots(tix, tres, circm, opref='sda')
multivariate_plots(rownames(rol)[1:4], tres, circm, opref='rolsda')
multivariate_plots(rownames(sqt)[1:4], tres, circm, opref='sqtsda')



