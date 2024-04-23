# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
# important functions are: all_measures, from_measures_to_table, from_any_table_to_rmarkdown, from_measure_to_calibration_plot, from_measure_to_ROC_curve, from_measure_to_PR_curve, table_one_way_with_percentage, wcmc.

## Calculate true precision for a threshold under binormal distribution.
binormal.precision <- function(c,pi=0.5,pos.mean=1.0,pos.sd=1.0)
{
  r = (pi*pnorm(c,pos.mean,pos.sd,lower.tail=FALSE))/(pi*pnorm(c,pos.mean,pos.sd,lower.tail=FALSE) + (1-pi)*pnorm(c,mean=0,sd=1,lower.tail=FALSE))
  r[is.nan(r)] = 1 ## nan from a denominator of 0 should be precision 1
  return(r)
}


## Calculate true recall for a threshold under binormal distribution.
binormal.recall <- function(c,pi=0.5,pos.mean=1.0,pos.sd=1.0)
{
  pnorm(c,pos.mean,pos.sd,lower.tail=FALSE)
}

## Calculate true precision for a particular recall under binormal distribution.
binormal.precision.at.recall <- function(recall, pi=0.5,pos.mean=1.0,pos.sd=1.0) {
  return(binormal.precision(qnorm(recall,pos.mean,pos.sd,lower.tail=FALSE),pi=pi,pos.mean=pos.mean,pos.sd=pos.sd))
}

# Plot the true binormal PR curve.
binormal.plot <- function(pi=0.5,pos.mean=1.0,pos.sd=1.0)
{
  plot(function(x) {binormal.precision(qnorm(x,pos.mean,pos.sd,lower.tail=FALSE),pi=pi,pos.mean=pos.mean,pos.sd=pos.sd)},xlim=c(0,1),ylim=c(0,1),xlab="Recall",ylab="Precision",sub=paste("pi=",pi,", pos.mean=",pos.mean,", pos.sd=",pos.sd))
}


# Calculate true area under binormal PR curve using numeric integration.
binormal.area <- function(pi=0.5,pos.mean=1.0,pos.sd=1.0)
{
  f <- function(q) { binormal.precision(qnorm(q,pos.mean,pos.sd,lower.tail=FALSE),pi=pi,pos.mean=pos.mean,pos.sd=pos.sd) }
  ## monte carlo integration
  #mean(f(runif(samples)))

  ## integrate method of R, much faster (although might have problems with x=0)
  r = integrate(f,0,1)
  return(r$value)
}


## Generate sample from binormal distribution. When pi.exact is false,
## the number of positive examples is distributed according to
## Binomial(n,pi). If pi.exact is false, the number of positive
## examples is always pi*n.
binormal.sample <- function(n=1, pi=0.5, pos.mean=1.0, pos.sd=1.0, pi.exact = TRUE) {
  if (pi.exact) {
    num.pos = as.integer(pi*n)
  }
  else {
    num.pos = rbinom(n=1,prob=pi,size=n)
  }
  num.neg = n-num.pos

  pos.values = rnorm(num.pos,mean=pos.mean,sd=pos.sd)
  neg.values = rnorm(num.neg,mean=0,sd=1)
  return(list(pos.values=pos.values,neg.values=neg.values))
}




######################################################################
### Bibeta
### Assumes scores are from two Beta distributions.
### Y ~ Beta(pos.a,pos.b)
### X ~ Beta(neg.a,neg.b)


## Calculates true precision for a threshold under bibeta distribution.
bibeta.precision <- function(c,pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  ## X ~ Beta(neg.a,neg.b), Y ~ Beta(pos.a,pos.b)
  r = (pi*pbeta(c,pos.a,pos.b,lower.tail=FALSE))/(pi*pbeta(c,pos.a,pos.b,lower.tail=FALSE) + (1-pi)*pbeta(c,neg.a,neg.b,lower.tail=FALSE))
  r[is.nan(r)] = 1 ## nan from a denominator of 0 should have precision of 1
  return(r)
}

## Calculates true recall for a threshold under bibeta distribution.
bibeta.recall <- function(c,pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  pbeta(c,pos.a,pos.b,lower.tail=FALSE)
}

## Calculate true precision for a particular recall under bibeta distribution.
bibeta.precision.at.recall <- function(recall,pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  return(bibeta.precision(qbeta(recall,pos.a,pos.b,lower.tail=FALSE),pi=pi,pos.a=pos.a,pos.b=pos.b,neg.a=neg.a,neg.b=neg.b))
}

## Plot the true bibeta PR curve.
bibeta.plot <- function(pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  plot(function(x) {bibeta.precision(qbeta(x,pos.a,pos.b,lower.tail=FALSE),pi=pi,pos.a=pos.a,pos.b=pos.b,neg.a=neg.a,neg.b=neg.b)},
       xlim=c(0,1),
       ylim=c(0,1),
       xlab="Recall",
       ylab="Precision",
       sub=paste("pi=",pi,",X (neg) ~ Beta(",neg.a,",",neg.b,"), Y (pos) ~ Beta(",pos.a,",",pos.b,")")
  )

}

## Calculate the true area under the bibeta PR curve using numeric integration.
bibeta.area <- function(pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1) {
  f <- function(q) { bibeta.precision(qbeta(q,pos.a,pos.b,lower.tail=FALSE),pi=pi,pos.a=pos.a,pos.b=pos.b,neg.a=neg.a,neg.b=neg.b) }

  ## use R's integrate method
  res = integrate(f,0,1)
  return(res$value)
}

## Generate sample from bibeta distribution. When pi.exact is false,
## the number of positive examples is distributed according to
## Binomial(n,pi). If pi.exact is false, the number of positive
## examples is always pi*n.
bibeta.sample <- function(n=1,pi=0.5,pos.a=2,pos.b=1,neg.a=1,neg.b=1,pi.exact = TRUE) {
  if (pi.exact) {
    num.pos = as.integer(pi*n)
  }
  else {
    num.pos = rbinom(n=1,prob=pi,size=n)
  }

  num.neg = n-num.pos

  ## Generates a sample of observed values for neg.values (X) ~
  ## Beta(neg.a,neg.b) and pos.values (Y) ~ Beta(pos.a,pos.b).

  pos.values = rbeta(num.pos,pos.a,pos.b)
  neg.values = rbeta(num.neg,neg.a,neg.b)

  return(list(pos.values=pos.values,neg.values=neg.values))
}


######################################################################
### Offset Uniform
### Assumes the scores are drawn from two uniform distributions
### spanning different ranges

### X ~ Uniform(0,1)
### Y ~ Uniform(a,a+b)


## Calculate true precision for a threshold under offset uniform distributions.
offsetuniform.precision <- function(c,pi=0.5,a=0.5,b=1) {
  r = (pi*punif(c,min=a,max=a+b,lower.tail=FALSE))/(pi*punif(c,min=a,max=a+b,lower.tail=FALSE) + (1-pi)*punif(c,min=0,max=1,lower.tail=FALSE))
  r[is.nan(r)] = 1 ## nan from denominator of 0 should have precision of 1
  return(r)
}

## Calculate true recall for a threshold under offset uniform distributions.
offsetuniform.recall <- function(c,pi=0.5,a=0.5,b=1) {
  punif(c,min=a,max=a+b,lower.tail=FALSE)
}

## Calculate true precision for a particular recall under offset uniform distribution.
offsetuniform.precision.at.recall <- function(recall, pi=0.5,a=0.5,b=1) {
  return(offsetuniform.precision(qunif(recall,a,a+b,lower.tail=FALSE),pi=pi,a=a,b=b))
}

## Plot the true offset uniform PR curve.
offsetuniform.plot <- function(pi=0.5,a=0.5,b=1) {
  plot(function(x) {offsetuniform.precision(qunif(x,min=a,max=a+b,lower.tail=FALSE),
                                            pi=pi,
                                            a=a,
                                            b=b)},
       xlim=c(0,1),
       ylim=c(0,1),
       xlab="Recall",
       ylab="Precision",
       sub=paste("pi=",pi," a=",a,"b=",b))
}

## Calculate true area under offset uniform PR curve using numeric integration.
offsetuniform.area <- function(pi=0.5,a=0.5,b=1) {
  f <- function(q) { offsetuniform.precision(qunif(q,min=a,max=a+b,lower.tail=FALSE),
                                             pi=pi,
                                             a=a,
                                             b=b)}
  r = integrate(f,0,1)
  return (r$value)
}

## Generate sample from offset uniform distribution. When pi.exact is false,
## the number of positive examples is distributed according to
## Binomial(n,pi). If pi.exact is false, the number of positive
## examples is always pi*n.
offsetuniform.sample <- function(n=1,pi=0.5,a=0.5,b=1,pi.exact = TRUE) {
  if (pi.exact) {
    num.pos = as.integer(pi*n)
  }
  else {
    num.pos = rbinom(n=1,prob=pi,size=n)
  }
  num.neg = n-num.pos

  ## X (neg) ~ Uniform(0,1)
  ## Y (pos) ~ Uniform(a,a+b)

  pos.values = runif(n=num.pos,min=a,max=a+b)
  neg.values = runif(n=num.neg,min=0,max=1)

  return(list(pos.values=pos.values,neg.values=neg.values))
}




aucpr.conf.int <- function(estimate,pos.values,neg.values,prcurve.estimator,conf.level=0.95,method="binomial",bootstrap.replicates=1000) {
  ## Calculates confidence interval for an AUCPR estimate.
  ## method=c("binomial","expit")
  m = match.arg(method,c("binomial","expit","bootstrap"))

  if (m=="binomial") {
    return(aucpr.conf.int.binomial(estimate,num.pos=length(pos.values),num.neg=length(neg.values),conf.level=conf.level))
  }
  else if (m=="expit") {
    return(aucpr.conf.int.expit(estimate,num.pos=length(pos.values),num.neg=length(neg.values),conf.level=conf.level))
  }
  else if (m=="bootstrap") {
    return(aucpr.conf.int.bootstrap(estimate,pos.values,neg.values,prcurve.estimator,conf.level,bootstrap.replicates))
  }
}

aucpr.conf.int.binomial <- function(estimate,num.pos,num.neg,conf.level=0.95) {
  ## Calculates confidence interval for an AUCPR estimate under
  ## binomial assumptions. Uses num.pos as the sample size (might not
  ## be best?).

  ci = estimate+qnorm(c((1-conf.level)/2,(1+conf.level)/2))*sqrt(estimate*(1-estimate)/num.pos)
  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "binomial"
  return(ci)
}

aucpr.conf.int.expit <- function(estimate,num.pos,num.neg,conf.level=0.95) {
  ## Calculates confidence interval for an AUCPR estimate using expit.

  ## convert to logit scale
  est.logit = log(estimate/(1-estimate))
  ## standard error (from Kevin Eng)
  se.logit = sqrt(estimate*(1-estimate)/num.pos)*(1/estimate + 1/(1-estimate))
  ## confidence interval in logit
  ci.logit = est.logit+qnorm(c((1-conf.level)/2,(1+conf.level)/2))*se.logit

  ## back to original scale
  ci = exp(ci.logit)/(1+exp(ci.logit))
  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "expit"
  return(ci)
}

aucpr.conf.int.bootstrap <- function(estimate,pos.values,neg.values,prcurve.estimator,conf.level=0.95,replicates=1000) {
  areas = rep(0,replicates)

  for (i in 1:replicates) {
    p.v = sample(pos.values,replace=TRUE)
    n.v = sample(neg.values,replace=TRUE)

    res = prcurve.estimator(p.v,n.v)
    areas[i] = res$area
  }

  q = quantile(areas,c((1-conf.level)/2,(1+conf.level)/2))

  ci = c(q[1],q[2])
  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "bootstrap"
  attr(ci,"median") = median(areas)
  attr(ci,"mean") = mean(areas)

  return (ci)
}

aucpr.conf.int.crossvalidation <- function(estimate,pos.values,neg.values,prcurve.estimator,conf.level=0.95,folds=10) {
  areas = rep(0,folds)

  pos = sample(pos.values)
  neg = sample(neg.values)

  pos.counts = rep(c(as.integer(length(pos.values)/folds)+1,as.integer(length(pos.values)/folds)),c(length(pos.values)%%folds,10-length(pos.values)%%folds))
  neg.counts = rep(c(as.integer(length(neg.values)/folds)+1,as.integer(length(neg.values)/folds)),c(length(neg.values)%%folds,10-length(neg.values)%%folds))

  pos.index = 1
  neg.index = 1
  for (k in 1:folds) {
    p = pos[pos.index:(pos.index+pos.counts[k]-1)]
    n = neg[neg.index:(neg.index+neg.counts[k]-1)]

    areas[k] = prcurve.estimator(p,n)$area

    pos.index = pos.index + pos.counts[k]
    neg.index = neg.index + neg.counts[k]
  }


  ## normal approximation
  ##ci = mean(areas) + qnorm(c((1-conf.level)/2,(1+conf.level)/2))*sd(areas)/sqrt(length(areas))
  ## use t-distribution
  ci = mean(areas) + qt(c((1-conf.level)/2,(1+conf.level)/2),df=length(areas)-1)*sd(areas)/sqrt(length(areas))

  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "crossvalidation, t-dist"
  attr(ci,"values") = areas
  attr(ci,"mean") = mean(areas)

  return(ci)
}





# return confusion matrix info for >= threshold
# assume labels are 0 and 1
confusion.matrix <- function(values,labels,threshold)
{
  counts=table(factor(values>=threshold,c(FALSE,TRUE)),labels)
  ## print(counts)

  TP = counts[2,2]
  TN = counts[1,1]
  FP = counts[2,1]
  FN = counts[1,2]

  if (TP+FP == 0)
    return(list(tp=TP,tn=TN,fp=FP,fn=FN,precision=NaN,sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),accuracy=(TP+TN)/(TP+TN+FP+TN)))
  else
    return(c(tp=TP,tn=TN,fp=FP,fn=FP,precision=TP/(TP+FP),sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),accuracy=(TP+TN)/(TP+TN+FP+TN)))
}


## Create all confusion matrices from set of pos.values and neg.values
make.confusion.matrices <- function(pos.values,neg.values) {
  ## thresholds, sorting to have increasing recall
  thresholds = sort(c(-Inf,unique(c(pos.values,neg.values)),Inf),decreasing=TRUE)

  ## 0 - negative (control)
  ## 1 - positive (case)
  labels = c(rep(1,length(pos.values)),rep(0,length(neg.values)))
  values = c(pos.values,neg.values)

  # create set of confusion matrices
  d=sapply(thresholds,function(t) { confusion.matrix(values,labels,t) })

  return(d)
}

precisions.recalls.fast <- function(pos.values,neg.values) {
  ## 0 for negative, 1 for positive
  l = rep(c(0,1),c(length(neg.values),length(pos.values)))
  v = c(neg.values,pos.values)

  ## sort by descending value
  indices = order(v,decreasing=TRUE)
  labels = l[indices]
  values = v[indices]

  tp = 0
  fp = 0

  z = 1
  precisions = rep(0,length(indices)+1)
  recalls = rep(0,length(indices)+1)

  tp = 0
  fp = 0

  for (i in 1:length(indices)) {
    if (labels[i]==1) {
      tp = tp + 1
    }
    else {
      fp = fp + 1
    }

    ## make sure not between tied values
    if (i==length(indices) || values[i]>values[i+1]) {
      ## can put a threshold

      ## if first, insert (r=0,p=1) point
      if (z==1) {
        recalls[z] = 0
        precisions[z] = 1
        z = z + 1
      }
      recalls[z] = tp/(length(pos.values))
      precisions[z] = tp/(tp+fp)
      z = z + 1
    }
  }

  ## truncate recals and precisions
  recalls = recalls[1:(z-1)]
  precisions = precisions[1:(z-1)]

  return (list(recalls=recalls,precisions=precisions))
}


# return c(v1[1],v2[1],v1[2],v2[2],...)  works with uneven vectors,
# appending the extra elements from the longer vector so return always
# has length of length(v1)+length(v2) from
# http://tolstoy.newcastle.edu.au/R/help/06/03/22717.html
interleave <- function(v1,v2) {
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}


#################################################################
## Estimators
#################################################################
# Under, connects lowest precision on left side to highest precision
# on the right
prcurve.lowertrap <- function(pos.values,neg.values,conf.level=0.95,conf.int.method="binomial") {
  t = precisions.recalls.fast(pos.values,neg.values)

  recalls = t$recalls
  precisions = t$precisions


  # get max and min precision for each unique recall
  r.unique = unique(recalls)
  p.min = sapply(r.unique,function(r) { min(precisions[recalls==r])})
  p.max = sapply(r.unique,function(r) { max(precisions[recalls==r])})

  # use min on left side and max on right side of each area between known recalls
  ids=2:length(r.unique)
  area = as.double((r.unique[ids]-r.unique[ids-1]) %*% (p.max[ids] + p.min[ids-1]) / 2)

  rs = rep(r.unique,each=2)
  ps = interleave(p.max,p.min)

  ci = aucpr.conf.int(area,
                      pos.values=pos.values,
                      neg.values=neg.values,
                      prcurve.estimator=prcurve.lowertrap,
                      conf.level=conf.level,
                      method=conf.int.method)

  return(list(area=area,
              x=rs,
              y=ps,
              conf.int=ci))
}

# connects highest precisions at each recall, definitely an
# over-estimate
prcurve.uppertrap <- function(pos.values,neg.values,conf.level=0.95,conf.int.method="binomial") {
  t = precisions.recalls.fast(pos.values,neg.values)

  recalls = t$recalls
  precisions = t$precisions

  # get max and min precision for each unique recall
  r.unique = unique(recalls)
  p.min = sapply(r.unique,function(r) { min(precisions[recalls==r])})
  p.max = sapply(r.unique,function(r) { max(precisions[recalls==r])})

  # use max precision for all recalls (NOT LEGITIMATE!)

  ids=2:length(r.unique)
  area = as.double((r.unique[ids]-r.unique[ids-1]) %*% (p.max[ids] + p.max[ids-1]) / 2)


  ci = aucpr.conf.int(area,
                      pos.values=pos.values,
                      neg.values=neg.values,
                      prcurve.estimator=prcurve.uppertrap,
                      conf.level=conf.level,
                      method=conf.int.method)


  return(list(area=area,
              x=r.unique,
              y=p.max,
              conf.int=ci))
}


## estimate PR curve and area under PR curve using average precision
## method, mean precision at each positive example
prcurve.ap.slow <- function(pos.values, neg.values, conf.level=0.95, conf.int.method="binomial") {
  ## thresholds, sorting to have increasing recall
  ## only use thresholds for positives
  thresholds = sort(pos.values,decreasing=TRUE)

  labels = c(rep(1,length(pos.values)),rep(0,length(neg.values)))
  values = c(pos.values,neg.values)

  ## create set of confusion matrices
  d=sapply(thresholds,function(t) { confusion.matrix(values,labels,t) })
  recalls = unlist(d[6,])
  precisions = unlist(d[5,])

  area = mean(precisions)


  y = rep(precisions,each=2)
  x = c(0,rep(recalls[1:(length(recalls)-1)],each=2),recalls[length(recalls)])

  ci = aucpr.conf.int(area,
                      pos.values=pos.values,
                      neg.values=neg.values,
                      prcurve.estimator=prcurve.ap.slow,
                      conf.level=conf.level,
                      method=conf.int.method)

  return (list(area=area,
               x=x,
               y=y,
               conf.int=ci))
}

prcurve.ap <- function(pos.values, neg.values, conf.level=0.95, conf.int.method="binomial") {
  ## 0 for negative, 1 for positive
  l = rep(c(0,1),c(length(neg.values),length(pos.values)))
  v = c(neg.values,pos.values)

  ## sort by descending value
  indices = order(v,decreasing=TRUE)
  labels = l[indices]
  values = v[indices]


  tp = 0
  fp = 0

  rs = rep(0,length(pos.values))
  ps = rep(0,length(pos.values))
  z = 1
  cur.tp = 0
  cur.fp = 0
  for (i in 1:length(indices)) {
    if (labels[i]==1) {
      cur.tp = cur.tp + 1
    }
    else {
      cur.fp = cur.fp + 1
    }
    ## make sure not between tied values
    if (i==length(indices) || values[i]>values[i+1]) {
      ## can put a threshold here
      tp = tp + cur.tp
      fp = fp + cur.fp

      if (cur.tp > 0) {
        for (j in 1:cur.tp) {
          rs[z] = tp/length(pos.values)
          ps[z] = tp/(tp+fp)
          z = z + 1
        }
      }
      cur.tp = 0
      cur.fp = 0
    }

  }


  if (z != length(pos.values)+1) {
    cat("WARNING: did not fill recall and precision arrays correctly (z: received=",z,", expected=",length(pos.values)+1,")\n",sep="")
  }
  area = mean(ps)

  y = rep(ps,each=2)
  x = c(0,rep((rs[1:(length(rs)-1)]),each=2),rs[length(rs)])


  ci = aucpr.conf.int(area,
                      pos.values=pos.values,
                      neg.values=neg.values,
                      prcurve.estimator=prcurve.ap,
                      conf.level=conf.level,
                      method=conf.int.method)

  return (list(area=area,
               x=x,
               y=y,
               conf.int=ci))
}

## estimate PR curve and area under PR curve with MLE parameters for
## positive and negative normal distributions
prcurve.binormal <- function(pos.values,neg.values,conf.level=0.95,conf.int.method="binomial") {
  ## Assume X ~ N(0,1), i.e. negatives have standard normal distribution.
  ## So only have 2 df to calculate.

  if (length(pos.values)<=1) {
    cat("ERROR: cannot use binormal estimate with fewer than 2 positive samples\n")
  }
  if (length(neg.values)<=1) {
    cat("ERROR: cannot use binormal estimate with fewer than 2 negative samples\n")
  }
  mean.hat = (mean(pos.values)-mean(neg.values))/sd(neg.values)
  sd.hat = sd(pos.values)/sd(neg.values)

  pi = length(pos.values)/(length(pos.values)+length(neg.values))

  area = binormal.area(pi=pi,pos.mean=mean.hat,pos.sd=sd.hat)

  ci = aucpr.conf.int(area,
                      pos.values=pos.values,
                      neg.values=neg.values,
                      prcurve.estimator=prcurve.binormal,
                      conf.level=conf.level,
                      method=conf.int.method)

  return(list(area=area,pi.hat=pi,pos.mean.hat=mean.hat,pos.sd.hat=sd.hat, conf.int=ci))
}


## Use linear interpolation in ROC space, translate to PR space and
## find the precision for the specified recall (r) for interpolating
## between PR points (r1,p1) and (r2,p2)
interpolate.point <- function(r1,p1,r2,p2,r) {
  res = r / (r*(1 + (1-p2)*r2/(p2*(r2-r1)) - (1-p1)*r1/(p1*(r2-r1))) + (1-p1)*r1/p1 - r1*(1-p2)*r2/(p2*(r2-r1)) + r1*(1-p1)*r1/(p1*(r2-r1)))
  res[r==0] = p2
  return(res)
  ## if (r==0) {
  ##   ## hmm, use limit which is just p2
  ##   return(p2)
  ## }
  ## else {
  ##   return(r / (r*(1 + (1-p2)*r2/(p2*(r2-r1)) - (1-p1)*r1/(p1*(r2-r1))) + (1-p1)*r1/p1 - r1*(1-p2)*r2/(p2*(r2-r1)) + r1*(1-p1)*r1/(p1*(r2-r1))))
  ## }
}

interpolate.curve <- function(recalls,precisions,num.samples=1000) {

  ## sort by increasing recall
  indices = order(recalls,decreasing=FALSE)
  rs = recalls[indices]
  ps = precisions[indices]

  xs = rep(0,num.samples)
  ys = rep(0,num.samples)
  ## point lies between rs[index-1] and rs[index]
  index = 1
  for (i in 1:num.samples) {
    ## evenly spaced in [0,1], including end-points
    recall = (i-1)/(num.samples-1)
    xs[i] = recall

    while (rs[index]<recall) {
      index = index + 1
    }

    if (index==1) {
      ## use (0,1) as beginning interp point (creates horizontal line
      ## from smallest (r,p) at y=p
      ys[i] = interpolate.point(0,1,rs[index],ps[index],recall)
    }
    else {
      ys[i] = interpolate.point(rs[index-1],ps[index-1],rs[index],ps[index],recall)
    }
  }
  return(list(x=xs,y=ys))
}

## Use ROC interpolation to calculate areas between PR points
## Assumes recalls has no duplicates and sorted in ascending order
## only calculates area from min(recalls) to max(recalls)
interpolate.area <- function(recalls, precisions) {
  area = 0.0

  for (i in 2:length(recalls)) {
    ## calculate area for recalls[i-1] to recalls[i]

    r1 = recalls[i-1]
    p1 = precisions[i-1]
    r2 = recalls[i]
    p2 = precisions[i]

    if (r1 == r2) {
      cat("WARNING: recalls passed to interpolate.area are not unique\n")
    }
    if (r1 > r2) {
      cat("WARNING: recalls passed to interpolate.area are not sorted ascending\n")
      cat("recalls: ",recalls,"\n")
      cat("precisions: ",precisions,"\n")
    }

    if (r1==0) {
      ## definite integral is undefined
      ## use rectangle with height p2 (as interpolate.curve creates)
      area = area + p2*(r2-r1)

    }
    else {
      ## formula between these is p' = r' / (a + b*r')
      a = (1-p1)*r1/p1 - r1*(1-p2)*r2/(p2*(r2-r1)) + r1*(1-p1)*r1/(p1*(r2-r1))
      b = 1 + (1-p2)*r2/(p2*(r2-r1)) - (1-p1)*r1/(p1*(r2-r1))

      ## indefinite integral is (bx - a log (a + bx))/b^2
      ## area contribution is definite integral from r1 to r2
      area = area + (b*r2 - a*log(a+b*r2))/(b*b) - (b*r1 - a*log(a+b*r1))/(b*b)
    }

  }
  return(area)
}

## estimate PR curve and AUCPR by interpolating between the max
## precisions for each recall
prcurve.interpolate <- function(pos.values,neg.values,conf.level=0.95,conf.int.method="binomial",aggregator = max) {
  t = precisions.recalls.fast(pos.values,neg.values)

  recalls = t$recalls
  precisions = t$precisions

  r.unique = unique(recalls)
  p.max = sapply(r.unique,function(r) { aggregator(precisions[recalls==r])})

  ## sort
  indices = order(r.unique,decreasing=FALSE)
  rs = r.unique[indices]
  ps = p.max[indices]

  area = interpolate.area(rs,ps)

  ci = aucpr.conf.int(area,
                      pos.values=pos.values,
                      neg.values=neg.values,
                      prcurve.estimator= function(p.v,n.v,c.l,c.i.m) { prcurve.interpolate(pos.values=p.v,neg.values=n.v,conf.level=c.l,conf.int.method=c.i.m,aggregator) },
                      conf.level=conf.level,
                      method=conf.int.method)

  curve = interpolate.curve(rs,ps)

  return(list(area=area,x=curve$x,y=curve$y,conf.int=ci))

}


roccurve.points <- function(pos.values,neg.values) {
  ## 0 for negative, 1 for positive
  l = rep(c(0,1),c(length(neg.values),length(pos.values)))
  v = c(neg.values,pos.values)

  ## sort by descending value
  indices = order(v,decreasing=TRUE)
  labels = l[indices]
  values = v[indices]

  tp = 0
  fp = 0

  z = 1
  tprs = rep(0,length(indices)+1)
  fprs = rep(0,length(indices)+1)

  for (i in 1:length(indices)) {
    if (labels[i]==1) { tp = tp + 1 }
    else { fp = fp + 1 }

    ## make sure not between tied values
    if (i==length(indices) || values[i]>values[i+1]) {
      ## can place threshold

      ## if first, insert the all negative point (tpr=0,fpr=0)
      if (z==1) {
        tprs[z] = 0
        fprs[z] = 0
        z = z + 1
      }
      tprs[z] = tp/length(pos.values)
      fprs[z] = fp/length(neg.values)
      z = z + 1
    }
  }

  ## truncate to used indices
  tprs = tprs[1:(z-1)]
  fprs = fprs[1:(z-1)]

  return (list(tprs=tprs,fprs=fprs))
}

## Estimate PR curve by taking convex hull in ROC space, convert only
## points on convex hull in ROC to PR, then use interpolation to go
## between.
prcurve.interpolate.convex <- function(pos.values, neg.values, conf.level=0.95,conf.int.method="binomial") {
  pi = length(pos.values)/(length(pos.values)+length(neg.values))

  t = roccurve.points(pos.values,neg.values)

  ## add the point (1,0) to make convex hull calculation clean (won't
  ## pick up any worse than random points)
  points = matrix(c(t$fprs,1,t$tprs,0),ncol=2)

  hull.indices = chull(points)

  ## remove the dummy point (0,1) we added, it's always the last element
  ## with index length(points)/2
  indices = hull.indices[!hull.indices==length(points)/2]
  points = points[indices,]

  ## due to bugs in chull and just to be safe sort points ascending by
  ## recall (second column)
  points = points[order(points[,2],decreasing=FALSE),]

  ## convert to PR space
  precisions = rep(0,length(indices))
  recalls = rep(0,length(indices))
  for (i in 1:length(indices)) {
    if (points[i,1]==0 & points[i,2]==0) {
      ## all negative
      precisions[i] = 1.0
      recalls[i] = 0.0
    }
    else {
      recalls[i] = points[i,2]
      precisions[i] = pi*points[i,2]/(pi*points[i,2] + (1-pi)*points[i,1])
    }
  }

  ## remove duplicate recalls, most likely having multiple recall=1.0
  ## points coming from the y=1 horizontal line in ROC space
  recalls.unique = unique(recalls)
  precisions.unique = sapply(recalls.unique,function(r) { max(precisions[recalls==r])})

  ## check that recalls are sorted ascending, having problems with non-sorted recalls
  if (is.unsorted(recalls.unique)) {
    cat("pos.values = ",pos.values,"\n")
    cat("neg.values = ",neg.values,"\n")
    cat("recalls = ",recalls,"\n")
    cat("precisions = ",precisions,"\n")
    cat("recalls.unique = ",recalls.unique,"\n")
    cat("precisions.unique = ",precisions.unique,"\n")

  }

  area = interpolate.area(recalls.unique,precisions.unique)

  ci = aucpr.conf.int(area,
                      pos.values=pos.values,
                      neg.values=neg.values,
                      prcurve.estimator=prcurve.interpolate.convex,
                      conf.level=conf.level,
                      method=conf.int.method)

  curve = interpolate.curve(recalls.unique,precisions.unique)

  return(list(area=area,x=curve$x,y=curve$y,conf.int=ci))

}

## stratified bootstrap
## resample positive and negative separately
bootstrap <- function(pos.values,neg.values,prcurve.function,replicates=1000) {
  areas = rep(0,replicates)

  for (i in 1:replicates) {
    p.v = sample(pos.values,replace=TRUE)
    n.v = sample(neg.values,replace=TRUE)

    res = prcurve.function(p.v,n.v)
    areas[i] = res$area
  }

  return(areas)
}


get_performance_measures = function(score = data_BED_PLANNING_training$TOTAL_SCORE, truth = data_BED_PLANNING_training$ADMIT_FLAG=="Admitted", thresholds = 50, round = TRUE){
  # pred <- prediction(score, true)
  # cutoffs = pred@cutoffs[[1]]
  # sens <- performance(pred, "sens")@y.values[[1]]
  # spec <- performance(pred, "spec")@y.values[[1]]
  # balanced_accuracy = (sens + spec)/2
  # ppv <- performance(pred, "ppv")@y.values[[1]]
  # npv <- performance(pred, "npv")@y.values[[1]]
  # return(list(sens = sens, spec = spec, balanced_accuracy = balanced_accuracy, ppv = ppv, npv = npv, cutoffs = cutoffs))

  data.table::data.table(t(sapply(thresholds, function(threshold){
    predict <- score >= threshold
    TNs = sum((truth==0)&(predict==0))
    TPs = sum((truth==1)&(predict==1))
    FPs = sum((truth==0)&(predict==1))
    FNs = sum((truth==1)&(predict==0))



    sens = TPs/(TPs + FNs)
    sens[is.nan(sens)] = NA
    spec = TNs/(TNs + FPs)
    spec[is.nan(spec)] = NA
    balanced_accuracy = (sens + spec)/2
    accuracy = (TNs + TPs)/length(truth)

    prevalence = sum(truth)/length(truth)
    ppv = (sens * prevalence)/(sens * prevalence + (1 - spec) * (1 - prevalence))
    ppv[is.nan(ppv)] = NA
    npv = (spec * (1-prevalence))/((1 - sens) * prevalence + spec * (1 - prevalence))
    npv[is.nan(npv)] = NA
    f1_score = (ppv * npv)*2/(ppv + npv)
    f1_score[is.nan(f1_score)] = NA

    youden = sens + spec - 1
    youden[is.nan(youden)] = NA


    if(round){

      sens = round(sens, 3)
      spec = round(spec, 3)
      balanced_accuracy = round(balanced_accuracy, 3)
      accuracy = round(accuracy, 3)
      f1_score = round(f1_score, 3)
      youden = round(youden, 3)

    }

    return(
      c(sens = sens, spec = spec, balanced_accuracy = balanced_accuracy, accuracy = accuracy, ppv = ppv, npv = npv, f1_score = f1_score, youden = youden, TNs = TNs, TPs = TPs, FPs = FPs, FNs = FNs, Ps = TPs + FNs, Ns = TNs + FPs)
    )
  })))







}





all_measures = function(score = data_BED_PLANNING_training$TOTAL_SCORE,
                        truth = data_BED_PLANNING_training$ADMIT_FLAG == 'Admitted',
                        thresholds = 1:100,
                        n_Boot = 1000,
                        use_Boot = FALSE,
                        BootID = 1:length(score),
                        confidence_level = 0.95){

  alpha = 1 - confidence_level

  pacman::p_load(ROCR)

  # score = data_BED_PLANNING_training$TOTAL_SCORE
  # truth = data_BED_PLANNING_training$ADMIT_FLAG == 'Admitted'
  # thresholds = 1:100
  # n_Boot = 1000
  # source("supporting functions.R")



  # get sens spec balanced_accuracy accuracy ppv npv f1_score youden
  performance_measures = get_performance_measures(
    score,
    truth,
    thresholds,
    round = FALSE
  )

  performance_measurements_lower_bound = performance_measurements_upper_bound = list()
  if(!use_Boot){
    # https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_One-Sample_Sensitivity.pdf
    var_sens = performance_measures$sens * (1 - performance_measures$sens) / performance_measures$Ps
    var_spec = performance_measures$spec * (1 - performance_measures$spec) / performance_measures$Ns

    z = qnorm(1-alpha/2)
    performance_measurements_lower_bound$sens = pmax(0,performance_measures$sens - z * sqrt(var_sens))
    performance_measurements_lower_bound$spec = pmax(0,performance_measures$spec - z * sqrt(var_spec))
    performance_measurements_lower_bound$ppv = pmax(0, performance_measures$ppv - z * sqrt(performance_measures$ppv * (1 - performance_measures$ppv) / (performance_measures$TPs + performance_measures$FPs)))
    performance_measurements_lower_bound$npv = pmax(0, performance_measures$npv - z * sqrt(performance_measures$npv * (1 - performance_measures$npv) / (performance_measures$TNs + performance_measures$FNs)))




    performance_measurements_upper_bound$sens = pmin(1, performance_measures$sens + z * sqrt(var_sens))
    performance_measurements_upper_bound$spec = pmin(1, performance_measures$spec + z * sqrt(var_spec))
    performance_measurements_upper_bound$ppv = pmin(1, performance_measures$ppv + z * sqrt(performance_measures$ppv * (1 - performance_measures$ppv) / (performance_measures$TPs + performance_measures$FPs)))
    performance_measurements_upper_bound$npv = pmin(1, performance_measures$npv + z * sqrt(performance_measures$npv * (1 - performance_measures$npv) / (performance_measures$TNs + performance_measures$FNs)))



    # https://stats.stackexchange.com/questions/475052/calculate-the-confidence-interval-of-a-balanced-accuracy-by-taking-the-mean-of-t
    performance_measurements_lower_bound$balanced_accuracy = pmax(0,performance_measures$balanced_accuracy - z * sqrt((var_sens + var_spec)/4))
    performance_measurements_upper_bound$balanced_accuracy = pmin(1, performance_measures$balanced_accuracy + z * sqrt((var_sens + var_spec)/4))


    performance_measurements_lower_bound$accuracy = pmax(0,performance_measures$accuracy - z * sqrt(performance_measures$accuracy * (1-performance_measures$accuracy)/(performance_measures$TNs[1] + performance_measures$TPs[1] + performance_measures$FNs[1] + performance_measures$FPs[1])))
    performance_measurements_upper_bound$accuracy = pmin(1,performance_measures$accuracy + z * sqrt(performance_measures$accuracy * (1-performance_measures$accuracy)/(performance_measures$TNs[1] + performance_measures$TPs[1] + performance_measures$FNs[1] + performance_measures$FPs[1])))





  }else{
    sample_indexes = list()
    for(b in 1:n_Boot){
      set.seed(b)

      sample_indexes[[b]] = sample(unique(BootID), replace = TRUE)

      # sample_indexes[[b]] = sample(1:length(score), replace = T)
    }
    # get performance for each bootstrap sample
    performance_measures_boot = list()


    for(b in 1:n_Boot){
      # print(b/n_Boot)

      boot_score = boot_truth = c()
      for(s in 1:length(sample_indexes[[b]])){
        boot_score = c(boot_score, score[BootID %in% sample_indexes[[b]][s]])
        boot_truth = c(boot_truth, truth[BootID %in% sample_indexes[[b]][s]])
      }

      performance_measures_boot[[b]] =
        get_performance_measures(boot_score, boot_truth, thresholds, round = FALSE)
    }
    # calculate the 2.5% and 97.5% confidence bounds
    measurements = colnames(performance_measures)
    for(m in 1:length(measurements)){
      measure = measurements[m]
      measure_per_threshold_for_each_boot = sapply(performance_measures_boot, function(x){
        unlist(x[[measure]])
      })

      performance_measurements_lower_bound[[measure]] = apply(measure_per_threshold_for_each_boot, 1, function(x){
        quantile(x, alpha/2, na.rm = TRUE)
      })

      performance_measurements_upper_bound[[measure]] = apply(measure_per_threshold_for_each_boot, 1, function(x){
        quantile(x, 1-alpha/2, na.rm = TRUE)
      })

    }
  }



  # calculate AUC and confidence bounds
  auc = performance(prediction(score, truth), "auc")@y.values[[1]]

  if(!use_Boot){
    ciAUC = pROC::ci(pROC::roc(truth, score), conf.level = confidence_level)
    if(auc<0.5){ # this performance() function does not distinguish control vs case. However, pROC::ci does. By doing this, pROC::ci always returns greater than 0.5. This part is to correct this.
      auc_lower_bound = 1 - as.numeric(ciAUC)[3]
      auc_upper_bound= 1 - as.numeric(ciAUC)[1]
    }else{
      auc_lower_bound = as.numeric(ciAUC)[1]
      auc_upper_bound= as.numeric(ciAUC)[3]
    }

  }else{
    auc_boot = c()
    for(b in 1:n_Boot){
      auc_boot[b] = performance(prediction(score[BootID %in% sample_indexes[[b]]], truth[BootID %in% sample_indexes[[b]]]), "auc")@y.values[[1]]
    }
    auc_lower_bound = quantile(auc_boot, alpha/2)
    auc_upper_bound = quantile(auc_boot, 1-alpha/2)
  }



  # calculate PR auc and confidence bounds
  pr_temp = prcurve.ap(score[truth==1], score[truth==0])
  pr = pr_temp$area
  if(!use_Boot){
    ciPR = pr_temp
    pr_lower_bound = ciPR$conf.int[1]
    pr_upper_bound = ciPR$conf.int[2]
  }else{
    pr_boot = c()
    for(b in 1:n_Boot){
      pr_boot[b] = prcurve.ap(score[BootID %in% sample_indexes[[b]]][truth[BootID %in% sample_indexes[[b]]]==1], score[BootID %in% sample_indexes[[b]]][truth[BootID %in% sample_indexes[[b]]]==0])$area
    }
    pr_lower_bound = quantile(pr_boot, alpha/2)
    pr_upper_bound = quantile(pr_boot, 1-alpha/2)
  }


  # calibration curve
  score_bin = cut(score, breaks = thresholds)
  percentage_of_truth_in_each_bin = by(truth, score_bin, function(x){sum(x)/length(x)})



  return(list(

    performance_measures = performance_measures,
    performance_measurements_lower_bound = performance_measurements_lower_bound,
    performance_measurements_upper_bound = performance_measurements_upper_bound,
    auc = auc,
    auc_lower_bound = auc_lower_bound,
    auc_upper_bound= auc_upper_bound,
    pr = pr,
    pr_lower_bound = pr_lower_bound,
    pr_upper_bound = pr_upper_bound,
    percentage_of_truth_in_each_bin = percentage_of_truth_in_each_bin,
    thresholds = thresholds


  ))
}

detach_package = function(package = "ctsc"){
  detach(paste0("package:",package), unload=TRUE)
}

from_measures_to_table = function(measures = all_measures, caption = "Performance Measurements using Different Threshold", include_TP_TF_FN_FP =FALSE){

  pacman::p_load(dplyr, kableExtra)

  thresholds = measures$thresholds

  sens = measures$performance_measures$sens
  sens_ci_low = measures$performance_measurements_lower_bound$sens
  sens_ci_up = measures$performance_measurements_upper_bound$sens

  spec = measures$performance_measures$spec
  spec_ci_low = measures$performance_measurements_lower_bound$spec
  spec_ci_up = measures$performance_measurements_upper_bound$spec

  ppv = measures$performance_measures$ppv
  ppv_ci_low = measures$performance_measurements_lower_bound$ppv
  ppv_ci_up = measures$performance_measurements_upper_bound$ppv

  npv = measures$performance_measures$npv
  npv_ci_low = measures$performance_measurements_lower_bound$npv
  npv_ci_up = measures$performance_measurements_upper_bound$npv

  if(!is.null(measures$performance_measurements_lower_bound$balanced_accuracy)){
    balanced_accuracy = measures$performance_measures$balanced_accuracy
    balanced_accuracy_ci_low = measures$performance_measurements_lower_bound$balanced_accuracy
    balanced_accuracy_ci_up = measures$performance_measurements_upper_bound$balanced_accuracy
  }else{
    balanced_accuracy = 0
    balanced_accuracy_ci_low = 0
    balanced_accuracy_ci_up = 0
  }
  if(!is.null(measures$performance_measurements_lower_bound$accuracy)){
    accuracy = measures$performance_measures$accuracy
    accuracy_ci_low = measures$performance_measurements_lower_bound$accuracy
    accuracy_ci_up = measures$performance_measurements_upper_bound$accuracy
  }else{
    accuracy = 0
    accuracy_ci_low = 0
    accuracy_ci_up = 0
  }



  TNs = measures$performance_measures$TNs
  TPs = measures$performance_measures$TPs
  FNs = measures$performance_measures$FNs
  FPs = measures$performance_measures$FPs

  # auc = measures$auc
  # auc_ci_low = measures$auc_lower_bound
  # auc_ci_up = measures$auc_upper_bound


  if(include_TP_TF_FN_FP){
    result = data.frame(thresholds,
                        TPs = TPs,
                        TNs = TNs,
                        FPs = FPs,
                        FNs = FNs,
                        sens = round(sens, digits = 3),
                        sens_ci = paste0("[", round(sens_ci_low, digits = 3), ", ", round(sens_ci_up, digits = 3), "]"),
                        spec = round(spec, digits = 3),
                        spec_ci = paste0("[", round(spec_ci_low, digits = 3), ", ", round(spec_ci_up, digits = 3), "]"),
                        ppv = round(ppv, digits = 3),
                        ppv_ci = paste0("[", round(ppv_ci_low, digits = 3), ", ", round(ppv_ci_up, digits = 3), "]"),
                        npv = round(npv, digits = 3),
                        npv_ci = paste0("[", round(npv_ci_low, digits = 3), ", ", round(npv_ci_up, digits = 3), "]"),
                        bal_acc = round(balanced_accuracy, digits = 3),
                        bal_acc_ci = paste0("[", round(balanced_accuracy_ci_low, digits = 3), ", ", round(balanced_accuracy_ci_up, digits = 3), "]"),
                        acc = round(accuracy, digits = 3),
                        acc_ci = paste0("[", round(accuracy_ci_low, digits = 3), ", ", round(accuracy_ci_up, digits = 3), "]"))

    colnames(result) = c("Thresholds",
                         "True Positive",
                         "True Negative",
                         "False Positive",
                         "False Negative",
                         "Sensitivity",
                         "Sensitivity CI",
                         "Specificity",
                         "Specificity CI",
                         "Positive Predictive Value",
                         "Positive Predictive Value CI",
                         "Negative Predictive Value",
                         "Negative Predictive Value CI",
                         "Balanced Accuracy",
                         "Balanced Accuracy CI",
                         "Accuracy",
                         "Accuracy CI")

    result = result %>% kable(format='html',align="ccc", caption=caption) %>% kable_classic(full_width = F)





  }else{
    result = data.frame(thresholds, sens = round(sens, digits = 3),
                        sens_ci = paste0("[", round(sens_ci_low, digits = 3), ", ", round(sens_ci_up, digits = 3), "]"),
                        spec = round(spec, digits = 3),
                        spec_ci = paste0("[", round(spec_ci_low, digits = 3), ", ", round(spec_ci_up, digits = 3), "]"),
                        ppv = round(ppv, digits = 3),
                        ppv_ci = paste0("[", round(ppv_ci_low, digits = 3), ", ", round(ppv_ci_up, digits = 3), "]"),
                        npv = round(npv, digits = 3),
                        npv_ci = paste0("[", round(npv_ci_low, digits = 3), ", ", round(npv_ci_up, digits = 3), "]"),
                        bal_acc = round(balanced_accuracy, digits = 3),
                        bal_acc_ci = paste0("[", round(balanced_accuracy_ci_low, digits = 3), ", ", round(balanced_accuracy_ci_up, digits = 3), "]"),
                        acc = round(accuracy, digits = 3),
                        acc_ci = paste0("[", round(accuracy_ci_low, digits = 3), ", ", round(accuracy_ci_up, digits = 3), "]"))


    colnames(result) = c("Thresholds",
                         "Sensitivity",
                         "Sensitivity CI",
                         "Specificity",
                         "Specificity CI",
                         "Positive Predictive Value",
                         "Positive Predictive Value CI",
                         "Negative Predictive Value",
                         "Negative Predictive Value CI",
                         "Balanced Accuracy",
                         "Balanced Accuracy CI",
                         "Accuracy",
                         "Accuracy CI")

    result = result %>% kable(format='html',align="ccc", caption=caption) %>% kable_classic(full_width = F)
  }




  # if(!is.null(measures$performance_measurements_lower_bound$balanced_accuracy)){
  #   result$bal_acc = result$bal_acc_ci = "-"
  # }
  #
  # if(!is.null(measures$performance_measurements_lower_bound$accuracy)){
  #   result$acc = result$acc_ci = "-"
  # }

  return(result)





}


from_any_table_to_rmarkdown = function(table, caption = ""){
  pacman::p_load(kableExtra, dplyr)
  as.data.frame(table) %>%
    kbl(caption = caption) %>%
    kable_paper("hover", full_width = F)
}

from_any_data_frame_to_rmarkdown_using_DT = function(df, caption = "", pageLength = 20){
  # pacman::p_load(kableExtra, dplyr)
  # as.data.frame(table) %>%
  #   kbl(caption = caption) %>%
  #   kable_paper("hover", full_width = F)

  if(!identical(sum(abs(as.numeric(df[[1]]) - 1:nrow(df))),0)){

    df = cbind(data.frame(o = 1:nrow(df)), df)

  }

  DT::datatable(df, rownames = FALSE, options = list(
    pageLength = pageLength
  ), caption = caption)
}



from_measure_to_calibration_plot = function(measures = BED_PLANNING_training_measures, title = "Calibration Curve"){
  pacman::p_load(ggplot2)
  truth_in_each_bin = as.numeric(measures$percentage_of_truth_in_each_bin)
  truth_in_each_bin[is.na(truth_in_each_bin)] = 0
  df = data.frame(y = truth_in_each_bin,
                  x = factor(as.character(1:length(truth_in_each_bin)), levels = as.character(1:length(truth_in_each_bin))),
                  xlabels = names(measures$percentage_of_truth_in_each_bin))


  ggplot(df, aes(x=x, y=y)) +
    geom_col() +
    scale_x_discrete(labels = df$xlabels)+
    labs(y = "Percentage of prevalence", x = 'Bins', title = title) + theme_bw() + ylim(0,1)

}

from_measure_to_ROC_curve = function(measures = BED_PLANNING_training_measures, title = "ROC Curve", legend = TRUE, legend_position = 'bottom', xlab = '1 - Specificity', ylab = 'Sensitivity'){

  pacman::p_load(ggplot2)

  df = data.frame(x = 1-measures$performance_measures$spec,
                  y = c(measures$performance_measures$sens,
                        measures$performance_measurements_lower_bound$sens,
                        measures$performance_measurements_upper_bound$sens),
                  type = rep(c('values',"lower_bound","upper_bound"), each = length(measures$performance_measures$spec)))
  df$type = factor(df$type, levels = c('values',"lower_bound","upper_bound"))

  ggplot(data=df, aes(x=x, y=y, group = type)) +
    geom_line(aes(linetype = type))+
    geom_point(aes(size = type))+
    scale_linetype_manual(values=c("solid","dotted","dotted"))+
    scale_size_manual(values = c(1, 0, 0))+
    theme_bw() +
    ylim(0,1)+
    xlim(0,1)+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = 0.1, linetype = "dashed")+
    labs(title = title, x = xlab, y = ylab)+
    annotate(geom="text", x=0.5, y=0.1, label=paste0("AUC = ",round(measures$auc,3)," (",round(measures$auc_lower_bound,3),", ",round(measures$auc_upper_bound,3),")"),
             color="red") +
    theme(legend.position=ifelse(legend, legend_position,"none"))


}


from_measure_to_PR_curve = function(measures = BED_PLANNING_training_measures, title = "PR Curve", legend = TRUE, legend_position = 'bottom', xlab = 'Recall', ylab = 'Precision'){

  pacman::p_load(ggplot2)

  df = data.frame(x = measures$performance_measures$sens,
                  y = c(measures$performance_measures$ppv,
                        measures$performance_measurements_lower_bound$ppv,
                        measures$performance_measurements_upper_bound$ppv),
                  type = rep(c('values',"lower_bound","upper_bound"), each = length(measures$performance_measures$sens)))
  df$type = factor(df$type, levels = c('values',"lower_bound","upper_bound"))

  ggplot(data=df, aes(x=x, y=y, group = type)) +
    geom_line(aes(linetype = type))+
    geom_point(aes(size = type))+
    scale_linetype_manual(values=c("solid","dotted","dotted"))+
    scale_size_manual(values = c(1, 0, 0))+
    theme_bw() +
    ylim(0,1)+
    xlim(0,1)+
    # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = 0.1, linetype = "dashed")+
    labs(title = title, x = xlab, y = ylab)+
    annotate(geom="text", x=0.5, y=0.1, label=paste0("AUC = ",round(measures$pr,3)," (",round(measures$pr_lower_bound,3),", ",round(measures$pr_upper_bound,3),")"),
             color="red")+
    theme(legend.position=ifelse(legend, legend_position,"none"))


}




table_one_way_with_percentage = function(factor = data_BED_PLANNING_training$ADMIT_YN){
  number = paste0(table(factor), " (",round((table(factor)/sum(table(factor)))*100,1),"%)")
  names(number) = names(table(factor))
  class(number) = "table"
  number
}


read_data  = function(path = "D:\\Jennly Zhang MetaboliteSD_ver03\\Raw_Phenotype_UCDavis.xlsx", sheet  = 1){
  library(data.table)


  # path = "C:\\Users\\Sili\\Documents\\Github\\Bajaj_2_5_2019\\Serum\\Bajaj complete urine and serum data transposed 12.31.18pm 90day 6month.csv"

  if(grepl("xlsx",path)){
    data = readxl::read_excel(path, col_names = FALSE, sheet  = sheet )
    data = data.table(data)
  }else{
    data = fread(path)
    data[data=='']=NA
  }



  sample_col_range = min(which(!is.na(data[1,]))):ncol(data)
  sample_row_range = 1:min(which(!is.na(data[[1]])))
  compound_col_range = 1:(min(which(!is.na(data[1,]))))
  compound_row_range = (min(which(!is.na(data[[1]])))):nrow(data)

  p = t(data[sample_row_range,sample_col_range,with=F])
  colnames(p) = p[1,]
  p = p[-1,]
  p = p[,c(ncol(p),1:(ncol(p)-1))]
  p = data.table(p)

  p = as.data.table(p)

  colnames(p) = make.unique(colnames(p), sep = "_")
  if(!"label"%in%colnames(p)){
    stop("Cannot find 'label' in your data. Please check the data format requirement.")
  }
  if(sum(is.na(p$label))>0){
    p$label[is.na(p$label)] = "na"
  }




  f = data[compound_row_range,compound_col_range,with=F]
  colnames(f) = as.character(f[1,])
  f = f[-1,]
  f = f[,c(ncol(f),1:(ncol(f)-1)),with=F]

  f = as.data.table(f)
  colnames(f) = make.unique(colnames(f), sep = "_")
  if(sum(is.na(f$label))>0){
    f$label[is.na(f$label)] = "na"
  }


  e = data[compound_row_range, sample_col_range, with = F]
  colnames(e) = as.character(e[1,])
  colnames(e)[is.na(colnames(e))] = "na"
  e = e[-1,]

  e_cat = e
  colnames(e_cat) = make.unique(colnames(e_cat), sep = "_")
  e_cat$label[is.na(e_cat$label)] = "na"
  e_cat$label = f$label
  colnames(e_cat) = c("label",p$label)

  e_cat_matrix = as.matrix(e_cat[,-1,with=F])


  e = data.table(label = e$label, sapply(e[,-1,with=F], function(x){
    as.numeric(x)
  }))

  colnames(e) = make.unique(colnames(e), sep = "_")
  e$label[is.na(e$label)] = "na"
  e$label = f$label
  colnames(e) = c("label",p$label)


  e_matrix = data.matrix(e[,-1,with=F])

  return(list(p = p, f = f, e = e, e_matrix = e_matrix,e_cat_matrix = e_cat_matrix))
}


summary_stat = function(x, digits = 2, include_5_and_95 = FALSE){

  mean = mean(x, na.rm = TRUE)
  sd = sd(x, na.rm = TRUE)
  min =  min(x, na.rm = TRUE)

  first_quantile = quantile(x, .25, na.rm = TRUE)
  median = median(x, na.rm = TRUE)
  third_quantile = quantile(x, .75, na.rm = TRUE)
  max =  max(x, na.rm = TRUE)
  range = max - min


  if(include_5_and_95){
    five_perc = quantile(x, .05, na.rm = TRUE)
    ninty_five_perc = quantile(x, .95, na.rm = TRUE)

    result = round(c(Mean = mean, SD = sd, Min = min, "5th Quantile" = five_perc, "25th Quantile" = first_quantile, Median = median, "75th Quantile" = third_quantile, "95th Quantile" = ninty_five_perc, Max = max, Range = range), 2)

  }else{
    result = round(c(Mean = mean, SD = sd, Min = min, "25th Quantile" = first_quantile, Median = median, "75th Quantile" = third_quantile, Max = max, Range = range), 2)
  }

  names(result) = c("Mean", "SD", "Min", "25th Quantile", "Median", "75th Quantile", "Max", "Range")

  return(result)
}





table_two_way_with_sum = function(x,y, x_name = "", y_name = ""){
  executable_text = paste0("table(",ifelse(x_name=="", "x",paste0(x_name, " = x")),",",ifelse(y_name=="", "y)",paste0(y_name, " = y)")))
  tbl = eval(parse(text = executable_text))

  value = as.numeric(tbl)
  matrix = matrix(value, nrow = nrow(tbl))
  matrix = rbind(matrix, colSums(matrix))
  matrix = cbind(matrix, rowSums(matrix))
  data.frame = as.data.frame(matrix, check.names = F)

  rownames(data.frame) = c(paste0(ifelse(x_name=="", "", paste0(x_name, " - ")), rownames(tbl)), "Sum")
  colnames(data.frame) = c(paste0(ifelse(y_name=="", "", paste0(y_name, " - ")), colnames(tbl)), "Sum")

  return(data.frame)
}

rmarkdowncolortext <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}
