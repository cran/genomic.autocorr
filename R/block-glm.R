##' @import data.table
##' @import magrittr
##' @importFrom parallel mclapply
##' @importFrom reshape melt.matrix
##' @importFrom stats acf as.formula coef glm median quantile rnorm sd ts var
NULL

globalVariables(c(".","value","beta.025","beta.975","variable"))

##' acf.summary
##'
##' summarize the autocorrelation in 
##' @inheritParams block.glm
##' @param data data.table containing variables named in `variables` and `order.by`
##' @param variables character vector listing columns of `data` to be explored for autocorrelation
##' @param order.by optionally, order `data` by variables in character vector `order.by`
##' @param lag.max maximum block size to explore (default=100)
##' @export
##' @examples
##' ## simulate data with 10 repeated observations in a row - ie there
##' ## should be autocorrelation only within windows <= 10
##' library(data.table)
##' data <- genomic.autocorr:::.sim.data() 
##' summ <- acf.summary(data,c("x","y0","y1"),lag.max=20)
##'
##' ## plot it
##' df <- melt(summ,c("lag","variable"),variable.name="acf")
##' par(mfrow=c(2,1))
##' matplot(matrix(df[acf=="full",]$value,ncol=3),
##'         main="full",
##'         pch=c("x","o","+"),
##'         type="b")
##' abline(h=0,lty=2)
##' legend("bottomright",
##'        c("x","y0","y1"),
##'        pch = "xo+", col = 1:3)
##' matplot(matrix(df[acf=="partial",]$value,ncol=3),
##'         main="partial",
##'         pch=c("x","o","+"),
##'         type="b")
##' abline(h=0,lty=2)
##' legend("bottomright",
##'        c("x","y0","y1"),
##'        pch = "xo+", col = 1:3)
acf.summary <- function(data,variables,order.by=NULL,lag.max=100) {
    if(!is.null(order.by))
        data <- data[ order(data[[order.by]]), ]
    dropna <- function(x) x[ !is.na(x) ]
    SUMM <- mclapply(variables, function(x) {
        if(var(data[[x]],na.rm=TRUE)==0)
            return(NULL)
        my.ts=ts(dropna(data[[x]]),start=1)
        acf.a <- acf(my.ts,lag.max=lag.max,plot=FALSE) 
        acf.p <- acf(my.ts,type="p",lag.max=lag.max,plot=FALSE)
        dt <- data.table(lag=acf.p$lag,full=acf.a$acf[-1],partial=acf.p$acf)
        dt[,variable:=x]
        dt })
    nulls <- sapply(SUMM,is.null)
    if(any(nulls)) {
        warning("dropped ",sum(nulls)," non-varying variable(s):\n",paste(variables[which(nulls)],collapse=" "))
        SUMM <- SUMM[!nulls]
    }
    do.call("rbind", SUMM)
}
##' internal function to simulate data for examples
##'
##' @title internal function to simulate data for examples
##' @param n number of independent observations
##' @param m group size - number of times each independent observation is repeated
##' @param beta Y1 ~ N( beta * X, 1)
##' @return data.table with
##' Chr (always 1, possibly needed for bootstrap),
##' x (explanatory variable),
##' y1 (response variable related to x),
##' y0 (response variable unrelated to x)
##' name (unique name for each independent observation)
##' @author Chris Wallace
##' @rdname internalfunctions
.sim.data <- function(n=500,m=10,beta=0.2) {
    n <- 500 # 
    m <- 10 #
    message("simulating ",n," independent observations repeated in blocks of ",m)
    X <- rnorm(n) #arima.sim(model=list(ar=c(0.1, 0.2, 0.1)),n=n)
    Y1 <- rnorm(n,mean=beta * X)
    Y0 <- rnorm(n)
    data.table(Chr=1,
                       name=rep(1:n,each=m),
                       x=rep(X,each=m),
                       y1=rep(Y1,each=m),
                       y0=rep(Y0,each=m))
}


## ggplot(melt(acf,c("lag","mark")),
##        aes(x=lag,y=value,col=mark)) + geom_path() + facet_grid(variable~.)

##' block.glm
##' 
##' Regression models for genomic data often assume there is
##' independence between neighbouring genomic elements when, in
##' reality, there is spatial dependence.  This function implements a
##' block bootstrap method for estimating correct variances of
##' parameter estimates.
##'
##' Note that this function uses `mclapply` to parallelise the
##' bootstrapping.  Please set `mc.cores` to something sensible, eg
##' \code{options(mc.cores=10)}
##' if you have 10 cores.
##'
##' @param f.lhs character vector, left hand side of a formula, the
##'     model(s) to be fit will be defined by
##'     `paste(f.lhs, f.rhs, sep=" ~ ")`
##' @param f.rhs character string, right hand side of a formula
##' @param data data.table containing the columns referred to in f.lhs
##'     and f.rhs
##' @param order.by if not `NULL`, the name of a column in `data` on
##'     which it should be sorted
##' @param strat.by if not `NULL`, the name of a column in `data` on
##'     which it should be stratified before block sampling.  Eg, if
##'     you are considering genomic data, you should stratify by
##'     chromosome as there should be no spatial correlation between
##'     chromosomes
##' @param block.size size of blocks of contiguous observations that
##'     will be sampled for bootstrap estimation of variance of
##'     parameter estimates
##' @param B number of bootstrap estimates
##' @param ... other arguments passed to `glm()` (eg
##'     `family="binomial"`)
##' @return data.table giving the estimated effect ("beta") of each
##'     item in f.rhs on each item in f.lhs, together with block
##'     bootstrap estimates of confidence interval (beta.025,
##'     beta.975) and standard error (se.beta) and the number of
##'     bootstraps on which those estimates are based.
##' @author Chris Wallace and Oliver Burren
##' @export
##' @examples
##'
##' ## simulate data with 10 repeated observations in a row - ie there
##' ## should be autocorrelation only within windows <= 10
##' library(data.table)
##' data <- genomic.autocorr:::.sim.data(beta=0.2) 
##'
##' ## suppose we ignored the autocorrelation and look at the
##' ## confidence interval for the effect of x on y1
##' r1<-summary(glm(y1 ~ x, data=data))$coefficients
##' r1
##'
##' ## if we know the block structure, as here, we can see the
##' ## confidence interval is (inappropriately) much tighter than
##' ## if we used just independent observations
##' r2<-summary(glm(y1 ~ x, data=data[!duplicated(name),]))$coefficients
##' r2
##'
##' ## use block bootstrap - x should only have a significant effect
##' ## on y1 and the confidence interval around its effect should be
##' ## closer to r2, above
##' r <- block.glm(f.lhs=c("y0","y1"), f.rhs="x",data=data,block.size=20,B=200)
##' r
##'
##' ## compare the block bootstrap and model based confidence intervals for x on y1
##' results <- rbind(c(r1[2,1], r1[2,1]-1.96*r1[2,2], r1[2,1]+1.96*r1[2,2]),
##' c(r2[2,1], r2[2,1]-1.96*r2[2,2], r2[2,1]+1.96*r2[2,2]),
##' as.numeric(r[4,.(beta,beta.025,beta.975)]))
##' dimnames(results) <- list(c("standard, ignore blocked","standard, independent obs","bootstrap"),
##' c("beta","LCI","UCI"))
##' results
##' 
##' with(as.data.frame(results), {
##' plot(1:nrow(results), beta,ylim=c(min(c(-0.01,LCI)),max(UCI)),axes=FALSE,xlab="Method",
##' main="Comparison of confidence intervals around coefficient estimates")
##' segments(x0=1:nrow(results),y0=LCI,y1=UCI)
##' abline(h=c(0,0.2),lty="dotted")
##' axis(1,1:nrow(results),rownames(results))
##' axis(2)
##' text(x=c(3,3),y=c(0,0.2),labels=c("null","true"),adj=c(1.1,0))
##' box()
##' })
##' 
block.glm<-function(f.lhs,
                    f.rhs,
                      data,
                      order.by=NULL,
                      strat.by=NULL,
                      block.size=20,
                      B=200,
                      ...){

#  message(f)

    ## compute valid starting points for blocks
    if(!is.null(order.by))
        data<-data[order(data[[order.by]]),]
    if(!is.null(strat.by)) {
        dsplit <- split(1:nrow(data),as.character(data[[strat.by]]))
        dsplit <- dsplit[ sapply(dsplit,length)>=block.size ]
    } else {
        dsplit <- list(1:nrow(data))
    }
    dsplit.idx<- lapply(dsplit,function(i){
        ret<-i[1:(length(i)-block.size)]
        names(ret)<-NULL
        ret
    }) %>% unlist()

  ##   effect size estimates  
    f1 <- paste(f.lhs,f.rhs,sep=" ~ ")
    fun <- function(data,...) {
        tmp <- lapply(f1, function(f) {
            glm(as.formula(f), data=data, ...) %>% coef()
        })  %>% do.call("rbind",.)
        rownames(tmp) <- f.lhs
        tmp
    }
    effects <- fun(data,...) #, ...)


  ## we want to select n/block.size blocks
  sample.no<-floor(nrow(data)/block.size)
    ## bootstrap B
#    dsplit.idx <- seq(block.size,nrow(data),by=block.size)
  RESULTS<-mclapply(1:B,function(z){
#    cat(".")
      samp<-sample(dsplit.idx,sample.no,replace=TRUE)
    idx<-lapply(samp,function(x){
      seq(from=x,to=x+block.size-1,by=1)
    }) %>% unlist()
    fun(data[idx,],...)
  })

    byvars <- c("y","x")
    RESULTS %<>% lapply(., reshape::melt.matrix, varnames=byvars) %>%
    do.call("rbind",.)
    results <- as.data.table(RESULTS)
    setkeyv(results,byvars)
    results <- results[,.(beta.med=median(value),
                          beta.025=quantile(value,0.975),
                          beta.975=quantile(value,0.025),
                          se.beta=sd(value),
                          nboot=.N),
                       by=byvars]
  
  beta <- reshape::melt.matrix(effects,varnames=byvars)
  beta <- as.data.table(beta)
  setnames(beta,"value","beta")
    results <- merge(results,beta,by=byvars)
    results[,beta.025:=2*beta - beta.025]
    results[,beta.975:=2*beta - beta.975]
  return(results)  
}
