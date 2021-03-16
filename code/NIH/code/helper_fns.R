## Function which calculates (un)weighted AUC given a set of predictions, labels, cutpoints, and weights.
## Takes the arguments
##    response: predictions associated with each label
##    labels: Vector of 0/1 indicating event label (i.e. alive vs died, etc.)
##    cutpts: Vector of numeric values at which to estimate the receiver operating characteristic curve.
##    weights: Vector of weights for each observation. If unspecified, this will result in equal weights for all observations,
##             calculating un-weighted AUC.
calc_weighted_AUC <- function(response, labels, cutpts, weights=NULL){
    stopifnot(length(response) == length(labels))
    stopifnot(all(labels %in% c(0,1)))

    N      <- length(response)
    u_resp <- unique(sort(response))
    n_cuts <- length(cutpts)

    x1 <- x2 <- rep(NA, n_cuts)

    ## if weights argument is unspecified, set equal weights for all response/label combinations
    if(is.null(weights)) weights = rep(1,N)

    ## obtain normalized weights and response values for those individuals who did not experience the event
    f_0    <- weights[labels==0]/mean(weights[labels==0])
    resp_0 <- response[labels==0]

    ## obtain normalized weights and response values for those individuals who did experience the event
    f_1    <- weights[labels==1]/mean(weights[labels==1])
    resp_1 <- response[labels==1]

    ## evaluate the weighted TPR and FPR at each cutpoint
    for(j in 1:n_cuts){
        x1[j] <- mean((resp_1 >= cutpts[j])*f_1)
        x2[j] <- mean((resp_0 >= cutpts[j])*f_0)
    }

    ## estimate the weighted AUC as the integral of the estimated ROC curve
    ## using the trapezoidal rule
    sum((x1[1:(n_cuts-1)] + x1[2:(n_cuts)])/2 * (x2[2:(n_cuts)] - x2[1:(n_cuts-1)]))
}
