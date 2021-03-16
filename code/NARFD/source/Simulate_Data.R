# Simulate data using two functional prototypes
 # Inputs
  # N is the number of functional observations
  # D is the number of grid points
  # type can be "GFPCA" or "NARFD"
  # seed is the random seed
 # Outputs
  # returns a list with elements
    # observations: matrix of dimension DxN of N functional observations
    #               on a grid of length D
    # trajectories: matrix of dimension DxN of means used to simulate
    #               observations using the Poisson distribution
    # protoypes: functional prototypes, matrix of dimension Dx2
    # scores: scores, matrix of dimension 2xN

Simulate <- function(N = 20, D = 144, type = "NARFD", seed){
  set.seed(seed)
  Kp <- 2
  s <- seq(from=0, to=2*pi, length.out=D)
  if (type=="NARFD"){      # for NARFD everything needs to be positive
    FPC1 <- sin(2 * s) + 1
    FPC2 <- cos(s) + 1
  } else if (type=="GFPCA"){
    FPC1 <- sin(2 * s)
    FPC2 <- cos(s)
  } else {
    stop("Only type=NARFD or GFPCA supported")
  }
  FPCs <- cbind(FPC1, FPC2)
  colnames(FPCs) <- NULL
  if (type=="NARFD"){
    score.sds <- c(4, 3)
    overall.mean <- matrix(s * 0) 
  } else if (type=="GFPCA"){
    score.sds <- c(1.5, 1.0)
    overall.mean <- matrix(rep(3, length(s)))
  }
  scores <- matrix(rnorm(N * Kp, mean=0, sd=score.sds), nrow=2)
  if (type=="NARFD"){
    scores <- scores ^ 2
  }
  overall.mean.matrix <- overall.mean %*% matrix(rep(1, N), nrow=1)
  FPCs.matrix <- FPCs %*% scores
  subject.mean.matrix <- overall.mean.matrix + FPCs.matrix
  if (type=="GFPCA"){
    subject.mean.matrix <- exp(subject.mean.matrix)
  }
  observations <- matrix(rpois(N * D, lambda=subject.mean.matrix), nrow=D, ncol=N)
  list(observations=observations, trajectories=subject.mean.matrix, 
       scores=scores, prototypes=FPCs)
}

# converts wide data frame with observations in columns
# and converts to long data frame with columns .id, .observed and .index
# .id is character-valued, starting with I
# .index is the time index, starting with 1 and going to number of rows
# in wide
MakeLong <- function(wide){
  D <- nrow(wide)
  
  long <- as.data.frame(wide) %>% 
    mutate(.index = 1:D) %>% 
    gather(.id, .observed, -.index) %>%
    mutate(.id = paste0("I", .id)) %>% 
    as_tibble()
  
  long
}