### generate X
ls()
n <- 100

set.seed(201)
X <- runif(n)

### generate Y

probY_given_x <- function(x){

    if(x < 0.2){
       return(0.9)
   }
    else{
        if(x < 0.8){
           return(0.2)
       }
        else{
            return(0.9)
        }

    }
}


probs <- sapply(X, probY_given_x)

Y <- (2*(probs > runif(n)))-1

plot(X, Y, col = Y+2, pch = 16)
