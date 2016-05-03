multihyper <- function(x,n,M,N)
{
	return (pmultihyper(c(x,n-sum(x)),n,c(M,N-sum(M)),N))
}