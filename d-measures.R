# From https://github.com/tischieber/Quantifying-Network-Structural-Dissimilarities/blob/master/d_distance.r
# and this paper - https://www.nature.com/articles/ncomms13928

# functions from hypoweavr by Grames - https://github.com/elizagrames/hypoweavr/blob/main/R/functions.R

#' Calculate graph dissimilarity
#' @description Calculates dissimilarity between two graphs based on topological distance
#' @param x a graph object or list of graph objects
#' @param y a graph object (if x is a graph object)
#' @return the dissimilarity metric from Scheiber et al. 2017
calc_dissimilarity <- function(x, y=NULL){
    if(!is.null(y)){
        x <- list(x, y)
    }
    changes <- rep(NA, length(x))
    used <- which(!unlist(lapply(x, is.null)))
    for (i in used[2]:length(x)) {
        if(i %in% used){
            previous_valid <- max((1:(i-1))[1:(i-1) %in% used])
            changes[i] <- gDist(x[[i]], x[[previous_valid]], 0.45, 0.45, 1)
        }
    }
    if(length(x)==2){
        changes <- changes[2]
    }
    
    return(changes)
}







#' Graph distance
#' @description See Schieber et al. (2017) for details
#' @param g a graph object to compare to h
#' @param h a graph object to compare to g
#' @param w1 weight of term1
#' @param w2 weight of term2
#' @param w3 weight of term3
#' @return graph distance
gDist <-function(g,h,w1,w2,w3){
    
    # Network dissimilarity functions from Schieber et al. (2017)
    # https://github.com/tischieber/Quantifying-Network-Structural-Dissimilarities
    
    # Shannon entropy
    entropia<-function(a) {
        a<-a[which(a>0)];
        -sum(a*log(a));
    }
    
    # Node distance
    node_distance<-function(g){
        n<-length(igraph::V(g))
        if(n==1){
            retorno=1
        }
        
        if(n>1){
            a<-Matrix::Matrix(0,nrow=n,ncol=n,sparse=TRUE)
            m<-igraph::shortest.paths(g,algorithm=c("unweighted"))
            m[which(m=="Inf")]<-n
            quem<-setdiff(intersect(m,m),0)
            
            for(j in (1:length(quem))){
                l<-which(m==quem[j])/n
                linhas<-floor(l)+1
                posicoesm1<-which(l==floor(l))
                if(length(posicoesm1)>0){
                    linhas[posicoesm1]<-linhas[posicoesm1]-1
                }
                a[1:n,quem[j]]<-graphics::hist(linhas,plot=FALSE,breaks=(0:n))$counts
            }
            #m<-c()
            retorno=(a/(n-1))
            retorno <- as.matrix(retorno) # was returning a sparse matrix
        }
        return(retorno)
    }
    
    # Network node dispersion
    nnd<-function(g){
        N<-length(igraph::V(g))
        nd<-node_distance(g)
        pdfm<-colMeans(nd)
        norm<-log(max(c(2,length(which(pdfm[1:(N-1)]>0))+1)))
        return(c(pdfm,max(c(0,entropia(pdfm)-entropia(nd)/N))/norm))
    }
    
    # Alpha centrality
    alpha<-function(g){
        N<-length(igraph::V(g))
        r<-sort(igraph::alpha_centrality(g,exo=igraph::degree(g)/(N-1),alpha=1/N))/((N^2))
        return(c(r,max(c(0,1-sum(r)))))
    }
    
    first<-0
    second<-0
    third<-0
    g<-g
    h<-h
    N<-length(igraph::V(g))
    M<-length(igraph::V(h))
    PM<-matrix(0,ncol=max(c(M,N)))
    
    if(w1+w2>0){
        pg=nnd(g)
        PM[1:(N-1)]=pg[1:(N-1)]
        PM[length(PM)]<-pg[N]
        ph=nnd(h)
        PM[1:(M-1)]=PM[1:(M-1)]+ph[1:(M-1)]
        PM[length(PM)]<-PM[length(PM)]+ph[M]
        PM<-PM/2
        
        first<-sqrt(max(c((entropia(PM)-(entropia(pg[1:N])+entropia(ph[1:M]))/2)/log(2),0)))
        second<-abs(sqrt(pg[N+1])-sqrt(ph[M+1]))
    }
    
    if(w3>0){
        pg<-alpha(g)
        ph<-alpha(h)
        m<-max(c(length(pg),length(ph)))
        Pg<-matrix(0,ncol=m)
        Ph<-matrix(0,ncol=m)
        Pg[(m-length(pg)+1):m]<-pg
        Ph[(m-length(ph)+1):m]<-ph
        third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))/2
        g<-igraph::complementer(g)
        h<-igraph::complementer(h)
        
        pg<-alpha(g)
        ph<-alpha(h)
        m<-max(c(length(pg),length(ph)))
        Pg<-matrix(0,ncol=m)
        Ph<-matrix(0,ncol=m)
        Pg[(m-length(pg)+1):m]<-pg
        Ph[(m-length(ph)+1):m]<-ph
        third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))/2
    }
    
    return(w1*first+w2*second+w3*third)
    
}


