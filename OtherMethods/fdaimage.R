#! /usr/bin/env Rscript
rm(list=ls())

library(FDAimage)

# Y: image
# X: covariate
# Z: coordinates (2D only)
# V: vertices
# Tr: triangles (list of triplets of vertices)
fdaimage = function(Y, X, Z, V=NULL, Tr=NULL, boundary="chull",
         d=5, r=1, nfold=10, lambda=10^(seq(-6,6,by=1)),
         nboot=100, alpha0=0.05, return_sigma=F, verbose=F, debug=F){
    if(debug){
        save(Y, X, Z, file="tmp.RData")
    }
    n = dim(Y)[1]
    N = dim(Y)[2]
    if(is.null(V) | is.null(Tr)){
        c1 = 0.5
        nT = floor(min(c1 * n^(1/(2*d+2)) * N^(1/2), N / 10))
        nT = max(8, nT)
        nMesh = floor(sqrt(nT/2))
        if(boundary == "rectangle"){
            umin = min(Z[,1])
            umax = max(Z[,1])
            vmin = min(Z[,2])
            vmax = max(Z[,2])
            us = seq(umin, umax, length.out=nMesh+1)
            vs = seq(vmin, vmax, length.out=nMesh+1)
            boundary = rbind(cbind(umin, vs), cbind(umax, vs), cbind(us, vmin), cbind(us, vmax))
        } else if(boundary == "chull"){
            # use the convex hull of all the points
            boundary = Z[chull(Z),]
        }
        #pngfile = ifelse(debug, 'tmp.png', './dev/null')
        #png(pngfile)
        VTr = TriMesh(boundary, nMesh)
        #dev.off()
        V = VTr$V 
        Tr = VTr$Tr
        if(verbose){
            print(nT)
            print(nMesh)
            print(dim(boundary))
            print(dim(V))
            print(dim(Tr))
        }
    }
    cv = cv.FDAimage(Y, X, Z, V, Tr, d, r, lambda, nfold)
    lamc = cv$lamc
    if(return_sigma){
        fdaout = scc.FDAimage(Y, X, Z, V, Tr, d, r, lamc, alpha0, nboot)
        beta = fdaout$beta.hat
        sigma = fdaout$Sigma
        ci_lower = fdaout$cc.l
        ci_upper = fdaout$cc.u
        alpha_adj = fdaout$alpha.adj
    } else {
        fdaout = fit.FDAimage(Y, X, Z, V, Tr, d, r, lamc)
        beta = fdaout$beta
        sigma = NA
        ci_lower = NA
        ci_upper = NA
        alpha_adj = NA
    }
    out = list(
        beta=beta, sigma=sigma,
        ci_lower=ci_lower, ci_upper=ci_upper,
        alpha_adj=alpha_adj)
    return(out)
}
