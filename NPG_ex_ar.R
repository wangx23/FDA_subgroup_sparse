###### NPGex and NPGar #####

#### two steps #### first use NPG to estimate beta and then update working correlation matrix to run the algrithm again 

### NPGobj is the NPG obj; if obj is null, run NPG again 

NPGex = function(NPGobj = NULL, datind, dattime, X, y, diagD, B_0, index, gamma1, gamma2, 
                 theta, tau, n, p, 
                 max_iter, eps_abs, eps_rel)
{
  if(is.null(NPGobj)) ### if NPGobj is not provided, rerun the NPG
  {
    NPGobj = prclust_admm(X, y, diagD, B_0, index,
                 gamma1 = gamma1, gamma2 = gamma2, 
                 theta=theta, tau = tau, 
                 n = n, p = p,  max_iter= max_iter,
                 eps_abs=eps_abs, eps_rel=eps_rel)
  }
  
  Best = NPGobj$B   #### estimate of B
  
  #### residuals #####
  
  residualsdf = data.frame(ind = dat$ind, residuals = y - X %*% c(Best))
##### some problem with unblanced setting ####
  
}