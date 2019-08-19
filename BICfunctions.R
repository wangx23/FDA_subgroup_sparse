## first step 
BIC1 = function(Best, y, id, X, basis, n, gamma1, D)
{
  # Best is from the output a matrix
  yest = X %*% c(Best)
  ni = table(id)
  Dg1 = gamma1 * D
  dfv = rep(0,n)
  for(i in 1:n)
  {
    Bmi = basis[[i]]
    dfv[i] = sum(diag(Bmi %*% solve(t(Bmi) %*% Bmi + Dg1) %*% t(Bmi)))
  }
  
  BICvalue1 = sum(log(as.numeric(by((yest - y)^2,id,sum))/ni) + log(ni)/ni * dfv)
  return(BICvalue1)
}


### second step ###

# n is the number of subjects
BIC2 = function(obj,y, id, X, basis, n, gamma1, D)
{
  yest = X %*% c(obj$B)
  ni = table(id)
  Dg1 = gamma1 * D
  
  dfv = rep(0,n)
  for(i in 1:n)
  {
    Bmi = basis[[i]]
    dfv[i] = sum(diag(Bmi %*% solve(t(Bmi) %*% Bmi + Dg1) %*% t(Bmi)))
  }
  
  Ad_final <- create_adjacency(obj$V, n);
  G_final <- graph.adjacency(Ad_final, mode = 'upper')
  #clustering membership
  cls_final <- components(G_final);
  #number of clusters
  k_final <- cls_final$no;
  
  BICvalue2 = log(sum((yest-y)^2)/n) + log(n)/n * k_final*mean(dfv)
  return(BICvalue2)
}