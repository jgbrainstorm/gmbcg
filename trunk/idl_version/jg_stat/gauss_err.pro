; this is the same as gauss.pro for historical reason
function gauss_err,x,x_err,mu,sigma

     tau2=x_err^2.+sigma^2.
     res=1./sqrt(2.*!pi*tau2)*exp(-(x-mu)^2./tau2/2.)

     return,res
end
