FUNCTION gauss,x,x_err,mu,sigma
           
         gau = 1.0/sqrt(2*!pi*(sigma^2.+x_err^2.))*exp(-(x-mu)^2/(2.*(sigma^2+x_err^2)))
         
         RETURN,gau
END

