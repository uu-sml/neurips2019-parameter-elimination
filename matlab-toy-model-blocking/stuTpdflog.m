function [p] = stuTpdflog(z,mu,sig2,nu)
% Computes the log pdf of a student t distribution
p = log(gamma( (nu+1)/2 )) - log(gamma(nu/2)) - 0.5*log(pi*nu*sig2) -(0.5*(nu+1))*log(1 + (z-mu).^2./(nu*sig2));
end

