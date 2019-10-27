function r = acf(d,maxlag)
% Computes the autocorrelation function of a sequence d 
% Input: Sequence (d), lagest lag (maxlag)
%Output: Autocorrelation of the seqence d

mu = mean(d);   
r = xcorr((d-mu),maxlag,'coeff');
r = r(maxlag+1:end);


