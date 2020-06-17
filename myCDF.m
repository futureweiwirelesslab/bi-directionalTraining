function [x x_cdf] = myCDF(inVect)
%% generate cdf data
%% [x x_cdf] = myCDF(inVect)
% FUTUREWEI 2019, MIT license

inVect = sort(inVect);
x = unique(inVect);

if length(x)==1
  x_cdf = 1;
else
  x_pdf = hist(inVect,x);
  x_cdf = cumsum(x_pdf)/length(inVect);
end
