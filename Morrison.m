function [Cd] = Morrison(Re)
% Calculated Cd for a sphere for Re = 0 to Re = inf. From Morrison (2013)
Cd = 24 ./ Re + 2.6*Re./5.0./(1+(Re./5.0).^1.52)  +  0.411.*(Re./263000).^(-7.94)./(1+(Re./263000).^(-8.00))  + Re.^(0.80)./461000;