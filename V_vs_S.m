function [V] = V_vs_S(S, spline_fit_V, data)
% wrapper V vs S function, uses spline fit within the bounds of the
% original data and extrapolates beyond that assuming a cubic scaling of V
% vs S


persistent alpha_lo alpha_hi S_min S_max % These values are kept as they're constant

if isempty(alpha_lo) % Need to calculate alpha_lo, alpha_hi if this the first time the script is being run
    [S_min] = min(data.S);
    V_min = slmeval(S_min, spline_fit_V);  %this is not the measured value, but the value of the spline at S_min, so that the resulting V(S) function is continuous
    alpha_lo = V_min / S_min^3; % proportionality constant, i.e. V = alpha_lo * S^3 for S < S_min
    
    [S_max] = max(data.S);
    V_max = slmeval(S_max, spline_fit_V);  %this is not the measured value, but the value of the spline at S_min, so that the resulting V(S) function is continuous
    alpha_hi = V_max / S_max^3; % proportionality constant, i.e. V = alpha_lo * S^3 for S < S_min
end

% evaluate predicted V at input S
V = NaN(size(S));
V( S > S_min & S < S_max ) = slmeval(S( S > S_min & S < S_max ), spline_fit_V); % above S_min and below S_max, use spline
V( S <= S_min ) =  alpha_lo * S( S <= S_min ).^3; % below left limit of data
V( S >= S_max ) =  alpha_hi * S( S >= S_max ).^3; % above right limit of data