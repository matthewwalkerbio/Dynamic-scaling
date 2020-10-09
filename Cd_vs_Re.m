function [Cd_Inf, spline_fit] = Cd_vs_Re(Re,data, fitting, varargin)

% Re is Reynolds number
% data is data.Re and data.Cd, either both empty [] (in which case sphere curve is used) or vectors of
% experimental values (in which case a monotonically decreasing cubic spline is fit and sampled at the
% input Re)

% fitting.n_knots is number of knots for fitted spline - the higher, the better fit
% to the data points but the more possibility for unwanted fitting behavior between
% data points. These options are specified in main script.

% optional varargin can contain preexisting polynomial fit coeffs, to save time
% needlessly recomputing it for the same input data as before

% outputs:

% Cd is Cd sampled from either sphere curve or empirical polynomial at the
% input Re
% coeffs is either empty or the polynomial coefficients from the fit to the
% data, for possible reuse later

spline_fit = [];  %initialize so it's always defined

if  ~isempty(data) && length(data.Re) >= 3  % we have at least 3 data points, fit a spline
    if isempty(varargin)
        if fitting.n_knots >= 3
            knot_placement = fitting.knot_placement;
        else
            knot_placement = 'fixed'; %for 2 knots, no choice
        end
        options = slmset('decreasing','on','extrapolation','linear','knots',fitting.n_knots,'degree',fitting.degree,'interiorknots',knot_placement,'concaveup',fitting.concave_up);
%         options = slmset('decreasing','on','extrapolation','linear','knots',fitting.n_knots,'degree','linear','interiorknots',knot_placement,'concaveup',fitting.concave_up);
%         spline_fit = slmengine(data.Re,data.Cd_Inf,options);
           spline_fit = slmengine(data.Re,data.Cd_Inf,options);     
    else
        spline_fit = varargin{1};
    end
    
    if ~isnan(Re)
        Cd_Inf = slmeval(Re, spline_fit);
    else
        Cd_Inf = NaN;
    end
    
else % Don't have any data, so just use Morrison's (2013) sphere curve, should be valid for Re  = 0 to Re  = inf
    
    Cd_Inf = Morrison(Re);
    
end

