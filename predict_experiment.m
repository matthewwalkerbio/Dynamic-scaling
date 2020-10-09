function [S, U] = predict_experiment(Re, Cd, data, constants, spline_fit_V, varargin)

% this is used for a prediction of tank S and U, either from the sphere
% curve if we don't yet have at least 2 data points, or for
% future predictions once you have an experimental Cd vs Re curve

% always input an Re that you want to have in tank
% Cd is either [] if this is the first ever prediction (in which case
% sphere curve is assumed) or, the Cd that you want to have in tank
% (probably output from predict_real_life since Re, Cd should match between
% real life, tank)

% if varargin is input, it is true if we want to display stuff, false if
% not.  defaults to true if not input.

% before you have any data, need to make a guess for first experiment
% input a guessed Re


if isempty(Cd)  %use sphere equation to guess Cd if its not input
    data.Re = [];  data.Cd = [];
    Cd = Cd_vs_Re(Re, data);
    
    disp('Using sphere eq');
end

if isempty(spline_fit_V)  % no spline available, use V_model = V_real * S^3 like we used to do
    coeffs = fliplr([ 1/2*Cd*data.real.A*Re^2*constants.mu^2/constants.rho_f/data.real.L^2   ,    0,   0,    data.real.V*constants.g*(constants.rho_f - constants.rho_m)   ]);  %shatlab wants it from high to low order terms
    temp = roots(coeffs);
    answer_inds = find(imag(temp) == 0);  %hopefully exactly one real root
    if length(answer_inds) ~= 1
        disp('shit, problem with root finding')
    else
        S = temp(answer_inds);
    end
else % use empirical relationship between V_model and S to predict V_model from S
    S_guess = 15;  % guess for S prediction for model to match real life Re, Cd
    fun = @(S) 1/2*Cd*data.real.A*Re^2*constants.mu^2/constants.rho_f/data.real.L^2   +   constants.g*(constants.rho_f - constants.rho_m) * slmeval(S, spline_fit_V)   ;
    
    S =  fzero(fun, S_guess);
    
end



U = Re*constants.mu/constants.rho_f/S/data.real.L;

if isempty(varargin) || varargin{1}
    disp(['Given input Re = ',num2str(Re,5),' and Cd = ',num2str(Cd,5),',     suggested scale factor S = ',num2str(S,5),'   and predicted tank speed U = ',num2str(U,5),'   m/s']);
end