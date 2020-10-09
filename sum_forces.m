function [obj_zero] = sum_forces(Re, Cd, S_guess, data, constants, predict_V, varargin)

% Generates a prediction of the scale and sinking speed (S and U) for tank
% conditions, either using the sphere curve (Morrison, 2013) if there are
% fewer than 3 data points (i.e. scales). If there is an experimental Cd vs
% Re curve, then this is used instead.

% The input Re is always what is wanted in the tank (i.e. the output value
% from predict_real_life). In the absence of data Cd (Cd = [])is predicted 
% from sphere curve, otherwise Cd is taken from predict_real_life, as Re 
% and Cd should match between tank and real life conditions 


% if varargin is input, it is true if we want to display stuff, false if
% not.  defaults to true if not input.

% if there is no data at all, then an estimated Re should be input



if isempty(Cd)  % Use sphere equation to estimate Cd, if no value is input
    data.Re = [];  data.Cd = [];
    Cd = Cd_vs_Re(Re, data);
    
    disp('Using sphere eq');
end

if isempty(predict_V)  % no spline available, use V_model = V_real * S^3, this will probably negatively impact results
    obj_zero =  1/2*Cd*data.real.A*Re^2*constants.mu^2/constants.rho_f/data.real.L^2   +   data.real.V*constants.g*(constants.rho_f - constants.rho_m)  * S_guess^3;
else % use empirical relationship between V_model and S to predict V_model from S
    % here we continue to assume that L and A do scale with S and S^2
    % respectively
    obj_zero =  1/2*Cd*data.real.A*Re^2*constants.mu^2/constants.rho_f/data.real.L^2   +   constants.g*(constants.rho_f - constants.rho_m) * predict_V(S_guess)   ;
    
    
    
end