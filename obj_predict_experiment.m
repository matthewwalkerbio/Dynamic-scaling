function [obj_zero, lambda, Cd_tank, U] = obj_predict_experiment(Cd_Inf, Re, S_guess, data, constants, predict_V, do_correction)
% This predicts the outcome of the next experiment to be ran, taking the
% inputs from the real life predictions (predict_real_life.m)


if ~isempty(Cd_Inf) && do_correction
% need to correct Cd_Inf for wall effects in tank
lambda = S_guess * data.real.L ./ constants.tank_dia;
[Cd_tank] = Cd_Inf2Tank(Cd_Inf, Re, lambda);
 
else % using sphere eq or ignoring wall effects
    lambda = NaN;
    Cd_tank = Cd_Inf;
end

obj_zero = sum_forces(Re, Cd_tank, S_guess, data, constants, predict_V, false); % Inputs from predict_real_life
 U = Re*constants.mu/constants.rho_f/S_guess/data.real.L; % predict the sinking speed (U) of the model

