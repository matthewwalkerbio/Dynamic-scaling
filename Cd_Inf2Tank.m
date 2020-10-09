function [Cd_tank] = Cd_Inf2Tank(Cd_Inf, Re, lambda)
% Convert Cd with no walls (Cd_Inf) to Cd measured with walls(Cd_tank)
Cd_tank = Cd_Inf + 24./Re * (K_HS(lambda) - 1);