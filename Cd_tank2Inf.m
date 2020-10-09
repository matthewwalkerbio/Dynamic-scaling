function [Cd_Inf, K] = Cd_tank2Inf(Cd_tank, Re, lambda)
K = K_HS(lambda) - 1; % Obtain K values, if required
Cd_Inf = Cd_tank  - 24./Re .* K; % from Bubbles, Drops & Particles pg [226] eq.[9-14] and K_HS from pg[225] table [9.2]. adds the effect of walls to the values measured