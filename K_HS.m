function [K] = K_HS(lambda)
% Correction factor for wall effects, from Fayon & Happel (1960; summarised in Clift et al., 1978).
% Derived from Bubbles Drops & Particles pg [226] eq.[9-14] and K_HS from pg[225] table [9.2]. adds the effect of walls to the values measured
K = (1 - 0.75857 * lambda.^5) ./ (1 - 2.1050.*lambda + 2.0865.*lambda.^3 - 1.7068.*lambda.^5 + 0.72603.*lambda.^6);