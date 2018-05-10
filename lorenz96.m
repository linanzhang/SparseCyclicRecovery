function udot = lorenz96(t,u,F)

% ============================================================
% Inputs:
%   t = time
%   u = state of the ODE system
%   F = forcing term (constant)
%
% Outputs:
%	udot_{i} = ( u_{i+1} - u_{i-2} ) * u_{i-1} - u_{i} + F
%
% Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
% Date: May 9, 2018
% ============================================================

udot = (circshift(u,-1)-circshift(u,2)) ...
    .* circshift(u,1) - u + F;

end
