% =================================================================
% Target equation: the Lorenz 96 equation
%       udot_{i} = ( u_{i+1} - u_{i-2} ) * u_{i-1} - u_{i} + F
%
%
% Auxiliary functions:
%   - Dictionary.m (The Candidate Functions)
%   - DouglasRachford.m (Optimization Routine) 
%   - leg2mon.m (Legendre to Monomial Transform)
%   - lorenz96.m (ODE to Test)
%
% Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
% Date: May 9, 2018
% =================================================================

clc; clear all; close all;

% =================================================================
% Steps 1-3: Construct the data matrix U and the velocity vector V.
% =================================================================

% parameters
n = 128; % number of equations in the ODE system
F = 8; % constant in the ODE system
dt = 0.001; % step size
T = dt; % terminal time -- one time step
K = 2; % number of bursts
var = 0.0001; % noise level
nb = 55; % size of the block (required: nb odd, nb+2 <= n, and nb+ii < n)
ii = 1; % index of the first point of the block in the entire data, corner
%parameters for Step 7
tau = 1; mu = 1/2; MaxIt = 1e5; tol = 1e-7; %Optimization Parameters


% indices used in cyclic permutation and restriction of the data
r = (nb-1)/2;
indc = 1:n; indc = [ indc(1:r+1) indc(n-r+1:n)];
indc = sort(indc); % column indices
indr = 1:n; indr = indr(ii:ii+nb-1); % row indices
nc = length(indc);
nr = length(indr);

% Example: n=5, nb=3, ii=1, 
%           data matrix U = | u1 u2 u3 u4 u5 |
%           velocity vector V = | v1 v2 v3 v4 v5 |
%   cyclic permutation of the data =>
%               | u1 u2 u3 u4 u5 |
%               | u2 u3 u4 u5 u1 |
%           U = | u3 u4 u5 u1 u2 |
%               | u4 u5 u1 u2 u3 |
%               | u5 u1 u2 u3 u4 |
%   restriction of the data =>
%               | u1 u2 u5 |
%           U = | u2 u3 u1 |
%           	| u3 u4 u2 |
%   corresponding velocity vector =>
%               | v1 |
%           V = | v2 |
%           	| v3 |

% initialization
U = zeros(nr*K,nc); % data matrix
Udot = zeros(nr*K,1); % velocity vector -- "exact" velocity for testing
V = zeros(nr*K,1); % velocity vector -- computed velocity

% multiple burst
for k=1:K
    
    % Generate the initial data.
    u0 = 2*rand(n,1)-1;
    
    % Solve for u and udot on a finer grid.
    [t,u] = ode45(@(t,x) lorenz96(t,x,F), 0:T/2:T, u0);
    % Compute the "exact" derivative.
    udot = lorenz96(0,u0,F);  
    % Record only u(t0) and u(t1).
    u = u([1,end],:);
    
    % Add noise to data.
    u = u + var*randn(2,n);
    % Compute the numerical derivative, may be unstable if noise is large.
    udot1 = ( u(2,:) - u(1,:) )/dt;
    
    % Perform cyclic permutation on U.
    v = zeros(n,n); % temporary variable
    v(1,:) = u(1,:);
    for j=2:n
        v(j,:) = circshift(v(j-1,:),[0 -1]);
    end
    % Restrict the data onto the block.
    U((k-1)*nr+1 : k*nr, :) = v(indr,indc);
    % Construct Udot and V accordingly.
    Udot((k-1)*nr+1 : k*nr, :) = udot(indr);
    V((k-1)*nr+1 : k*nr, :) = udot1(1,indr)';
    
end

% =================================================================
% Step 4:  Transform U to be a*U + b so that each element in U
%           is valued in [-1, 1].
% =================================================================

% Notations: U = unscaled data
%            U1 = scaled data

% % Scale the data to be valued in [-1, 1].
% a = 2/( max(U(:)) - min(U(:)) );
% b = -2*min(U(:))/( max(U(:)) - min(U(:)) )-1;
% U1 = a*U + b;

% Alternatively, since U is close to be uniformly distributed in 
%   [-1,1] for this test, we can set a = 1 and b = 0.
a = 1; b = 0; U1 = U; %PARAMETERS

% =================================================================
% Step 5: Construct the dictionary matrix Aleg using U1.
% =================================================================

% Notations: Amon = dictionary of monomials of U1
%            Aleg = dictionary of Legendre polynomials of U1

p = 3; % degree of the basis element (p=2 or p=3)
r = 10; % localization of the dictionary (radius of the restricted subset)
% Create dictionaries, Amon and Aleg, of monomials and Legendre
%   polynomials, respectively, of degree at most p.
[Amon,Aleg,L,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111] = Dictionary(U1,p,r,indc);
N = size(Aleg,2);

% =================================================================
% Step 6: Normalize each column of Aleg to have unit l2 norm.
% =================================================================

% Notations: Aleg = dictionary of Legendre basis (from Step 5)
%            Aleg1 = scaled Alge, l2-norm of each column is 1

% Compute the norm of each column.
Aleg_cnorm = sqrt(sum(Aleg.^2,1)); 
Aleg_cnorm = Aleg_cnorm(:);
% Normalize columns of Alge.
Aleg1 = normc(Aleg);

% =================================================================
% Step 7: Apply the Douglas-Rachford algorithm to solve
%           min_cleg        || cleg ||_1
%           subject to      || Aleg1 * cleg - V ||_2 <= sigma
% =================================================================

sigma = 1.01 * norm(Udot-V,2); % For testing purposes, in practice must be determined.
%tau = 1; mu = 1/2; MaxIt = 1e5; tol = 1e-6; %Optimization Parameters
cleg = DouglasRachford(Aleg1,V,sigma,tau,mu,MaxIt,tol);

% =================================================================
% Step 8: Map the solution cleg to the coefficients wrt the standard
%           monomial basis on U1.
% =================================================================

% Notations: cleg = coefficient for Aleg i.e. Aleg(U1) * cleg = V
%            cmon = coefficient for Amon i.e. Amon(U1) * cmon = V

% Convert cleg to be the coefficient for Aleg.
cleg = cleg./Aleg_cnorm;
% Convert cleg to be the coefficient wrt the monomial basis.
cmon = leg2mon(cleg,p,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111);
%Optional Thresholding, See Proposition 3.5
cmon = cmon.*(abs(cmon)>1e-1);
% =================================================================
% Step 9: Create the exact coefficient vector ctrue.
%                   Amon(U1) * ctrue = Udot
% =================================================================

%   For ease of comparison, ctrue is the coefficient vector wrt the
% the dictionary Amon(U1) from Step 5, which is constructed using 
% monomials of the scaled data U1.
%   When we set a=1 and b=0, ctrue is then the coefficient 
% vector wrt the the dictionary Amon(U), which is constructed using 
% monomials of the unscaled data U.

ctrue = zeros(N,1);

% Find the non-zeros terms.
str_u1 = strcat('u', num2str(1)); % u_{1}
str_u2 = strcat('u', num2str(2)); % u_{2}
str_un = strcat('u', num2str(n)); % u_{n}
str_unm = strcat('u', num2str(n-1)); % u_{n-1}
str_u2un = strcat(str_u2,str_un); % u_{2} * u_{n}
str_unmun = strcat(str_unm,str_un); % u_{n-1} * u_{1,n}

% Find the support set.
supp = [1,                          find(ismember(L,str_u1)), ...
        find(ismember(L,str_u2)),   find(ismember(L,str_unm)),  ...
        find(ismember(L,str_u2un)), find(ismember(L,str_unmun))];

% Assign values.
ctrue(supp(1)) = F + b/a; % 1
ctrue(supp(2)) = -1/a; % u_{1}
ctrue(supp(3)) = -b/a/a; % u_{2}
ctrue(supp(4)) = b/a/a; % u_{n-1}
ctrue(supp(5)) = 1/a/a; % u_{2} * u_{n}
ctrue(supp(6)) = -1/a/a; % u_{n-1} * u_{1,n}

% =================================================================
% Step 10: Compare the results.
% =================================================================

% Display results for c.
disp('sparsity of the solution');
display(nnz(cmon));
disp('solution on the correct support set (L: exact soln, R: approx soln)');
display([ctrue(supp) cmon(supp)]);
disp('relative error');
disp(norm(cmon-ctrue,2)/norm(ctrue,2));
disp('least squares relative error');
disp(norm(Amon\V-ctrue,2)/norm(ctrue,2));