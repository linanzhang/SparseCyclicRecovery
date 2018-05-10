function [Dmon,Dleg,L,Ind1,Ind20,Ind11,Ind300,Ind210,Ind120,Ind111] = Dictionary(u,p,r,indc)

% ============================================================
% Description:
%   The function outputs the Legendre Dictionary Matrix.
%
% Inputs:
%   u = the state of a nonlinear system, i.e. given data
%   p = maximum degree of polynomials in the dictionary (p<=3)
%   r = number of left/right neighbors
%   indc = column index (restriction of the data)
%
% Outputs:
%   Dmon = dictionary of monomials of degree at most p
%   Dleg = dictionary of Legendre polynomials of degree at most p
%   L = legend of D
%   Ind1 = indices for terms of the form ui
%   Ind20 = indices for terms of the form (ui)^2
%   Ind11 = indices for terms of the form (ui)*(uj)
%   Ind300 = indices for terms of the form (ui)^3
%   Ind210 = indices for terms of the form (ui)^2*(uj)
%   Ind120 = indices for terms of the form (ui)*(uj)^2
%   Ind111 = indices for terms of the form (ui)*(uj)*(uk)
%
% Example: u=[u1,u2,u3,u8,u9], p=2, r=1, indc=[1,2,3,8,9]
%   Dmon = [1, u1, u2, u9, u1^2, u1u2, u1u9, u2^2, u2u9, u9^2]
%   Dleg = [1, sqrt(3)*u1, sqrt(3)*u2, sqrt(3)*u9, 
%           (3*u1^2-1)*sqrt(5)/2, 3*u1u2, 3*u1u9, 
%           (3*u2^2-1)*sqrt(5)/2, 3*u2u9, (3*u9^2-1)*sqrt(5)/2]
%
% Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
% Date: May 9, 2018
% ============================================================

% Compute the size of u.
%   m = number of time grid
%   n = number of states.
[m,n] = size(u);
% List the name of u1.
l1 = cell(1,n+1);
for ii=1:n
    l1{ii+1} = strcat('u',num2str(indc(ii)));
end
% Augment u to [1;u] and down-sample the columns.
if n>= 2*r+1
    umon = [ones(m,1) u(:,1:r+1) u(:,n-r+1:n)];
    uleg = [ones(m,1) u(:,1:r+1)*sqrt(3) u(:,n-r+1:n)*sqrt(3)];
    l1 = l1([1:r+2,n-r+2:n+1]);
else
    umon = [ones(m,1) u];
    uleg = [ones(m,1) u*sqrt(3)];
end

% Redefine n.
n = min(2*r+1,n);

% Find all polynomials of degree at most 3.
C = 1:(n+1);
for ii=1:p-1
    C = combvec(C,1:(n+1));
end
% Sort each column of C in ascending order.
C = sort(C,1);
% Remove duplicate rows in C' (i.e. duplicate columns in C).
C = unique(C','rows');
% number of rows of C = number of columns of D
% Note: number of columns of C = 3
nD = size(C,1);
% Construct D.
Dmon = ones(m,nD);
Dleg = ones(m,nD);
L = cell(nD,1);

%   P2: 3*ui*uj for i~=j
%   P3: sqrt(27)*ui*uj*uk
for ii=1:nD
    % Multiply the corresponding columns to make a polynomial of degree at most p.
    for jj = 1:p
        Dmon(:,ii) = Dmon(:,ii) .* umon(:,C(ii,jj));
        Dleg(:,ii) = Dleg(:,ii) .* uleg(:,C(ii,jj));
        L{ii} = strcat(L{ii},l1{C(ii,jj)});
    end
end
L{1} = num2str(1);

if p==2
    
    %   P2: (3*ui^2-1)*sqrt(5)/2
    Ind20 = [];
    for ii=2:n+1
        ind = find(ismember(C,[1 1]*ii,'rows'));
        Dleg(:,ind) = (Dleg(:,ind)-1)*sqrt(5)/2;
        Ind20 = [Ind20 ind];
    end
    Ind1 = 2:n+1;
    Ind11 = 1:nD;
    Ind11([1 Ind1 Ind20]) = [];
    %   P3
    Ind300 = []; Ind210 = []; Ind120 = []; Ind111 = [];
    
elseif p==3
    
    %   P2: (3*ui^2-1)*sqrt(5)/2
    N = n^2/2+3*n/2+1;
    Ind20 = [];
    for ii=2:n+1
        ind = find(ismember(C(1:N,2:3),[1 1]*ii,'rows'));
        Dleg(:,ind) = (Dleg(:,ind)-1)*sqrt(5)/2;
        Ind20 = [Ind20 ind];
    end
    Ind1 = 2:n+1;
    Ind11 = 1:N; Ind11([1 Ind1 Ind20]) = [];
    
    %   P3: sqrt(15)*(3*ui^2-1)*uj for i~=j
    Ind210 = [];
    for ii=2:n+1
        ind = find(ismember(C(N+1:end,1:2),[1 1]*ii,'rows'));
        for jj=1:length(ind)
            Dleg(:,ind(jj)+N) = ( 3*u(:,ii-1).^2 - 1 ) .* uleg(:,C(ind(jj)+N,3)) * sqrt(5)/2;
        end
        Ind210 = [Ind210; ind];
    end
    Ind210 = sort(Ind210); Ind210 = unique(Ind210') + N;
    Ind120 = [];
    for ii=2:n+1
        ind = find(ismember(C(N+1:end,2:3),[1 1]*ii,'rows'));
        for jj=1:length(ind)
            Dleg(:,ind(jj)+N) = ( 3*u(:,ii-1).^2 - 1 ) .* uleg(:,C(ind(jj)+N,1)) * sqrt(5)/2;
        end
        Ind120 = [Ind120; ind];
    end
    Ind120 = sort(Ind120); Ind120 = unique(Ind120') + N;
    
    %   P3: (5*ui^3-3*ui)*sqrt(7)/2
    Ind300 = [];
    for ii=2:n+1
        ind = find(ismember(C,[1 1 1]*ii,'rows'));
        Dleg(:,ind) = ( 5*u(:,ii-1).^3 - 3*u(:,ii-1) ) * sqrt(7)/2;
        Ind300 = [Ind300 ind];
    end
    Ind210(find(ismember(Ind210,Ind300))) = [];
    Ind120(find(ismember(Ind120,Ind300))) = [];
    
    % P3: sqrt(27)*ui*uj*uk
    Ind111 = 1:nD; 
    Ind111([1 Ind1 Ind20 Ind11 Ind300 Ind210 Ind120]) = [];
end
end