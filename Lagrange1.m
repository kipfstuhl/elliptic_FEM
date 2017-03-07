clear('ivec', 'jvec', 'Mvec', 'Svec')


% Define basis functions

% basis
phi0 = [-1 1];
phi1 = [ 1 0];

% differential
Dphi0 = polyder(phi0);
Dphi1 = polyder(phi1);

%% Calculations on the reference element
% i.e. the interval [0,1]

% use convolution of vectors for polynomial multiplication

% mass matrix
ptemp = polyint(conv(phi0, phi0));
M00 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(phi1, phi0));
M10 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(phi1, phi1));
M11 = polyval(ptemp, 1) - polyval(ptemp, 0);

% stiffness matrix
ptemp = polyint(conv(Dphi0, Dphi0));
S00 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dphi1, Dphi0));
S10 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dphi1, Dphi1));
S11 = polyval(ptemp, 1) - polyval(ptemp, 0);

%% preparations for the assembly

% cell volumes
C = x(2:end) - x(1:end-1);
% dof
N = 2*M - (M+1);

nodes = x(2:end-1);

%% assembly of system matrix

nnz = 1;                                % counter for number of
                                        % non-zeros
                                        
% use i as counter
% indices for C are shifted one compared to lecture
i = 1;
ivec(nnz) = i;                          % row indices
jvec(nnz) = i;                          % column indices
Mvec(nnz) = C(i)*M11 + C(i+1)*M00;
Svec(nnz) = 1/C(i)*S11 + 1/C(i+1)*S00;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i+1;
Mvec(nnz) = C(i+1)*M10;
Svec(nnz) = 1/C(i+1)*S10;
nnz = nnz+1;

% for inner element the procedure is the same
for i=2:N-1
    ivec(nnz) = i;
    jvec(nnz) = i-1;
    Mvec(nnz) = C(i)*M10;
    Svec(nnz) = 1/C(i)*S10;
    nnz = nnz+1;
    
    ivec(nnz) = i;
    jvec(nnz) = i;
    Mvec(nnz) = C(i)*M11 + C(i+1)*M00;
    Svec(nnz) = 1/C(i)*S11 + 1/C(i+1)*S00;
    nnz = nnz+1;
    
    ivec(nnz) = i;
    jvec(nnz) = i+1;
    Mvec(nnz) = C(i+1)*M10;
    Svec(nnz) = 1/C(i+1)*S10;
    nnz = nnz+1;
    
end

i = N;
ivec(nnz) = i;
jvec(nnz) = i-1;
Mvec(nnz) = C(i)*M10;
Svec(nnz) = 1/C(i)*S10;
nnz = nnz+1;
    
ivec(nnz) = i;
jvec(nnz) = i;
Mvec(nnz) = C(i)*M11 + C(i+1)*M00;
Svec(nnz) = 1/C(i)*S11 + 1/C(i+1)*S00;
nnz = nnz+1;

% compose A
A = sparse(ivec, jvec, Mvec+Svec);

% display structure of A
figure
spy(A)
title('A')

%% Definition of right-hand side vector

switch rhs_calculation
  case 'exact'
    % calculate the integral for both basis fcts and add the
    % solutions to receive the vector entry for this element.
    rhovec = zeros(1, N);
    syms xi
    for i=1:N
        clear z
        % first basis function, this is xi
        % evaluated at the left side of the node
        z = x(i) + xi*C(i);             % transformation for left side
        igrand1(xi) = eval(rho)*xi;     % rho == f, in notes
        han1 = matlabFunction(igrand1(xi));
        int1 = C(i)*integral(han1, 0, 1);

        % second basis function, this is 1-xi
        % evaluated at the right side of the node
        z = x(i+1) + xi*C(i+1);         % transformation for right side
        igrand2(xi) = eval(rho)*(1-xi);
        han2 = matlabFunction(igrand2(xi));
        int2 = C(i+1)*integral(han2, 0, 1);

        % sum into discretisation of rho
        rhovec(i) = int1 + int2;

    end
        
  case 'basis'
    z = nodes;
    Mass = sparse(ivec, jvec, Mvec);
    rhovec = (Mass * eval(rho)')';
  otherwise
    error(['rhs_calculation has to be ''exact'' or ''basis'', you ' ...
           'used %s'], rhs_calculation);
end

%% solve the resulting system

u = A \ rhovec';
u = [0 u' 0];                           % get the boundary
                                        % conditions right
                                        % a feature of Lagrange
                                        % basis


%% Prepare output
% construct the values needed for error integration
% the integration is done using Boole's rule according to the
% exercise sheet
% this means the values at 0, 0.25, 0.5, 0.75 and 1.0 are needed

% in this case the values of the basis function are easyly seen,
% i.e. 4/4 3/4 2/4 1/4 0/4 for the points.
u_000 = u(1:end-1);
u_025 = 3.0/4.0*u(1:end-1) + 1.0/4.0*u(2:end);
u_050 = 2.0/4.0*u(1:end-1) + 2.0/4.0*u(2:end);
u_075 = 1.0/4.0*u(1:end-1) + 3.0/4.0*u(2:end);
u_100 = u(2:end);

% Differential for the H^1 error estimate
% the differential is constant over the cell interiors, because the
% basis functions are linear
% the values have to be divided by the cell volumes to get the
% right result
Du = (u(2:end) - u(1:end-1))./C;

Du_000 = Du;
Du_025 = Du;
Du_050 = Du;
Du_075 = Du;
Du_100 = Du;