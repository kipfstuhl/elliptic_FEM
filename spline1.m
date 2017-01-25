
%% FEM with b-splines order 2
% use of 2nd order splines for C1 regularity of the solution

% basis
ref_spl = bspline([0 1 2]);

phi0 = ref_spl.coefs(1,:);
phi1 = ref_spl.coefs(2,:);

Dphi0 = polyder(phi0);
Dphi1 = polyder(phi1);

%% Calculations on the reference element

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

%% preparation for the assembly

% DOFs
N = M-1;                                % M+p-2=M+1-2; p-1 regularity

% cell volumes
C = x(2:end) - x(1:end-1);
% save some computations
Cinv = 1./C;

% nodes; only inner ones
nodes = x(2:end-1);                     % size(nodes) = 1, N

%% assembly of matrix
nnz = 1;

% left boundary
i=1;
ivec(nnz) = i;
jvec(nnz) = i;
Mvec(nnz) = C(i)*M00 + C(i+1)*M11;
Svec(nnz) = Cinv(i)*S00 + Cinv(i+1)*S11;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i+1;
Mvec(nnz) = C(i+1)*M01;
Svec(nnz) = Cinv(i+1)*S01;
nnz = nnz+1;

% inner nodes
for i=2:N-1
    ivec(nnz) = i;
    jvec(nnz) = i-1;
    Mvec(nnz) = C(i)*M01;
    Svec(nnz) = Cinv(i)*S01;
    nnz = nnz+1;
    
    ivec(nnz) = i;
    jvec(nnz) = i;
    Mvec(nnz) = C(i)*M00 + C(i+1)*M11;
    Svec(nnz) = Cinv(i)*S00 + Cinv(i+1)*S11;
    nnz = nnz+1;

    ivec(nnz) = i;
    jvec(nnz) = i+1;
    Mvec(nnz) = C(i+1)*M01;
    Svec(nnz) = Cinv(i+1)*S01;
    nnz = nnz+1;
end

% right boundary
i = N;
ivec(nnz) = i;
jvec(nnz) = i-1;
Mvec(nnz) = C(i)*M01;
Svec(nnz) = Cinv(i)*S01;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i;
Mvec(nnz) = C(i)*M00 + C(i+1)*M11;
Svec(nnz) = Cinv(i)*S00 + Cinv(i+1)*S11;
nnz = nnz+1;

% deallocate memory
clear Cinv

%% assembly

A = sparse(ivec, jvec, Mvec+Svec);
% display system matrix
figure
spy(A);
title('System matrix for spl1');

%% compute right-hand side
rhovec = zeros(1,N);

switch rhs_calculation
  case 'exact'
    syms xi
    % no need for distinction
    for i=1:N
        clear z

        % first part of the spline
        z = x(i) + xi*C(i);             % transform
        igrand = eval(rho)*xi;
        handle = matlabFunction(igrand(xi));
        int_left = C(i)*integral(handle, 0, 1);

        % second part of the spline
        z = x(i+1) + xi*C(i+1);
        igrand = eval(rho)*(1-xi);
        handle = matlabFunction(igrand(xi));
        int_right = C(i)*integral(handle, 0, 1);

        rhovec(i) = int_left + int_right;
    end
  case 'basis'
    z = nodes;
    Mass = sparse(ivec, jvec, Mvec);
    rhovec = (Mass * eval(rho)')';
  otherwise
    error('rhs_calculation must be ''basis'' or ''exact''');
end

%% solve the system

u = A\rhovec';
u = [0 u' 0];