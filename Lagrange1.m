
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