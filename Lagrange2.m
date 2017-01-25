
% Define basis functions

% basis
phi0 = conv([2 -1], [ 1 -1]);
phi1 = conv([4  0], [-1  1]);
phi2 = conv([2 -1], [ 1  0]);

% differential
Dphi0 = polyder(phi0);
Dphi1 = polyder(phi1);
Dphi2 = polyder(phi2);


%% Calculations on the reference element

% mass matrix
pt = polyint(conv(phi0, phi0));
M00 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(phi0, phi1));
M01 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(phi0, phi2));
M02 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(phi1, phi1));
M11 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(phi1, phi2));
M12 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(phi2, phi2));
M22 = polyval(pt, 1) - polyval(pt, 0);

% stiffness matrix
pt = polyint(conv(Dphi0, Dphi0));
S00 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(Dphi0, Dphi1));
S01 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(Dphi0, Dphi2));
S02 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(Dphi1, Dphi1));
S11 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(Dphi1, Dphi2));
S12 = polyval(pt, 1) - polyval(pt, 0);

pt = polyint(conv(Dphi2, Dphi2));
S22 = polyval(pt, 1) - polyval(pt, 0);

%% preparations for the assembly

% dof
N = 3*M - (M+1)                         % 2M-1 => odd number

% cell volumes
C = zeros(1,N);
C(1:2:end)= x(2:end) - x(1:end-1);

% avoid expensive computations of the inverse.
% spend some memory to gain speed; this can't be optimized by the
% JIT Compiler
Cinv = zeros(1,N);                      % inverse of C
Cinv(1:2:end) = 1./C(1:2:end);          % avoid division by 0

% nodes
% first and last node is not needed. The value of the function at
% this points is known due to boundary conditions.
% This results in a node vector that has two entries less than expected
nodes = zeros(1,N);                      % even: vertex; odd: midpoint
nodes(1:2:end) = (x(2:end) + x(1:end-1))/2;
nodes(2:2:end-1) = x(2:end-1);

%% assembly of system matrix

nnz = 1;                                % counter of non-zeros

% boundary conditions
% use i as counter for convenience

i=1;                                    % odd
ivec(nnz) = i;
jvec(nnz) = i;
Mvec(nnz) = C(i)*M11;
Svec(nnz) = Cinv(i)*S11;
nnz = nnz + 1;

ivec(nnz) = i;
jvec(nnz) = i+1;
Mvec(nnz) = C(i)*M12;
Svec(nnz) = Cinv(i)*S12;
nnz = nnz+1;

i=2;                                    % even
ivec(nnz) = i;
jvec(nnz) = i-1;
Mvec(nnz) = C(i-1)*M12;
Svec(nnz) = Cinv(i-1)*S12;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i;
Mvec(nnz) = C(i-1)*M22 + C(i+1)*M00;
Svec(nnz) = Cinv(i-1)*S22 + Cinv(i+1)*S00;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i+1;
Mvec(nnz) = C(i+1)*M01;
Svec(nnz) = Cinv(i+1)*S01;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i+2;
Mvec(nnz) = C(i+1)*M02;
Svec(nnz) = Cinv(i+1)*S02;
nnz = nnz+1;


% even i
for i=4:2:N-3                           % 3 and N-2 are odd
    ivec(nnz) = i;
    jvec(nnz) = i-2;
    Mvec(nnz) = C(i-1)*M02;
    Svec(nnz) = Cinv(i-1)*S02;
    nnz = nnz+1;

    ivec(nnz) = i;
    jvec(nnz) = i-1;
    Mvec(nnz) = C(i-1)*M12;
    Svec(nnz) = Cinv(i-1)*S12;
    nnz = nnz+1;

    ivec(nnz) = i;
    jvec(nnz) = i;
    Mvec(nnz) = C(i-1)*M22 + C(i+1)*M00;
    Svec(nnz) = Cinv(i-1)*S22 + Cinv(i+1)*S00;
    nnz = nnz+1;

    ivec(nnz) = i;
    jvec(nnz) = i+1;
    Mvec(nnz) = C(i+1)*M01;
    Svec(nnz) = Cinv(i+1)*S01;
    nnz = nnz+1;

    ivec(nnz) = i;
    jvec(nnz) = i+2;
    Mvec(nnz) = C(i+1)*M02;
    Svec(nnz) = Cinv(i+1)*S02;
    nnz = nnz+1;
end

% odd i
for i=3:2:N-2
    ivec(nnz) = i;
    jvec(nnz) = i-1;
    Mvec(nnz) = C(i)*M01;
    Svec(nnz) = Cinv(i)*S01;
    nnz = nnz+1;

    ivec(nnz) = i;
    jvec(nnz) = i;
    Mvec(nnz) = C(i)*M11;
    Svec(nnz) = Cinv(i)*S11;
    nnz = nnz+1;

    ivec(nnz) = i;
    jvec(nnz) = i+1;
    Mvec(nnz) = C(i)*M12;
    Svec(nnz) = Cinv(i)*S12;
    nnz = nnz+1;
end

i = N-1;                                % even i
ivec(nnz) = i;
jvec(nnz) = i-2;
Mvec(nnz) = C(i-1)*M02;
Svec(nnz) = Cinv(i-1)*S02;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i-1;
Mvec(nnz) = C(i-1)*M12;
Svec(nnz) = Cinv(i-1)*S12;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i;
Mvec(nnz) = C(i-1)*M22 + C(i+1)*M00;
Svec(nnz) = Cinv(i-1)*S22 + Cinv(i+1)*S00;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i+1;
Mvec(nnz) = C(i+1)*M01;
Svec(nnz) = Cinv(i+1)*S01;
nnz = nnz+1;

i = N;                                  % odd i
ivec(nnz) = i;
jvec(nnz) = i-1;
Mvec(nnz) = C(i)*M01;
Svec(nnz) = Cinv(i)*S01;
nnz = nnz+1;

ivec(nnz) = i;
jvec(nnz) = i;
Mvec(nnz) = C(i)*M11;
Svec(nnz) = Cinv(i)*S11;
nnz = nnz+1;

% remove unused variables
clear Cinv

%% actual assembly of the matrix

A = sparse(ivec, jvec, Mvec+Svec);
% display system matrix
figure
spy(A)
title('System matrix for Lagrange2')


rhovec = zeros(1,N);                    % allocate memory
switch rhs_calculation
  case 'exact'
    syms xi;

    % handle evens and odds dirrferent

    % i odd
    for i=1:2:N
        clear z igrand handle int
        % (i+1)/2 is the corresponding index (i) in vector x
        z = x((i+1)/2) + xi*C(i);     % variable transform R_{i-1}
        igrand(xi) = eval(rho)*poly2sym(phi1, xi); % symbolic representation
        handle = matlabFunction(igrand(xi)); % create function handle
        int = C(i)*integral(handle, 0, 1); % integrate transformed
                                        % function
        rhovec(i) = int;
    end

    % i even
    for i=2:2:N-1
        clear z igrand handle int
        % (i)/2 is the corresponding index (i-1) in vector x
        ind = (i)/2;
        z = x(ind) + xi*C(i-1);       % variable transform R_{i-2}
        igrand(xi) = eval(rho)*poly2sym(phi2, xi);
        handle = matlabFunction(igrand(xi)); % handle for integration
        int = C(i-1)*integral(handle, 0, 1);

        rhovec(i) = int;

        clear z igrand handle int
        % (i+2)/2 is the corresponding index (i+1) in vector x
        ind = (i+2)/2;
        z = x(ind) + xi*C(i+1);       % variable transform R_{i-2}
        igrand(xi) = eval(rho)*poly2sym(phi0, xi);
        handle = matlabFunction(igrand(xi)); % handle for integration
        int = C(i+1)*integral(handle, 0, 1);

        rhovec(i) = rhovec(i) + int;
    end
        
  case 'basis'
    z = nodes;
    Mass = sparse(ivec, jvec, Mvec);
    rhovec = (Mass * eval(rho)')';
  otherwise
    error('rhs_calculation must be ''exact'' or ''basis''');
end


%% solve the system

u = A\rhovec';
u = [0 u' 0];                           % boundary conds