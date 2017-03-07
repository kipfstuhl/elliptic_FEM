clear('ivec', 'jvec', 'Mvec', 'Svec')


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

% this "optimisation" improves the runtime by approximately
% nothing, because most of the time (over 90%) is spent evaluating
% the symbolic expressions.



% nodes
% first and last node is not needed. The value of the function at
% this points is known due to boundary conditions.
% This results in a node vector that has two entries less than expected
nodes = zeros(1,N);                      % even: vertex; odd: midpoint
nodes(1:2:end) = (x(2:end) + x(1:end-1))/2;
nodes(2:2:end-1) = x(2:end-1);

%% assembly of system matrix

nz = 1;                                % counter of non-zeros

% boundary conditions
% use i as counter for convenience

i=1;                                    % odd
ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i)*M11;
Svec(nz) = Cinv(i)*S11;
nz = nz + 1;

ivec(nz) = i;
jvec(nz) = i+1;
Mvec(nz) = C(i)*M12;
Svec(nz) = Cinv(i)*S12;
nz = nz+1;

i=2;                                    % even
ivec(nz) = i;
jvec(nz) = i-1;
Mvec(nz) = C(i-1)*M12;
Svec(nz) = Cinv(i-1)*S12;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i-1)*M22 + C(i+1)*M00;
Svec(nz) = Cinv(i-1)*S22 + Cinv(i+1)*S00;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+1;
Mvec(nz) = C(i+1)*M01;
Svec(nz) = Cinv(i+1)*S01;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+2;
Mvec(nz) = C(i+1)*M02;
Svec(nz) = Cinv(i+1)*S02;
nz = nz+1;


% even i
for i=4:2:N-3                           % 3 and N-2 are odd
    ivec(nz) = i;
    jvec(nz) = i-2;
    Mvec(nz) = C(i-1)*M02;
    Svec(nz) = Cinv(i-1)*S02;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i-1;
    Mvec(nz) = C(i-1)*M12;
    Svec(nz) = Cinv(i-1)*S12;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i;
    Mvec(nz) = C(i-1)*M22 + C(i+1)*M00;
    Svec(nz) = Cinv(i-1)*S22 + Cinv(i+1)*S00;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i+1;
    Mvec(nz) = C(i+1)*M01;
    Svec(nz) = Cinv(i+1)*S01;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i+2;
    Mvec(nz) = C(i+1)*M02;
    Svec(nz) = Cinv(i+1)*S02;
    nz = nz+1;
end

% odd i
for i=3:2:N-2
    ivec(nz) = i;
    jvec(nz) = i-1;
    Mvec(nz) = C(i)*M01;
    Svec(nz) = Cinv(i)*S01;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i;
    Mvec(nz) = C(i)*M11;
    Svec(nz) = Cinv(i)*S11;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i+1;
    Mvec(nz) = C(i)*M12;
    Svec(nz) = Cinv(i)*S12;
    nz = nz+1;
end

i = N-1;                                % even i
ivec(nz) = i;
jvec(nz) = i-2;
Mvec(nz) = C(i-1)*M02;
Svec(nz) = Cinv(i-1)*S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i-1;
Mvec(nz) = C(i-1)*M12;
Svec(nz) = Cinv(i-1)*S12;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i-1)*M22 + C(i+1)*M00;
Svec(nz) = Cinv(i-1)*S22 + Cinv(i+1)*S00;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+1;
Mvec(nz) = C(i+1)*M01;
Svec(nz) = Cinv(i+1)*S01;
nz = nz+1;

i = N;                                  % odd i
ivec(nz) = i;
jvec(nz) = i-1;
Mvec(nz) = C(i)*M01;
Svec(nz) = Cinv(i)*S01;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i)*M11;
Svec(nz) = Cinv(i)*S11;
nz = nz+1;

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


%% Prepare Output
% construct the values needed for error integration
% the integration is done using Boole's rule according to the
% exercise sheet
% this means the values at 0, 0.25, 0.5, 0.75 and 1.0 are needed

% exploit the fact, that the points 0.0, 0.5 and 1.0 are the
% interpolation points for the Lagrange polynomials, this means
% only exactly one basis function has value 1.0, the others have
% value 0.0

u_000 = u(1:2:end-2);

pv0 = polyval(phi0, 0.25);
pv1 = polyval(phi1, 0.25);
pv2 = polyval(phi2, 0.25);
u_025 = pv0*u(1:2:end-2) + pv1*u(2:2:end-1) + pv2*u(3:2:end);

u_050 = pv0*u(2:2:end-1);

pv0 = polyval(phi0, 0.75);
pv1 = polyval(phi1, 0.75);
pv2 = polyval(phi2, 0.75);
u_075 = pv0*u(1:2:end-2) + pv1*u(2:2:end-1) + pv2*u(3:2:end);

u_100 = u(3:2:end);


% computations for the differentials
% use that at 0.25, 0.5 and 0.75 one basis fct has zero derivative.
% note that the values have to be divided by the cell volume.

Dv0 = polyval(Dphi0, 0.00);
Dv1 = polyval(Dphi1, 0.00);
Dv2 = polyval(Dphi2, 0.00);
Du_000 = 1./C(1:2:end).*(Dv0*u(1:2:end-2) + Dv1*u(2:2:end-1) + ...
                         Dv2*u(3:2:end));

Dv0 = polyval(Dphi0, 0.25);
Dv1 = polyval(Dphi1, 0.25);
Dv2 = polyval(Dphi2, 0.25);
Du_025 = 1./C(1:2:end).*(Dv0*u(1:2:end-2) + Dv1*u(2:2:end-1));

Dv0 = polyval(Dphi0, 0.50);
Dv1 = polyval(Dphi1, 0.50);
Dv2 = polyval(Dphi2, 0.50);
Du_050 = 1./C(1:2:end).*(Dv0*u(1:2:end-2) + Dv2*u(3:2:end));

Dv0 = polyval(Dphi0, 0.75);
Dv1 = polyval(Dphi1, 0.75);
Dv2 = polyval(Dphi2, 0.75);
Du_075 = 1./C(1:2:end).*(Dv1*u(2:2:end-1) + Dv2*u(3:2:end));

Dv0 = polyval(Dphi0, 1.00);
Dv1 = polyval(Dphi1, 1.00);
Dv2 = polyval(Dphi2, 1.00);
Du_100 = 1./C(1:2:end).*(Dv0*u(1:2:end-2) + Dv1*u(2:2:end-1) + ...
                         Dv2*u(3:2:end));

