clear('ivec', 'jvec', 'Mvec', 'Svec')


%% FEM with b-splines of order 3
% use splines of degree 3 for C2 regularity

%% basis
ref_spl = bspline([0 1 2 3]);
ref_left = bspline([0 0 1 2]);
ref_right = bspline([1 2 3 3]);


% strange numbering due to definition
% if the numbering is chosen like this the according calculations
% become more obvious. For example the difference of the indices
% translate directly to the assigned indices.
% Thus the spline is defined like this:
% 
%           /  phi1(x)  x in [0,1]
% spl(x) = {   phi0(x)  x in [1,2]
%           \  phi2(x)  x in [2,3]
% 
phi0 = ref_spl.coefs(1, :);
phi1 = ref_spl.coefs(2, :);
phi2 = ref_spl.coefs(3, :);

left0 = ref_left.coefs(1, :);
left2 = ref_left.coefs(2, :);
% This numbering is a bit weird due to the chosen numbers for phi_i
% For the calculations this is no problem, but in the programming
% it may seem wrong at the first glance.

right0 = ref_right.coefs(2, :);
right1 = ref_right.coefs(1, :);
% only use the left part of the spline, i.e. number 0 and 1 of the polynomials

% diffferentials

Dphi0 = polyder(phi0);
Dphi1 = polyder(phi1);
Dphi2 = polyder(phi2);

Dl0 = polyder(left0);
Dl2 = polyder(left2);

Dr0 = polyder(right0);
Dr1 = polyder(right1);


%% Calculations on the reference element

% mass matrix
ptemp = polyint(conv(phi0, phi0));
M00 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(phi0, phi1));
M01 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(phi0, phi2));
M02 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(phi1, phi1));
M11 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(phi1, phi2));
M12 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(phi2, phi2));
M22 = polyval(ptemp, 1) - polyval(ptemp, 0);

% stiffness matrix
ptemp = polyint(conv(Dphi0, Dphi0));
S00 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dphi0, Dphi1));
S01 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dphi0, Dphi2));
S02 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dphi1, Dphi1));
S11 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dphi1, Dphi2));
S12 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dphi2, Dphi2));
S22 = polyval(ptemp, 1) - polyval(ptemp, 0);


%% calculations for the boundary references

% Mass matrix

% left
ptemp = polyint(conv(left0, left0));
l_M00 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(left0, phi1));
l_M01 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(left2, phi0));
l_M02 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(phi1, left2));
l_M12 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(left2, left2));
l_M22 = polyval(ptemp, 1) - polyval(ptemp, 0);

% right
ptemp = polyint(conv(right0, right0));
r_M00 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(right1, phi0));
r_M01 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(right0, phi2));
r_M02 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(right1, right1));
r_M11 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(right1, phi2));
r_M12 = polyval(ptemp, 1) - polyval(ptemp, 0);


% now for the Stiffness matrices as well

% left
ptemp = polyint(conv(Dl0, Dl0));
l_S00 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dl0, Dphi1));
l_S01 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dl2, Dphi0));
l_S02 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dl2, Dphi1));
l_S12 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dl2, Dl2));
l_S22 = polyval(ptemp, 1) - polyval(ptemp, 0);

% right
ptemp = polyint(conv(Dr0, Dr0));
r_S00 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dr1, Dphi0));
r_S01 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dr0, Dphi2));
r_S02 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dr1, Dr1));
r_S11 = polyval(ptemp, 1) - polyval(ptemp, 0);

ptemp = polyint(conv(Dr1, Dphi2));
r_S12 = polyval(ptemp, 1) - polyval(ptemp, 0);


%% preparation for the assembly

% DOF
N = M;   % M-p + 2p-2 = M+p-2

% Cell volumes
C = x(2:end) - x(1:end-1);
% inverse of Cell volumes
Cinv = 1./C;
% this does not improve much, because most time is spent in the
% evaluation of symbolic expressions.

% node Vector
nodes = x(2:end-1);


%% assembly of the matrix

% don't use nnz(number of non-zeros) as this is a matlab function
nz = 1;

% left boundary
i = 1;
ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i)*l_M00 + C(i+1)*l_M22;
Svec(nz) = Cinv(i)*l_S00 + Cinv(i+1)*l_M22;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+1;
Mvec(nz) = C(i)*l_M01 + C(i+1)*l_M02;
Svec(nz) = Cinv(i)*l_S01 + Cinv(i+1)*l_S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+2;
Mvec(nz) = C(i+1)*l_M12;
Svec(nz) = Cinv(i+1)*l_S12;
nz = nz+1;

i = 2;
ivec(nz) = i;
jvec(nz) = i-1;
Mvec(nz) = C(i-1)*l_M01 + C(i)*l_M02;
Svec(nz) = Cinv(i-1)*l_S01 + Cinv(i)*l_S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i-1)*M11 + C(i)*M00 + C(i+1)*M22;
Svec(nz) = Cinv(i-1)*S11 + Cinv(i)*S00 + Cinv(i+1)*S22;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+1;
Mvec(nz) = C(i)*M01 + C(i+1)*M02;
Svec(nz) = Cinv(i)*S01 + Cinv(i+1)*S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+2;
Mvec(nz) = C(i+1)*M12;
Svec(nz) = Cinv(i+1)*S12;
nz = nz+1;

i = 3;
ivec(nz) = i;
jvec(nz) = i-2;
Mvec(nz) = C(i-1)*l_M12;
Svec(nz) = Cinv(i-1)*l_S12;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i-1;
Mvec(nz) = C(i-1)*M01 + C(i)*M02;
Svec(nz) = Cinv(i-1)*S01 + Cinv(i)*S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i-1)*M11 + C(i)*M00 + C(i+1)*M22;
Svec(nz) = Cinv(i-1)*S11 + Cinv(i)*S00 + Cinv(i+1)*S22;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+1;
Mvec(nz) = C(i)*M01 + C(i+1)*M02;
Svec(nz) = Cinv(i)*S01 + Cinv(i+1)*S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+2;
Mvec(nz) = C(i+1)*M12;
Svec(nz) = Cinv(i+1)*S12;
nz = nz+1;

% now the inner ones
for i=4:N-3
    ivec(nz) = i;
    jvec(nz) = i-2;
    Mvec(nz) = C(i-1)*M12;
    Svec(nz) = Cinv(i-1)*S12;
    nz = nz+1;
    
    ivec(nz) = i;
    jvec(nz) = i-1;
    Mvec(nz) = C(i-1)*M01 + C(i)*M02;
    Svec(nz) = Cinv(i-1)*S01 + Cinv(i)*S02;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i;
    Mvec(nz) = C(i-1)*M11 + C(i)*M00 + C(i+1)*M22;
    Svec(nz) = Cinv(i-1)*S11 + Cinv(i)*S00 + Cinv(i+1)*S22;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i+1;
    Mvec(nz) = C(i)*M01 + C(i+1)*M02;
    Svec(nz) = Cinv(i)*S01 + Cinv(i+1)*S02;
    nz = nz+1;

    ivec(nz) = i;
    jvec(nz) = i+2;
    Mvec(nz) = C(i+1)*M12;
    Svec(nz) = Cinv(i+1)*S12;
    nz = nz+1;
end

% right boundary
i=N-2;
ivec(nz) = i;
jvec(nz) = i-2;
Mvec(nz) = C(i-1)*M12;
Svec(nz) = Cinv(i-1)*S12;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i-1;
Mvec(nz) = C(i-1)*M01 + C(i)*M02;
Svec(nz) = Cinv(i-1)*S01 + Cinv(i)*S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i-1)*M11 + C(i)*M00 + C(i+1)*M22;
Svec(nz) = Cinv(i-1)*S11 + Cinv(i)*S00 + Cinv(i+1)*S22;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+1;
Mvec(nz) = C(i)*M01 + C(i+1)*M02;
Svec(nz) = Cinv(i)*S01 + Cinv(i+1)*S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+2;
Mvec(nz) = C(i+1)*r_M12;
Svec(nz) = Cinv(i+1)*r_S12;
nz = nz+1;

i=N-1;
ivec(nz) = i;
jvec(nz) = i-2;
Mvec(nz) = C(i-1)*M12;
Svec(nz) = Cinv(i-1)*S12;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i-1;
Mvec(nz) = C(i-1)*M01 + C(i)*M02;
Svec(nz) = Cinv(i-1)*S01 + Cinv(i)*S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i-1)*M11 + C(i)*M00 + C(i+1)*M22;
Svec(nz) = Cinv(i-1)*S11 + Cinv(i)*S00 + Cinv(i+1)*S22;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i+1;
Mvec(nz) = C(i)*r_M01 + C(i+1)*r_M02;
Svec(nz) = Cinv(i)*r_S01 + Cinv(i+1)*r_S02;
nz = nz+1;

i = N;
ivec(nz) = i;
jvec(nz) = i-2;
Mvec(nz) = C(i-1)*r_M12;
Svec(nz) = Cinv(i-1)*r_S12;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i-1;
Mvec(nz) = C(i-1)*r_M01 + C(i)*r_M02;
Svec(nz) = Cinv(i-1)*r_S01 + Cinv(i)*r_S02;
nz = nz+1;

ivec(nz) = i;
jvec(nz) = i;
Mvec(nz) = C(i-1)*r_M11 + C(i)*r_M00;
Svec(nz) = Cinv(i-1)*r_S11 + Cinv(i)*r_S00;
nz = nz+1;


% remove unused variables
clear Cinv

%% actual assembly of the matrix

A = sparse(ivec, jvec, Mvec+Svec);
% display system matrix for debugging purposes
figure
spy(A)
title('System matrix for quadratic splines')

%% right hand side

rhovec = zeros(1,N);                    % allocate memory
switch rhs_calculation
  case 'exact'
    syms xi
    
    % left boundary
    i = 1;
    clear z igrand han int
    z = x(i) + xi*C(i);
    igrand(xi) = eval(rho)*poly2sym(left0, xi);
    han = matlabFunction(igrand(xi));
    int = C(i)*integral(han, 0, 1);
    tmp = int;

    clear z igrand han int
    z = x(i+1) + xi*C(i+1);
    igrand(xi) = eval(rho)*poly2sym(left2, xi);
    han = matlabFunction(igrand(xi));
    int = C(i+1)*integral(han, 0, 1);
    tmp = tmp + int;

    rhovec(i) = tmp;
    tmp = 0.0;

    % inner nodes
    for i=2:N-1
        clear z igrand han int
        z = x(i-1) + xi*C(i-1);
        igrand(xi) = eval(rho)*poly2sym(phi1, xi);
        han = matlabFunction(igrand(xi));
        int = C(i-1)*integral(han, 0, 1);
        tmp = int;
        
        clear z igrand han int
        z = x(i) + xi*C(i);
        igrand(xi) = eval(rho)*poly2sym(phi0, xi);
        han = matlabFunction(igrand(xi));
        int = C(i)*integral(han, 0, 1);
        tmp = tmp + int;
        
        clear z igrand han int
        z = x(i+1) + xi*C(i+1);
        igrand(xi) = eval(rho)*poly2sym(phi2, xi);
        han = matlabFunction(igrand(xi));
        int = C(i+1)*integral(han, 0, 1);
        tmp = tmp + int;

        rhovec(i) = tmp;
        tmp = 0.0;
    end

    % right boundary
    i = N;
    clear z igrand han int
    z = x(i-1) + xi*C(i-1);
    igrand(xi) = eval(rho)*poly2sym(right1, xi);
    han = matlabFunction(igrand(xi));
    int = C(i-1)*integral(han, 0, 1);
    tmp = int;

    clear z igrand han int
    z = x(i) + xi*C(i);
    igrand(xi) = eval(rho)*poly2sym(right0, xi);
    han = matlabFunction(igrand(xi));
    int = C(i)*integral(han, 0, 1);
    tmp = tmp + int;
    
    rhovec(i) = tmp;


  case 'basis'
    %z = nodes;

    %double check this
    z = (x(1:end-1) + x(2:end))./2.0;
    
    Mass = sparse(ivec, jvec, Mvec);
    rhovec = (Mass * eval(rho)')';
    
end


%% solve system

u = A\rhovec';