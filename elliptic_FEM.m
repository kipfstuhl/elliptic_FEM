
%clear all
clear('i', 'j')
close all

a = 0;                                  % left boundary
b = 2;                                  % right boundary
M = 16;                                 % number of elements

L = b - a;                              % length of interval

% 
% define data, i.e. function and rhs

% copied from reference
syms z;
u_ex(z) = -sin(2*pi/L*z)*exp(z);        % exact solution
u_ex(z) = (z-a)*(z-b)*(z-(a+b)/2)^4;
Du_ex(z) = diff(u_ex(z), 1);            % differential
rho(z) = -diff(u_ex(z), 2) + u_ex(z)    % rhs


% set solver properties

if ~exist('method')
    method = 'Lag1'
end
rhs_calculation = 'basis'
% rhs_calculation = 'exact'

%% define the grid

% type of grid calculation, numeric or analytical
grid_calculation = 'num'
%grid_calculation = 'ana'

%syms z;
%H(z) = sin(2*pi*z) + 1.2;

% for iteration and saving the errors
if ~exist('gridSpace')
    gridSpace = '0.1'
end
[x, M] = calculate_grid(a, b, M, grid_calculation, gridSpace);
figure
plot(x, zeros(size(x)), 'bx');
xlim([a b])
title('Grid from function');

switch method
  case 'Lag1'
    Lagrange1
  case 'Lag2'
    Lagrange2
  case 'spl1'
    spline1
  case 'spl2'
    spline2
  otherwise
    error(['Error!\nMethod %s not supported. Use Lag1, Lag2, spl1, ' ...
           'or spl2 instead'], method);
end


%% Error estimation
% the error is estimated via integrating the difference of the
% output u_xxx and u_ex, that is available in this case.

% prepare the x values and also the interval lengths
int_len = x(2:end) - x(1:end-1);
% parametrize the interval [x_1, x_2] by
% (1-lambda)*x_1 + lambda*x_2
x_25 = 0.75*x(1:end-1) + 0.25*x(2:end);
x_50 = 0.5*x(1:end-1) + 0.5*x(2:end);
x_75 = 0.25*x(1:end-1) + 0.75*x(2:end);

% first the L2 error
L2 = 0.0;
for i=1:M
    % assign z for symbolic evaluation
    z=x(i);
    e_000 = eval(u_ex) - u_000(i);
    z = x_25(i);
    e_025 = eval(u_ex) - u_025(i);
    z = x_50(i);
    e_050 = eval(u_ex) - u_050(i);
    z = x_75(i);
    e_075 = eval(u_ex) - u_075(i);
    z = x(i+1);
    e_100 = eval(u_ex) - u_100(i);

    % use Boole's rule for integration
    % squares needed to get a norm out of this
    L2 = L2 + ...
         int_len(i)/90*(7*e_000^2 + 32*e_025^2 + 12*e_050^2 + 32*e_075^2 + 7*e_100^2);
end
% so far the scalar product was calculated.
L2 = sqrt(L2);

% now the H1 seminorm to add to L2 norm for obtaining the H1 error
H1_semi = 0.0;
for i=1:M
    % assign z for symbolic evaluation
    z=x(i);
    e_000 = eval(Du_ex) - Du_000(i);
    z = x_25(i);
    e_025 = eval(Du_ex) - Du_025(i);
    z = x_50(i);
    e_050 = eval(Du_ex) - Du_050(i);
    z = x_75(i);
    e_075 = eval(Du_ex) - Du_075(i);
    z = x(i+1);
    e_100 = eval(Du_ex) - Du_100(i);

    % use Boole's rule for integration
    % squares needed to get a norm out of this
    H1_semi = H1_semi + ...
         int_len(i)/90*(7*e_000^2 + 32*e_025^2 + 12*e_050^2 + 32*e_075^2 + 7*e_100^2);
end

H1 = sqrt(L2^2 + H1_semi);