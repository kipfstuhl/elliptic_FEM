
clear all
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
Du_ex(z) = diff(u_ex(z), 1);            % differential
rho(z) = -diff(u_ex(z), 2) + u_ex(z)    % rhs


% set solver properties
method = 'spl1'
rhs_calculation = 'exact'
% rhs_calculation = 'basis'

%% define the grid

% type of grid calculation, numeric or analytical
grid_calculation = 'num'
%grid_calculation = 'ana'

%syms z;
%H(z) = sin(2*pi*z) + 1.2;
[x, M] = calculate_grid(a, b, M, grid_calculation, '0.1');
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
    Spline2
  otherwise
    error(['Error!\nMethod %s not supported. Use Lag1, Lag2, spl1, ' ...
           'or spl2 instead'], method);
end

