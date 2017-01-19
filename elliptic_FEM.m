
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
method = 'Lag1'
rhs = 'exact'

%% define the grid

% type of grid calculation, numeric or analytical
grid_calculation = 'num'
% grid_calculation = 'ana'


switch grid_calculation
  case 'ana'
    syms z
    H(z) = sin(2*pi*z) + 1.2;           % spacing function
    c = M/int(1/H, a, b);               % normalisaton constant
    g(z) = c * int(1/H);
    g(z) = g(z) - g(a);                 % set integration constant
    x = zeros(1,M+1);                   % allocate memory
    x(1) = a;                           % set first vertex
    for i=1:M-1
        vertex = solve(g(z) - i == 0, z);
        x(i+1) = vpa(vertex);
    end
    x(M+1) = b;

    % visualisation of the grid
    z = linspace(a,b,100);              % symbolic not needed
                                        % anymore
    plot(x, zeros(size(x)), 'bx', z, vpa(H(z)));
    xlim([a b]);
    
    
  case 'num'
    H = '0.1'
    count = 0;                          % number of vertices
    yvec = a;                           % first vertex
    y = a;
    while (y<b)
        count = count + 1;
        del = eval(H);
        y = y + del;
        yvec = [yvec y];
    end
    M = count;                          % has to be set here
    c = (b-a)/(y-a);                    % normalisaton factor
    x = zeros(1,M)+1;
    x(1) = a;
    for i=1:M
        y = yvec(i);                    % H may use y as a variable
        del = c*eval(H);
        x(i+1) = x(i) + del;
    end

    % visualisation of the grid
    clear y;
    y = linspace(a,b,100);
    plot(x,zeros(size(x)),'bx',y,eval(H),'go');
    xlim([a b]);
  otherwise
    error(['Error!\nGrid calculation only numerical (num) or ' ...
          'analytical (ana).\nYou entered %s'], grid_calculation);
end

switch method
  case 'Lag1'
    Lagrange1
  case 'Lag2'
    Lagrange2
  case 'spl1'
    Spline1
  case 'spl2'
    Spline2
  otherwise
    error(['Error!\nMethod %s not supported. Use Lag1, Lag2, spl1, ' ...
           'or spl2 instead'], method);
end
