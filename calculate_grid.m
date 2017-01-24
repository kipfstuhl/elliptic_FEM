
function [x, M] = calculate_grid (a, b, M, H, method)
% calculates a grid for one dimension
% a, b are the endpoints
% M is the desired number of elements
% H is the spacing function
% method is 'num' for numeric or 'ana' for analytical calculation


switch method
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
    x = zeros(1,M+1);
    x(1) = a;
    for i=1:M
        y = yvec(i);                    % H may use y as a variable
        del = c*eval(H);
        x(i+1) = x(i) + del;
    end
    clear y;
    clear yvec;
  otherwise
    error(['Error!\nGrid calculation only numerical (num) or ' ...
          'analytical (ana).\nYou entered %s'], grid_calculation);
end

