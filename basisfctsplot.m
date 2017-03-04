%% Plot basis functions for product spaces
%  polynomial of deg N+1 in x firection and N in y direction

N = 4;

[xnodes,~,~] = lglnodes(N);
[ynodes,~] = lgwt(N, -1, 1);           % returns also the
                                       % Vandermonde matrix

f = zeros(N,N+1);
for i=1:N+1
    fvals = zeros(N+1, 1);
    fvals(i) = 1.0;
    f(i,:) = polyfit(xnodes, fvals, N);
end

g = zeros(N-1, N);
for i=1:N
    gvals = zeros(N, 1);
    gvals(i) = 1.0;
    g(i,:) = polyfit(ynodes, gvals, N-1);
end

% make product functions f(i)*g(i)
%h = cell(N+1, N);
[xx, yy] = meshgrid(linspace(-1, 1, 100), linspace(-1, 1, 100));
%zz = zeros(N, N-1, 100, 100);
zz = cell(N+1, N);
for i=1:N+1
    for j=1:N
        %h_temp = @(x, y) polyval(f(i,:), x) * polyval(g(j,:), y);
        %h{i, j} = h_temp;
        zz{i, j} = polyval(f(i,:), xx) .* polyval(g(j,:), yy);
    end
end


for i=1:N+1
    for j=1:N
        surf(xx, yy, zz{i,j});
        xlabel('x');
        ylabel('y');
        waitforbuttonpress;
    end
end