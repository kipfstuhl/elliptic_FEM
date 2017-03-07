%% Plot basis functions for product spaces
%  polynomial of deg N+1 in x firection and N in y direction

if ~exist('N')
    N = 4;
end

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

% preparation for plotting the interpolation points
xvals = repmat(xnodes, size(ynodes, 1), 1);
yvals = repelem(ynodes, size(xnodes, 1));
zvals = zeros(size(xvals));

% speed up the plotting by updating just the necessary data
plot3(xvals, yvals, zvals, 'r.', 'MarkerSize', 15);
hold on
h = surf(xx, yy, zz{1,1});
xlabel('x');
ylabel('y');
xlim([-1 1]);
ylim([-1 1]);
hold off

for i=1:N+1
    for j=1:N
        set(h, 'ZData', zz{i, j});
        waitforbuttonpress;
    end
end