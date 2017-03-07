
%myarray = [0.5, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01];
myarray = [0.2, 0.1, 0.075, 0.05, 0.025, 0.01];
H1Arr = zeros(size(myarray));
L2Arr = zeros(size(myarray));

set(0,'DefaultFigureVisible','off');    % disable graphics
for k=1:size(myarray, 2)
    gridSpace = num2str(myarray(k));

    elliptic_FEM

    H1Arr(k) = H1;
    L2Arr(k) = L2;
    
end
set(0, 'DefaultFigureVisible', 'on');   % enable graphics again

figure
loglog(myarray, L2Arr, 'DisplayName', 'L2 Error');
hold on
loglog(myarray, H1Arr, 'DisplayName', 'H1 Error');
legend('show');
hold off

