function A = Plotmesh(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4)
    X = [X1, X2, X3, X4, X1]';
    Y = [Y1, Y2, Y3, Y4, Y1]';
    Z = [Z1, Z2, Z3, Z4, Z1]';

% fid = fopen('input.txt','r');
% data = fscanf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g', [18, inf]);
% data = data';

    hold on;
    A = plot3(X, Y, Z, 'Color','g')
    alpha 0.5;
    hold off;
end
