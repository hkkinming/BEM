function [X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4] = Grid(P1, P2, P3, P4, D)
    D = D + 1;
    XMax = max([P1(1), P2(1), P3(1), P2(1)]);
    XMin = min([P1(1), P2(1), P3(1), P2(1)]);
    YMax = max([P1(2), P2(2), P3(2), P2(2)]);
    YMin = min([P1(2), P2(2), P3(2), P2(2)]);
    ZMax = max([P1(3), P2(3), P3(3), P2(3)]);
    ZMin = min([P1(3), P2(3), P3(3), P2(3)]);
    
    X = linspace(XMin, XMax, D(1));
    Y = linspace(YMin, YMax, D(2));
    Z = linspace(ZMin, ZMax, D(3));
    
    if D(1) == 1
        [Z1, Y1] = meshgrid(Z(1 :end - 1), Y(1 :end - 1));
        [Z2, Y2] = meshgrid(Z(2 :end), Y(1 :end - 1));
        [Z3, Y3] = meshgrid(Z(2 :end), Y(2 :end));
        [Z4, Y4] = meshgrid(Z(1 :end - 1), Y(2 :end));
        Z1 = Z1(:);
        Y1 = Y1(:);
        X1 = repmat(X, size(Z1));
        Z2 = Z2(:);
        Y2 = Y2(:);
        X2 = repmat(X, size(Z2));
        Z3 = Z3(:);
        Y3 = Y3(:);
        X3 = repmat(X, size(Z3));
        Z4 = Z4(:);
        Y4 = Y4(:);
        X4 = repmat(X, size(Z4));
    elseif D(2) == 1
        [X1, Z1] = meshgrid(X(1 :end - 1), Z(1 :end - 1));
        [X2, Z2] = meshgrid(X(2 :end), Z(1 :end - 1));
        [X3, Z3] = meshgrid(X(2 :end), Z(2 :end));
        [X4, Z4] = meshgrid(X(1 :end - 1), Z(2 :end));
        X1 = X1(:);
        Z1 = Z1(:);
        Y1 = repmat(Y, size(X1));
        X2 = X2(:);
        Z2 = Z2(:);
        Y2 = repmat(Y, size(X2));
        X3 = X3(:);
        Z3 = Z3(:);
        Y3 = repmat(Y, size(X3));
        X4 = X4(:);
        Z4 = Z4(:);
        Y4 = repmat(Y, size(X4));
    elseif D(3) == 1
        [Y1, X1] = meshgrid(Y(1 :end - 1), X(1 :end - 1));
        [Y2, X2] = meshgrid(Y(2 :end), X(1 :end - 1));
        [Y3, X3] = meshgrid(Y(2 :end), X(2 :end));
        [Y4, X4] = meshgrid(Y(1 :end - 1), X(2 :end));
        Y1 = Y1(:);
        X1 = X1(:);
        Z1 = repmat(Z, size(Y1));
        Y2 = Y2(:);
        X2 = X2(:);
        Z2 = repmat(Z, size(Y2));
        Y3 = Y3(:);
        X3 = X3(:);
        Z3 = repmat(Z, size(Y3));
        Y4 = Y4(:);
        X4 = X4(:);
        Z4 = repmat(Z, size(Y4));
end
