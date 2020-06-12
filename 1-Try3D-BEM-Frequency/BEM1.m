clc;
clear;
close all;
addpath('C:\Users\Chan Kin Ming\Desktop\Current Working\Determination of Modulus of Pavement\BEM Matlab\1-Try3D-BEM-Frequency');



% Program to implement the Boundary Element Method for the Elastostatics Problem with given boundary conditions
% ======================================================


% ----- Pre-processing -----
% Input Data
fid = fopen('input.txt','r');
data = fscanf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g', [18, inf]);
data = data';
m = 6;
n = 21;
mn3 = m .* n .* 3;
m3 = m .* 3;

for r = 0: (mn3 - 1)
    a = mod(r, m);
    b = floor(mod(r, m .* n) ./ m);
    i = floor(r ./ (m .* n));
    index = mod(r, m) + floor(mod(r, m .* n) ./ m) .* m + 1;
    E(a + 1, b + 1, i + 1, :) = [data(index, 1: 12), data(index, 13: 14), data(index, 15: 15), data(index, 16 + i: 16 + i)];
    
end

Ef(:, :, :, 1: 15) = E(:, :, :, 1: 15);

for a = 0: mod(r, m)
    for i = 0: 2
        Ef(a + 1, :, i + 1, 16) = fft(E(a + 1, :, i + 1, 16));
    end
end


% ----- Processing -----
% b.*j = m.*3, u(x1,k,1)u(x2,k,1)...u(xm,k,1)u(x1,k,2)...u(xm-1,k,3)u(xm,k,3)

parfor k = 1: n
    A = zeros(m3, m3);
    B = zeros(m3);
    datatype = zeros(1, m3);
    datavalue = zeros(1, m3);
    UM = zeros(m3, m3);
    PM = zeros(m3, m3);
    for r = 0: (m3 - 1)% (a, i) is the origin
        a = mod(r, m) + 1;
        i = floor(r ./ m) + 1;
        datatype(r + 1) = Ef(a, k, i, 15);
        datavalue(r + 1) = Ef(a, k, i, 16);
        Ur = zeros(1, m3);
        Pr = zeros(1, m3);
        for c = 0: (m3 - 1)
            b = mod(c, m) + 1;
            j = floor(c ./ m) + 1;
            if((a ~= b) || (i == j))
                [temp1, temp2] = UP(E, b, j, reshape((Ef(a, k, i, 1: 3) + E(a, k, i, 4: 6) + E(a, k, i, 7: 9) + E(a, k, i, 10: 12)) ./ 4.0, [1,3]), i , k);
                Ur(c + 1) = temp1;
                Pr(c + 1) = temp2;
            end
        end
        UM(r + 1, :) = Ur;
        PM(r + 1, :) = Pr;
    end
    PM = PM - 0.5 .* eye(m3);
    A = PM .* datatype + UM .* (1.0 - datatype);
    B = - (PM .* (1.0 - datatype) + UM .* datatype) * datavalue';
    Z = A \ B;
    UMk(k,:,:) = UM;
    PMk(k,:,:) = PM;
    Ak(k,:,:) = A;
    Bk(k,:) = B;
    Zk(k,:) = Z;
    U(k,:) = (1 - datatype') .* datavalue' + datatype' .* Z;
    P(k,:) = (1 - datatype') .* Z + datatype' .* datavalue';
end
plot(real(ifft(U)))
plot(real(ifft(P)))
% % % ----- Post-Processing -----
% % def output(Data, m, n):
% %     data = []
% %     for a in range(m):
% %         datat = []
% %         for b in range(n):
% %             datati = []
% %             for i in range(3):
% %                 datati.append(Data[a + b .* m + i .* n .* m])
% %             datat.append(datati)
% %         data.append(datat)
% %     return data
% % u , p = output(U, m, n), output(P, m, n)
% %
% % np.savetxt('outputA.txt', A)
% % np.savetxt('outputB.txt', B)
% % np.savetxt('outputU.txt', U)
% % np.savetxt('outputP.txt', P)
% % np.savetxt('outputZ.txt', Z)
% % np.savetxt('outputUM.txt', UM)
% % np.savetxt('outputPM.txt', PM)
% %
% % print("done")
% %
% % % z = np.linalg.solve(A,B)
% % % t, v = np.array([e.t for e in E]), np.array([e.v for e in E])
% % % u, q = np.empty(n), np.empty(n)
% % % u = (1-t).*v+t.*z
% % % q = (1-t).*z+t.*v
% 
% filename = 'Result.xlsx';
% writematrix(U, filename, 'Sheet', 1, 'Range', 'A:A');
% writematrix(P, filename, 'Sheet', 2, 'Range', 'A:A');
% writematrix(UM, filename, 'Sheet', 3);
% writematrix(PM, filename, 'Sheet', 4);
