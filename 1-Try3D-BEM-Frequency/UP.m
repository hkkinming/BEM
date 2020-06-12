function [U, P] = UP(E, b, j, x, i ,k1)% (E->kl, x->alpha, t->beta, i, j)
    d0 = 1;
    mu = 1;
    upsilon = 1;
    rho = 1;
    c1 = 3000;
    c2 = c1 * ((1 - 2 * 0.35) / (2 - 2 *0.35)) ^ 0.5;
    e = reshape(E(b, k1, j, :), [16, 1])
    %Q = [[zeros(3), zeros(3), eye(3)]; [eye(3), zeros(3), eye(3)]; [eye(3), eye(3), eye(3)]; [zeros(3), eye(3), eye(3)]];
    %tranformCoeff = (Q' * Q) ^ (-1) * Q';
    %y1 = @(nu, xi) tranformCoeff(1: 1, :) * e(1: 12) .* nu + tranformCoeff(4: 4, :) * e(1: 12) .* xi + tranformCoeff(7: 7, :) * e(1: 12);
    %y2 = @(nu, xi) tranformCoeff(2: 2, :) * e(1: 12) .* nu + tranformCoeff(5: 5, :) * e(1: 12) .* xi + tranformCoeff(8: 8, :) * e(1: 12);
    %y3 = @(nu, xi) tranformCoeff(3: 3, :) * e(1: 12) .* nu + tranformCoeff(6: 6, :) * e(1: 12) .* xi + tranformCoeff(9: 9, :) * e(1: 12);
    y1 = @(nu, xi) (e(4) - e(1)) .* nu + (e(10) - e(1)) .* xi + e(1);
    y2 = @(nu, xi) (e(5) - e(2)) .* nu + (e(11) - e(2)) .* xi + e(2);
    y3 = @(nu, xi) (e(6) - e(3)) .* nu + (e(12) - e(3)) .* xi + e(3);
    r1 = @(nu, xi) y1(nu, xi) - x(1);
    r2 = @(nu, xi) y2(nu, xi) - x(2);
    r3 = @(nu, xi) y3(nu, xi) - x(3);
    r = @(nu, xi) (r1(nu, xi) .^ 2.0 + r2(nu, xi) .^ 2.0 + r3(nu, xi) .^ 2.0) .^0.5;
    dri = @(nu, xi) ((i == 1) .* r1(nu, xi) + (i == 2) .* r2(nu, xi) + (i == 3) .* r3(nu, xi)) ./ r(nu , xi);
    drj = @(nu, xi) ((j == 1) .* r1(nu, xi) + (j == 2) .* r2(nu, xi) + (j == 3) .* r3(nu, xi)) ./ r(nu , xi);
    yn = cross(e(7: 9) - e(4: 6), e(1: 3) - e(4: 6));
    dryn = @(nu, xi) (yn(1) .* r1(nu, xi) + yn(2) .* r2(nu, xi) +yn(3) .* r3(nu, xi)) ./ r(nu , xi);
    delta = @(nu, xi) nu == xi;
    deltaij = delta(i, j);
    
    k = 2.0 .* pi .* (k1 - 1);
    T = 2.0;
    
    if k == 0.0
        HA = @(nu, xi) 1.0 ./ c1 .^ 2.0 .* (T .^ 2.0 - (r(nu, xi) ./ c1) .^2.0) ./ 2.0 - 1.0 ./ c2 .^ 2.0 .* (T .^ 2.0 - (r(nu, xi) ./ c2) .^2.0) ./ 2.0;
    else
        HA = @(nu, xi) 1.0 ./ c1 .^ 2.0 .* (c1 ./ (1.0i .* r(nu , xi) .* k) + (c1 ./ (1.0i .* r(nu , xi) .* k)) .^ 2.0) .* exp (-1.0i .* r(nu , xi) .* k ./ c1) - 1.0 ./ c2 .^ 2.0 .* (c2 ./ (1.0i .* r(nu , xi) .* k) + (c2 ./ (1.0i .* r(nu , xi) .* k)) .^ 2.0) .* exp (-1.0i .* r(nu , xi) .* k ./ c2);
    end
    HB = @(nu, xi) 1.0 ./ c1 .^ 2.0 .* exp (-1.0i .* r(nu , xi) .* k ./ c1) - 1.0 ./ c2 .^ 2.0 .* exp (-1.0i .* r(nu , xi) .* k ./ c2);
    HC = @(nu, xi) 1.0 ./ c2 .^ 2.0 .* exp (-1.0i .* r(nu , xi) .* k ./ c2);
    HCp = @(nu, xi) (1.0i .* r(nu , xi) .* k ./ c2 + 1.0)  .* exp (-1.0i .* r(nu , xi) .* k ./ c2);
    HEp = @(nu, xi) (1.0i .* r(nu , xi) .* k ./ c1 + 1.0)  .* exp (-1.0i .* r(nu , xi) .* k ./ c1);   
    HDp = @(nu, xi) (1.0i .* r(nu , xi) .* k ./ c1 .^ 3.0)  .* exp (-1.0i .* r(nu , xi) .* k ./ c1) - (1.0i .* r(nu , xi) .* k ./ c2 .^ 3.0)  .* exp (-1.0i .* r(nu , xi) .* k ./ c2);
    
    Au = @(nu, xi) 1.0 ./ r(nu, xi) .* (3.0 .* dri(nu, xi) .* drj(nu, xi) - deltaij) .* HA(nu, xi);
    Bu = @(nu, xi) 1.0 ./ r(nu, xi) .* dri(nu, xi) .* drj(nu, xi) .* HB(nu, xi);
    Cu = @(nu, xi) 1.0 ./ r(nu, xi) .* deltaij .* HC(nu, xi);
    Ap = @(nu, xi) 1.0 ./ r(nu, xi) .^ 2.0 .* 6.0 .* c2 .^ 2.0 .* (yn(j) .* dri(nu, xi) + yn(i) .* drj(nu, xi) + (deltaij - 5.0 .* dri(nu, xi) .* drj(nu, xi)) .* dryn(nu, xi)) .* HA(nu, xi);
	Bp = @(nu, xi) 1.0 ./ r(nu, xi) .^ 2.0 .* 2.0 .* c2 .^ 2.0 .* (yn(j) .* dri(nu, xi) + yn(i) .* drj(nu, xi) + (deltaij - 6.0 .* dri(nu, xi) .* drj(nu, xi)) .* dryn(nu, xi)) .* HB(nu, xi);
    Cp = @(nu, xi) 1.0 ./ r(nu, xi) .^ 2.0 .* -1.0 .* (yn(i) .* drj(nu, xi) + deltaij .* dryn(nu, xi)) .* HCp(nu, xi);
    Dp = @(nu, xi) 1.0 ./ r(nu, xi) .^ 2.0 .* -2.0 .* c2 .^ 2.0 .* dri(nu, xi) .* drj(nu, xi) .* dryn(nu, xi) .* HDp(nu, xi);
    Ep = @(nu, xi) 1.0 ./ r(nu, xi) .^ 2.0 .* (2.0 .* c2 .^ 2.0 - c1 .^ 2.0) .* yn(j) .* dri(nu, xi) .* HEp(nu, xi);
    
    Uf = @(nu, xi) (Au(nu, xi) + Bu(nu, xi) + Cu(nu, xi)) .* 0.25 ./ pi ./ rho;
    Pf = @(nu, xi) (Ap(nu, xi) + Bp(nu, xi) + Cp(nu, xi) + Dp(nu, xi) + Ep(nu, xi)) .* 0.25 ./ pi;
    
    U = integral2(@(nu, xi) Uf(nu, xi), 0.0, 0.5, 0.0, 1.0, "AbsTol", 1e-3);
    P = integral2(@(nu, xi) Pf(nu, xi), 0.0, 1.0, 0.0, 1.0, "AbsTol", 1e-3);

end
