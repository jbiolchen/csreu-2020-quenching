function [u0, flag, iter] = jac_fsolve(f, Jf, u0, err, maxiter)
% Newton's Method with Jacobian; equivalent usage as fsolve with Jacobian

res = f(u0);
flag = 0; % initialize as "Newton's method didn't converge"
for iter = 0:maxiter
    if norm(res) < err
        flag = 1; % Newton's method converged
        break
    end
    jac = @(du0) Jf(du0, u0);
    step = gmres(jac, res, min(length(u0)/2, 1e4), 1e-8);
    u0 = u0 - step;
    res = f(u0);
end
end