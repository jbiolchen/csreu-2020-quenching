function [u0, flag, iter] = jac_fsolve(f, Jf, u0, err)
% Newton's Method with Jacobian; equivalent usage as fsolve with Jacobian

res = f(u0)
while norm(res) >= err
    jac = @(du0) Jf(du0, u0);
    [step, flag, ~, iter] = gmres(jac, res);
    
    pause
    step
    u0 = u0 - step
    res = f(u0)
end
flag = 1 - flag; % so that 1 = converged, 0 = didn't converge
iter = iter(1); % outer iteration
end