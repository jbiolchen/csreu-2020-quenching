function [u0, flag, iter] = jac_fsolve(f, Jf, u0, err, maxiter)
% Newton's Method with Jacobian; equivalent usage as fsolve with Jacobian
% gmres defaults to 10 iterations max; currently can't specify max Newton iterations

res = f(u0);
flag = 0; % initialize as "Newton's method didn't converge"
for iter = 0:maxiter
    if norm(res) < err
        flag = 1; % Newton's method converged
        break
    end
    jac = @(du0) Jf(du0, u0);
    step = gmres(jac, res);
    u0 = u0 - step;
    res = f(u0);
end
iter
norm(res)

end