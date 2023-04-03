function [alpha,iter,energy] = ArmijoLineSearch(phi,phi_der0,alpha_start)
    delta = 0.5;
    gamma = 1e-4;
    alpha = alpha_start;
    maxiter = 10^4;
    phi0 = phi(0);
    
    for iter = 1:maxiter
        energy = phi(alpha);
        if (energy <= phi0 + gamma*alpha*phi_der0)
            return;
        else
            alpha = alpha*delta;
        end    
    end
end

