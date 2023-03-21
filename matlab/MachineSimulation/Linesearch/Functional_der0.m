function [integral] = Functional_der0(dir,gamma,J_j,J_m,F1_vec,F2_vec,material)
    nu0 = material.nu0;
    integral = 0;
    integral = integral - J_j'*dir + nu0*J_m'*dir;
    integral = integral + gamma*F1_vec'*dir + F2_vec'*dir;   
end

