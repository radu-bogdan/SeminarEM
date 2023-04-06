function [integral] = Functional_der0_Nitsche(dir,gamma,J_j,J_m,F1_vec,F2_vec,Suh,material)

    integral = Functional_der0(dir,gamma,J_j,J_m,F1_vec,F2_vec,material) + Suh'*dir;   
    
end

