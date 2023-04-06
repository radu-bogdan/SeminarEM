function [integral] = Functional_Nitsche(t,uh,gamma,J_j,J_m,F1_vec,S,material,B)

integral = Functional(t,uh,gamma,J_j,J_m,F1_vec,material,B) + uh'*S*uh/2;

end
