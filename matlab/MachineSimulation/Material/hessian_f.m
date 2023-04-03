function [hessian_f] = hessian_f(indices,rot_duh,material)

    %Nonconductive hessian matrix
    nu0 = (10^7)/(4*pi); %reluctivity for air and copper
    non_cond_pos = ismember(indices,material.all_non_conductive);
    nt = length(non_cond_pos);
    non_cond_mat = zeros(2,2,nt);
    non_cond_mat(1,1,:) = nu0*non_cond_pos;
    non_cond_mat(2,2,:) = nu0*non_cond_pos;
    
    %Iron hessian matrix
    k1 = 49.4;
    k2 = 1.46;
    k3 = 520.6;
    norm = rot_duh(1,:).*rot_duh(1,:) + rot_duh(2,:).*rot_duh(2,:);
    iron_pos = ismember(indices,material.iron);
    factor = (k1*exp(k2*norm)).*iron_pos;
    mat1 = zeros(2,2,nt);
    mat1(1,1,:) = factor+ iron_pos*k3; mat1(2,2,:) = factor+iron_pos*k3;
    mat2 = zeros(2,2,nt);
    mat2(1,1,:) = factor.*rot_duh(1,:).*rot_duh(1,:);
    mat2(2,1,:) = factor.*rot_duh(2,:).*rot_duh(1,:);
    mat2(1,2,:) = factor.*rot_duh(1,:).*rot_duh(2,:);
    mat2(2,2,:) = factor.*rot_duh(2,:).*rot_duh(2,:);
    iron_mat = mat1 + k2*2*mat2;
    hessian_f = non_cond_mat + iron_mat;
end

