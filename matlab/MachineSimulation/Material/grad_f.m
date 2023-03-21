function [grad_f] = grad_f(indices,rot_duh,material)
nu0 = material.nu0;
k1 = material.k1;
k2 = material.k2;
k3 = material.k3;

norm = rot_duh(1,:).*rot_duh(1,:) + rot_duh(2,:).*rot_duh(2,:);
grad_f_mult = nu0*ismember(indices,material.all_non_conductive) + (k1*exp(k2*norm)+k3).*ismember(indices,material.iron);
grad_f = [grad_f_mult;grad_f_mult].*rot_duh;
end

