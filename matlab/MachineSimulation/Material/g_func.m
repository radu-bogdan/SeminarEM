function [ res ] = g_func(indices,grad_uh,material)
nu0 = material.nu0;
k1 = material.k1;
k2 = material.k2;
k3 = material.k3;

norm = grad_uh(1,:).*grad_uh(1,:) + grad_uh(2,:).*grad_uh(2,:);
res = (1/2)*nu0*ismember(indices,material.all_non_conductive).*norm + ( ((1/2)*(k1/k2)*exp(k2*norm)) + (1/2)*k3*norm).*ismember(indices,material.iron);
end

