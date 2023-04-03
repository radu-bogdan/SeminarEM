function [integral] = Functional(t,uh,gamma,J_j,J_m,F1_vec,material,B)
nu0 = material.nu0;
indices = t(4,:);
t1=t(1,:); t2=t(2,:); t3=t(3,:);
grad_uh_ref =  [uh(t2) - uh(t1), uh(t3) - uh(t1)]';
grad_uh = [B.iB11.*grad_uh_ref(1,:) + B.iB21.*grad_uh_ref(2,:); B.iB12.*grad_uh_ref(1,:) + B.iB22.*grad_uh_ref(2,:)];
g_val = g_func(indices,grad_uh,material);
integral = 0;    
w= [1/2];
vol = B.vol;
    for l=1:length(w) % numerical integration
        integral = integral + sum(g_val.*vol*w(1));
    end
integral = integral - J_j'*uh + nu0*J_m'*uh + (gamma/2)*F1_vec'*uh;
end
