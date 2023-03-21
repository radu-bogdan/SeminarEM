function [hessian_g] = hessian_g(indices,grad_uh,material)
    rotationMatrix1 = [0,1;-1,0];
    rot_duh = rotationMatrix1*grad_uh;
    hessian_f_res = hessian_f(indices,rot_duh,material);
    hessian_g = zeros(2,2,length(indices));
    hessian_g(1,1,:) = hessian_f_res(2,2,:);
    hessian_g(2,1,:) = -hessian_f_res(1,2,:);
    hessian_g(1,2,:) = -hessian_f_res(2,1,:);
    hessian_g(2,2,:) = hessian_f_res(1,1,:);
end

