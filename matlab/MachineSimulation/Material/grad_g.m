function [grad_g_rot] = grad_g(indices,grad_uh,material)
rotationMatrix1 = [0,1;-1,0];
rot_grad_uh = rotationMatrix1*grad_uh; %caculate roation of gradient
grad_f_vec = grad_f(indices,rot_grad_uh,material); %calculate grad_f of rotated gradient
rotationMatrix2 = [0,-1;1,0];
grad_g_rot = rotationMatrix2*grad_f_vec;   %rotates the results of grad_f
end

