function [F1,F2] = Assemble_F_Vec(p,t,e,uh,B,m,material)
% [R,N] = assemr_fast(p,e,d,h)
%    assemble Robin matrix and Neumann vector
%    for details see assemr_tb
%
%    This is the fast implementation avoiding 
%    loop over elements

np=size(p,2); ne=size(e,2); nt=size(t,2); 

% mapping 
e1=e(1,:);  e2=e(2,:);

% Neumann vector
elvecs = zeros(2,ne);           % every column has one element vector
irx=[0.0,0.5,1.0]; w=[1,4,1]/6;
for l=1:length(w)
    xh=irx(l);            % integration point
    phi1=1-xh; phi2=xh;   % value of basis functions
    hlw = ( uh(e1)' * (1-irx(l)) + uh(e2)' *  irx(l)) .*B.len*w(l);
    elvecs = elvecs + [ phi1*hlw; 
                        phi2*hlw];
end
ii=[e1;e2]; % row indices
F1 = accumarray(ii(:),elvecs(:),[np,1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[duhx,duhy] = pdegrad(p,t,uh);
grad_uh = [duhx;duhy];
g_grad_val = grad_g(t(4,:),grad_uh,material);
g_grad_valx = g_grad_val(1,:);
g_grad_valy = g_grad_val(2,:);
elvecs = zeros(3,nt);                            % element vector
irx=[1/3]; iry=[1/3]; w=[1/2];
for l=1:length(w) % numerical integration
    fvw = B.vol*w(l);
    elvecs = elvecs + [...
                  (g_grad_valx.*B.dphi1x + g_grad_valy.*B.dphi1y).*fvw; 
                  (g_grad_valx.*B.dphi2x + g_grad_valy.*B.dphi2y).*fvw; 
                  (g_grad_valx.*B.dphi3x + g_grad_valy.*B.dphi3y).*fvw];
end
t1=t(1,:); t2=t(2,:); t3=t(3,:);
ii=[t1;t2;t3];                          % row indices
F2=accumarray(ii(:),elvecs(:),[np,1]);   % assemble global vector
end

