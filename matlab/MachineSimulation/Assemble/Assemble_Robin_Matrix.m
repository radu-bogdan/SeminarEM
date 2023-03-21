function [R]= Assemble_Robin_Matrix(p,e,B)
% [R]= Assemble_EdgeElement_Matrix(p,e)
%    assemble Robin matrix
np=size(p,2); ne=size(e,2); 

% mapping 
e1=e(1,:);  e2=e(2,:);
len = B.len;

% Robin matrix
elmats = zeros(4,ne);               % every column stores one element matrix
irx=[0.0,0.5,1.0]; w=[1,4,1]/6;     % Simpson rule for 1d integral
for l=1:length(w)
    xh=irx(l);            % integration point
    phi1=1-xh; phi2=xh;   % value of basis functions
    dlw = len*w(l);
    elmats = elmats + [ phi1*phi1*dlw;
                        phi1*phi2*dlw; 
                        phi2*phi1*dlw;
                        phi2*phi2*dlw];
end
ii=[e1;e1;e2;e2]; % row indices
jj=[e1;e2;e1;e2]; % col indices
R = sparse(ii(:),jj(:),elmats(:),np,np);

