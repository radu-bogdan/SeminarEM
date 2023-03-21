function [N] = Assemble_Newton_Matrix(p,t,uh,B,material)
nt=size(t,2); np=size(p,2);
t1=t(1,:); t2=t(2,:); t3=t(3,:);
grad_uh_ref =  [uh(t2) - uh(t1), uh(t3) - uh(t1)]';
grad_uh = [B.iB11.*grad_uh_ref(1,:) + B.iB21.*grad_uh_ref(2,:); B.iB12.*grad_uh_ref(1,:) + B.iB22.*grad_uh_ref(2,:)];
       
hessian = hessian_g(t(4,:),grad_uh,material);
    
elmats = zeros(9,nt);          % every colum contains matrix for one element
w=[1/2];
a = ones(1,nt);
for l=1:length(w) % numerical integration
        dphi1x = B.dphi1x;   dphi1y = B.dphi1y;
        dphi2x = B.dphi2x;   dphi2y = B.dphi2y;
        dphi3x = B.dphi3x;   dphi3y = B.dphi3y;
        hessian11 = squeeze(hessian(1,1,:))';
        hessian12 = squeeze(hessian(1,2,:))'; 
        hessian21 = squeeze(hessian(2,1,:))'; 
        hessian22 = squeeze(hessian(2,2,:))'; 
        
        
        avw = a.*B.vol*w(l); 
        elmats = elmats + [ ...
                (dphi1x.*hessian11.*dphi1x + dphi1y.*hessian22.*dphi1y + dphi1x.*hessian12.*dphi1y + dphi1y.*hessian21.*dphi1x).*avw;     % K11   
                (dphi1x.*hessian11.*dphi2x + dphi1y.*hessian22.*dphi2y + dphi1x.*hessian12.*dphi2y + dphi1y.*hessian21.*dphi2x).*avw;     % K12                
                (dphi1x.*hessian11.*dphi3x + dphi1y.*hessian22.*dphi3y + dphi1x.*hessian12.*dphi3y + dphi1y.*hessian21.*dphi3x).*avw;     % K13                 
                (dphi2x.*hessian11.*dphi1x + dphi2y.*hessian22.*dphi1y + dphi2x.*hessian12.*dphi1y + dphi2y.*hessian21.*dphi1x).*avw;     % K21                
                (dphi2x.*hessian11.*dphi2x + dphi2y.*hessian22.*dphi2y + dphi2x.*hessian12.*dphi2y + dphi2y.*hessian21.*dphi2x).*avw;     % K22              
                (dphi2x.*hessian11.*dphi3x + dphi2y.*hessian22.*dphi3y + dphi2x.*hessian12.*dphi3y + dphi2y.*hessian21.*dphi3x).*avw;     % K23   
                (dphi3x.*hessian11.*dphi1x + dphi3y.*hessian22.*dphi1y + dphi3x.*hessian12.*dphi1y + dphi3y.*hessian21.*dphi1x).*avw;     % K31   
                (dphi3x.*hessian11.*dphi2x + dphi3y.*hessian22.*dphi2y + dphi3x.*hessian12.*dphi2y + dphi3y.*hessian21.*dphi2x).*avw;     % K32               
                (dphi3x.*hessian11.*dphi3x + dphi3y.*hessian22.*dphi3y + dphi3x.*hessian12.*dphi3y + dphi3y.*hessian21.*dphi3x).*avw];    % K33 
end
ii=[t1;t1;t1; t2;t2;t2; t3;t3;t3];            % row indices
jj=[t1;t2;t3; t1;t2;t3; t1;t2;t3];            % col indices
N=sparse(ii(:),jj(:),elmats(:),np,np);        % assemble 
end

