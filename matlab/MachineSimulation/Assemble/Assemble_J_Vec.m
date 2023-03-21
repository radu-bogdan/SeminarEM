function [J_j,J_m] = Assemble_J_Vec(p,t,J_vec,M_vec,B)
nt=size(t,2); np=size(p,2);
t1=t(1,:); t2=t(2,:); t3=t(3,:);
vol = B.vol;

% load vector
elvecs = zeros(3,nt);                            % element vector
irx=[0.5,0.5,0]; iry=[0,0.5,0.5]; w=[1,1,1]/6;   % edge midpoint rule
for l=1:length(w)                                % numerical integration
    xh=irx(l); yh=iry(l);                        % integration point
    phi1=1-xh-yh; phi2=xh; phi3=yh;              % value of basis functions
    fvw = J_vec.*vol*w(l);
    elvecs = elvecs + [...
                  phi1*fvw; 
                  phi2*fvw;
                  phi3*fvw];
end
ii=[t1;t2;t3];                          % row indices
J_j=accumarray(ii(:),elvecs(:),[np,1]);   % assemble global vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elvecs = zeros(3,nt);                            % element vector
irx=[1/3]; iry=[1/3]; w=[1/2];
mx = M_vec(1,:);
my = M_vec(2,:);
for l=1:length(w) % numerical integration
    fvw = B.vol*w(l);
    elvecs = elvecs + [...
                  (mx.*B.dphi1y - my.*B.dphi1x).*fvw; 
                  (mx.*B.dphi2y - my.*B.dphi2x).*fvw; 
                  (mx.*B.dphi3y - my.*B.dphi3x).*fvw];
end
t1=t(1,:); t2=t(2,:); t3=t(3,:);
ii=[t1;t2;t3];                          % row indices
J_m=accumarray(ii(:),elvecs(:),[np,1]);   % assemble global vector
end

