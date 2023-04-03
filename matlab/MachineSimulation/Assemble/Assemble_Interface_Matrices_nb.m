function [S1,S2] = Assemble_Interface_Matrices_nb(p,eI,t,a,gamma)
% S1:\int gamma/h [u][v] ds; S2: \int {a dn u}[v] ds;  

np=size(p,2); neI=size(eI,2); 


if length(gamma)==1, gamma=gamma*ones(1,neI); end
if length(a)==1, a=a*ones(1,neI); end

S1=sparse(np,np); S2=sparse(np,np); 

for i = 1:size(eI,2)
    e1 = eI(1,i); e2 = eI(2,i);
    p1 = p(:,e1); p2 = p(:,e2); 
    
    t1 = eI(3,i); t2 = eI(4,i);
    tr1 = t(:,t1); tr2 = t(:,t2); % triangles
    stPTr1 = eI(5,i); stPTr2 = eI(6,i); % starting points of triangles (first
    % second or third node in the triangle list)
    
    % position of the triangle points in p
    tp1 = tr1(stPTr1); tp2 = tr1(mod(stPTr1,3)+1); tp3 = tr1(mod(stPTr1+1,3)+1);
    tp4 = tr2(stPTr2); tp5 = tr2(mod(stPTr2,3)+1); tp6 = tr2(mod(stPTr2+1,3)+1);
    
    pt1 = p(:,tp1); pt2 = p(:,tp2); pt3 = p(:,tp3);
    pt4 = p(:,tp4); pt5 = p(:,tp5); pt6 = p(:,tp6);

    lambda1 = eI(7,i); lambda2 = eI(8,i); % relative positions of the
    lambda3 = eI(9,i); lambda4 = eI(10,i);% points on the two edges

    B1 = [pt2-pt1,pt3-pt1]; iB1t = inv(B1)'; % transposed!!
    B2 = [pt5-pt4,pt6-pt4]; iB2t = inv(B2)'; 

    dphih1=[-1;-1]; dphih2=[1;0]; dphih3=[0;1];           % reference basis functions
    dphi1=iB1t*dphih1; dphi2=iB1t*dphih2; dphi3=iB1t*dphih3;% mapped basis functions
    dphi4=iB2t*dphih1; dphi5=iB2t*dphih2; dphi6=iB2t*dphih3;% mapped basis functions
%         dphi1=dphih1'*iB1t; dphi2=dphih2'*iB1t; dphi3=dphih3'*iB1t;% mapped basis functions
%     dphi4=dphih1'*iB2t; dphi5=dphih2'*iB2t; dphi6=dphih3'*iB2t;% mapped basis functions
    
    orient = det([pt1,pt2,pt3;1,1,1]); % orientation
    tv = pt1 - pt2; % tangential vector
    if orient > 0 % counterclockwise
        n = [-tv(2);tv(1)];
    else % clockwise
        n = [tv(2);-tv(1)];
    end
    n = n/norm(n);

    dphin1 = dot(dphi1,n); dphin2 = dot(dphi2,n); dphin3 = dot(dphi3,n);
    dphin4 = dot(dphi4,n); dphin5 = dot(dphi5,n); dphin6 = dot(dphi6,n);

    p2_1 = p2-p1; 
    len_eI = sqrt(p2_1(1,1)^2+p2_1(2,1)^2);
    
    pt2_1 = pt2 - pt1;
    len_et1 = sqrt(pt2_1(1,1)^2+pt2_1(2,1)^2);
    
    %------------------------------

    gammai = gamma(i);
    ai = a(i);

    elmat1 = zeros(4,4);
    elmat2 = zeros(6,6);
    irx=[0.0,0.5,1.0]; w=[1,4,1]/6;      % Simpson rule for 1d integral
    for l=1:length(w)
        xh1=irx(l)*(lambda2-lambda1)+lambda1;         % integration point
        xh2=irx(l)*(lambda4-lambda3)+lambda3;         % integration point
        yh = 0;

        phi1=1-xh1-yh; phi2=xh1; phi3=yh;  % value of basis functions
        phi4=1-xh2-yh; phi5=xh2; phi6=yh;  % value of basis functions
        
        % \int gamma/h [u][v] ds
        elmat3 = [phi1*phi1, phi1*phi2, -phi1*phi4, -phi1*phi5; ...
                  phi2*phi1, phi2*phi2, -phi2*phi4, -phi2*phi5; ...
                 -phi4*phi1, -phi4*phi2, phi4*phi4, phi4*phi5; ...
                 -phi5*phi1, -phi5*phi2, phi5*phi4, phi5*phi5]*len_eI*w(l);  
        % \int {a dn u}[v] ds
        elmat4 = [ phi1*dphin1,   phi1*dphin2,  phi1*dphin3,   phi1*dphin4,  phi1*dphin5,  phi1*dphin6; ...
                   phi2*dphin1,   phi2*dphin2,  phi2*dphin3,   phi2*dphin4,  phi2*dphin5,  phi2*dphin6; ...
                   phi3*dphin1,   phi3*dphin2,  phi3*dphin3,   phi3*dphin4,  phi3*dphin5,  phi3*dphin6; ...
                   -phi4*dphin1, -phi4*dphin2, -phi4*dphin3,  -phi4*dphin4, -phi4*dphin5, -phi4*dphin6; ...
                   -phi5*dphin1, -phi5*dphin2, -phi5*dphin3,  -phi5*dphin4, -phi5*dphin5, -phi5*dphin6; ...
                   -phi6*dphin1, -phi6*dphin2, -phi6*dphin3,  -phi6*dphin4, -phi6*dphin5, -phi6*dphin6] ...
                     *len_eI*w(l)/2*ai; % /2 is coming from the average

        elmat1 = elmat1 + elmat3.*gammai/len_et1;
        elmat2 = elmat2 + elmat4;

    end

    ii=[tp1;tp2;tp4;tp5];
    S1(ii,ii) = S1(ii,ii) + elmat1;

    jj = [tp1;tp2;tp3;tp4;tp5;tp6];
    S2(jj,jj) = S2(jj,jj) + elmat2;

end
end

