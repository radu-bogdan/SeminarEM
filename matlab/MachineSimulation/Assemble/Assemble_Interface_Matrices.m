function [S1,S2] = Assemble_Interface_Matrices(p,eI,t,a,gamma)
% S1:\int gamma/h [u][v] ds; S2: \int {a dn u}[v] ds;  

np=size(p,2); neI=size(eI,2); 


if length(gamma)==1, gamma=gamma*ones(1,neI); end
if length(a)==1, a=a*ones(1,neI); end

e1 = eI(1,:); e2 = eI(2,:);
p1 = p(:,e1); p2 = p(:,e2);

t1 = eI(3,:); t2 = eI(4,:);
tr1 = t(:,t1); tr2 = t(:,t2); % triangles
P1t1 = eI(5,:); P1t2 = eI(6,:); % starting points of triangles (first,
% second or third node in the triangle list)
P2t1 = mod(P1t1,3)+1; P3t1 = mod(P2t1,3)+1;
P2t2 = mod(P1t2,3)+1; P3t2 = mod(P2t2,3)+1;
P1t1(P1t1 == 2) = neI+1;
P1t1(P1t1 == 3) = 2*neI+1;
P1t2(P1t2 == 2) = neI+1;
P1t2(P1t2 == 3) = 2*neI+1;

P2t1(P2t1 == 2) = neI+1;
P2t1(P2t1 == 3) = 2*neI+1;
P2t2(P2t2 == 2) = neI+1;
P2t2(P2t2 == 3) = 2*neI+1;

P3t1(P3t1 == 2) = neI+1;
P3t1(P3t1 == 3) = 2*neI+1;
P3t2(P3t2 == 2) = neI+1;
P3t2(P3t2 == 3) = 2*neI+1;

P1t1 = P1t1 + (0:(neI-1));
P2t1 = P2t1 + (0:(neI-1));
P3t1 = P3t1 + (0:(neI-1));
P1t2 = P1t2 + (0:(neI-1));
P2t2 = P2t2 + (0:(neI-1));
P3t2 = P3t2 + (0:(neI-1));
    
% position of the triangle points in p
tr1 = tr1';
tr2 = tr2';
tp1 = tr1(P1t1); tp2 = tr1(P2t1); tp3 = tr1(P3t1);
tp4 = tr2(P1t2); tp5 = tr2(P2t2); tp6 = tr2(P3t2);

pt1 = p(:,tp1); pt2 = p(:,tp2); pt3 = p(:,tp3);
pt4 = p(:,tp4); pt5 = p(:,tp5); pt6 = p(:,tp6);

lambda1 = eI(7,:); lambda2 = eI(8,:); % relative positions of the
lambda3 = eI(9,:); lambda4 = eI(10,:);% points on the two edges


dphih1=[-1;-1]; dphih2=[1;0]; dphih3=[0;1];

B11=p(1,tp2)-p(1,tp1);
B12=p(1,tp3)-p(1,tp1);
B21=p(2,tp2)-p(2,tp1);
B22=p(2,tp3)-p(2,tp1);
detB = B11.*B22-B12.*B21;
iB11= B22./detB;
iB12=-B12./detB; 
iB21=-B21./detB;
iB22= B11./detB;
dphi1x= iB11*dphih1(1) + iB21*dphih1(2);   dphi1y= iB12*dphih1(1)+ iB22*dphih1(2);
dphi2x= iB11*dphih2(1) + iB21*dphih2(2);   dphi2y= iB12*dphih2(1)+ iB22*dphih2(2);
dphi3x= iB11*dphih3(1) + iB21*dphih3(2);   dphi3y= iB12*dphih3(1)+ iB22*dphih3(2);

B11=p(1,tp5)-p(1,tp4);
B12=p(1,tp6)-p(1,tp4);
B21=p(2,tp5)-p(2,tp4);
B22=p(2,tp6)-p(2,tp4);
detB = B11.*B22-B12.*B21;
iB11= B22./detB;
iB12=-B12./detB; 
iB21=-B21./detB;
iB22= B11./detB;
dphi4x= iB11*dphih1(1) + iB21*dphih1(2);   dphi4y= iB12*dphih1(1)+ iB22*dphih1(2);
dphi5x= iB11*dphih2(1) + iB21*dphih2(2);   dphi5y= iB12*dphih2(1)+ iB22*dphih2(2);
dphi6x= iB11*dphih3(1) + iB21*dphih3(2);   dphi6y= iB12*dphih3(1)+ iB22*dphih3(2);


orient = pt1(1,:).*pt2(2,:) + pt2(1,:).*pt3(2,:) + pt3(1,:).*pt1(2,:) ...
       - pt1(1,:).*pt3(2,:) - pt2(1,:).*pt1(2,:) - pt3(1,:).*pt2(2,:);
tv = pt1 - pt2; % tangential vector
orient = orient > 0;
n = zeros(2,neI);
n(:,orient) = [-tv(2,orient);tv(1,orient)];
orient = ~orient;
n(:,orient) = [tv(2,orient);-tv(1,orient)];

n = n./vecnorm(n);

dphin1 = dot([dphi1x;dphi1y],n); dphin2 = dot([dphi2x;dphi2y],n); dphin3 = dot([dphi3x;dphi3y],n);
dphin4 = dot([dphi4x;dphi4y],n); dphin5 = dot([dphi5x;dphi5y],n); dphin6 = dot([dphi6x;dphi6y],n);

p2_1 = p2-p1; 
len_eI = sqrt(p2_1(1,:).^2+p2_1(2,:).^2);

pt2_1 = pt2 - pt1;
len_et1 = sqrt(pt2_1(1,:).^2+pt2_1(2,:).^2);

%------------------------------

elmats1 = zeros(16,neI);
elmats2 = zeros(36,neI);
irx=[0.0,0.5,1.0]; w=[1,4,1]/6;      % Simpson rule for 1d integral
for l=1:length(w)
    xh1=irx(l)*(lambda2-lambda1)+lambda1;         % integration point
    xh2=irx(l)*(lambda4-lambda3)+lambda3;         % integration point
    yh = 0;

    phi1=1-xh1-yh; phi2=xh1; phi3=yh;  % value of basis functions
    phi4=1-xh2-yh; phi5=xh2; phi6=yh;  % value of basis functions
    
    % \int gamma/h [u][v] ds
    elmats1 = elmats1 + ...
             [phi1.*phi1; phi1.*phi2; -phi1.*phi4; -phi1.*phi5; ...
              phi2.*phi1; phi2.*phi2; -phi2.*phi4; -phi2.*phi5; ...
             -phi4.*phi1; -phi4.*phi2; phi4.*phi4; phi4.*phi5; ...
             -phi5.*phi1; -phi5.*phi2; phi5.*phi4; phi5.*phi5].*len_eI*w(l).*gamma./len_et1; 

    % \int {a dn u}[v] ds
    elmats2 = elmats2 + ...
             [ phi1.*dphin1;   phi1.*dphin2;  phi1.*dphin3;   phi1.*dphin4;  phi1.*dphin5;  phi1.*dphin6; ...
               phi2.*dphin1;   phi2.*dphin2;  phi2.*dphin3;   phi2.*dphin4;  phi2.*dphin5;  phi2.*dphin6; ...
               phi3.*dphin1;   phi3.*dphin2;  phi3.*dphin3;   phi3.*dphin4;  phi3.*dphin5;  phi3.*dphin6; ...
               -phi4.*dphin1; -phi4.*dphin2; -phi4.*dphin3;  -phi4.*dphin4; -phi4.*dphin5; -phi4.*dphin6; ...
               -phi5.*dphin1; -phi5.*dphin2; -phi5.*dphin3;  -phi5.*dphin4; -phi5.*dphin5; -phi5.*dphin6; ...
               -phi6.*dphin1; -phi6.*dphin2; -phi6.*dphin3;  -phi6.*dphin4; -phi6.*dphin5; -phi6.*dphin6] ...
                 .*len_eI*w(l)/2.*a; % /2 is coming from the average

end

ii = [tp1;tp1;tp1;tp1;
      tp2;tp2;tp2;tp2;
      tp4;tp4;tp4;tp4;
      tp5;tp5;tp5;tp5];
jj = [tp1;tp2;tp4;tp5];
jj = [jj;jj;jj;jj];
S1 = sparse(ii(:),jj(:),elmats1(:),np,np);

ii = [tp1;tp1;tp1;tp1;tp1;tp1;
      tp2;tp2;tp2;tp2;tp2;tp2;
      tp3;tp3;tp3;tp3;tp3;tp3;
      tp4;tp4;tp4;tp4;tp4;tp4;
      tp5;tp5;tp5;tp5;tp5;tp5;
      tp6;tp6;tp6;tp6;tp6;tp6 ];
jj = [tp1;tp2;tp3;tp4;tp5;tp6];
jj = [jj;jj;jj;jj;jj;jj];

S2 = sparse(ii(:),jj(:),elmats2(:),np,np);


end

