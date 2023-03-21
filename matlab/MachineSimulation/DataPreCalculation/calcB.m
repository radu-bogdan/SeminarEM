function [B] = calcB(p,t,e)
%volume integral
t1=t(1,:); t2=t(2,:); t3=t(3,:);
B.B11=p(1,t2)-p(1,t1);
B.B12=p(1,t3)-p(1,t1);
B.B21=p(2,t2)-p(2,t1);
B.B22=p(2,t3)-p(2,t1);
B.detB = B.B11.*B.B22-B.B12.*B.B21;
B.iB11= B.B22./B.detB;
B.iB12=-B.B12./B.detB; 
B.iB21=-B.B21./B.detB;
B.iB22= B.B11./B.detB; 
B.vol = abs(B.detB);

%edge integrals
e1=e(1,:);  e2=e(2,:);
p1=p(:,e1); p2=p(:,e2);
B11e = p2(1,:)-p1(1,:); 
B21e = p2(2,:)-p1(2,:);
B.len = sqrt(B11e.^2+B21e.^2);

dphih1=[-1;-1]; dphih2=[1;0]; dphih3=[0;1];
B.iB11= B.B22./B.detB; B.iB12=-B.B12./B.detB; 
B.iB21=-B.B21./B.detB; B.iB22= B.B11./B.detB; 
B.dphi1x= B.iB11*dphih1(1) + B.iB21*dphih1(2);   B.dphi1y= B.iB12*dphih1(1)+ B.iB22*dphih1(2);
B.dphi2x= B.iB11*dphih2(1) + B.iB21*dphih2(2);   B.dphi2y= B.iB12*dphih2(1)+ B.iB22*dphih2(2);
B.dphi3x= B.iB11*dphih3(1) + B.iB21*dphih3(2);   B.dphi3y= B.iB12*dphih3(1)+ B.iB22*dphih3(2);
end

