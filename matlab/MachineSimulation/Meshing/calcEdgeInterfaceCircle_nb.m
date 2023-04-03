function eI = calcEdgeInterfaceCircle_nb(p,eI1,eI2,t)
% Notebook version
% calculates the partial edges on an interface.

% structure of eI in the end
% eI(1,:) ...   p-nr 1 
% eI(2,:) ...   p-nr 2 
% eI(3,:) ...   triangle at domain 1 belonging to the edge
% eI(4,:) ...   triangle at domain 2 belonging to the edge
% eI(5,:) ...   point of the triangle eI(3,:) such that 
%               the edge to the next point of the triangle in 
%               counterclockwise orientation is a superset of the interface 
%               edge.
% eI(6,:) ...   point of the triangle eI(4,:) such that 
%               the edge to the next point of the triangle in 
%               counterclockwise orientation is a superset of the interface 
%               edge.
% eI(7,:) ...   relative position of p(:,eI(1,:)) on the edge p1 - p2
% eI(8,:) ...   relative position of p(:,eI(2,:)) on the edge p1 - p2
% eI(9,:) ...   relative position of p(:,eI(1,:)) on the edge p3 - p4
% eI(10,:) ...  relative position of p(:,eI(2,:)) on the edge p3 - p4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- extremely inefficient version!!!! ----%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neI1 = size(eI1,2); neI2 = size(eI2,2);


% boundary between domains
eI = zeros(10,(neI1+neI2)*2);

k = 1; % index for eI
for i = 1:neI1
    e1 = eI1(1,i); e2 = eI1(2,i); 
    p1 = p(:,e1); p2 = p(:,e2); 
    
    for j = 1:neI2
        % e3 and e4 swapped because the edges are equally orientated
        e4 = eI2(1,j); e3 = eI2(2,j); 
        p3 = p(:,e3); p4 = p(:,e4);

        ptheta = cart2pol([p1(1),p2(1),p3(1),p4(1)], [p1(2),p2(2),p3(2),p4(2)]);

        if abs(ptheta(1) - ptheta(2)) > pi || abs(ptheta(3) - ptheta(4)) > pi
            ptheta = mod(ptheta+2*pi,2*pi);
            if abs(ptheta(1) - ptheta(2)) > pi || abs(ptheta(3) - ptheta(4)) > pi
                ptheta = mod(ptheta+pi/2,2*pi);
            end
        end
        lambda4 = (ptheta(4) - ptheta(1))/(ptheta(2) - ptheta(1));
        lambda3 = (ptheta(3) - ptheta(1))/(ptheta(2) - ptheta(1));

        intersect = 0;
        if lambda3 > 0
            if lambda4 < 0
                intersect = 1;
                eI(1,k) = e1;
                pI1 = ptheta(1);
                if lambda3 > 1
                    eI(2,k) = e2;
                    pI2 = ptheta(2);
                else
                    eI(2,k) = e3;
                    pI2 = ptheta(3);
                end
            elseif lambda4 < 1
                intersect = 1;
                eI(1,k) = e4;
                pI1 = ptheta(4);
                if lambda3 > 1
                    eI(2,k) = e2;
                    pI2 = ptheta(2);
                else
                    eI(2,k) = e3;
                    pI2 = ptheta(3);
                end
            end
        end
        
        if intersect 
            lambda1 = (pI1 - ptheta(1))/(ptheta(2) - ptheta(1));
            lambda2 = (pI2 - ptheta(1))/(ptheta(2) - ptheta(1));
            lambda3 = (pI1 - ptheta(3))/(ptheta(4) - ptheta(3));
            lambda4 = (pI2 - ptheta(3))/(ptheta(4) - ptheta(3));
            eI(7:10,k) = [lambda1; lambda2; lambda3; lambda4];

            tInd1 = eI1(4,i); tInd2 = eI2(4,j);
            eI(3:4,k) = [tInd1; tInd2];
    
            [eI(5,k),~] = find(t(1:3,tInd1)==e1);
            [eI(6,k),~] = find(t(1:3,tInd2)==e3);
    
            % if the two new edge points are too close, do not consider it
            % as an element
            if abs(lambda1-lambda2) > 1e-14
                k = k+1;
            end
        end
    end
end

eI = eI(:,1:k-1);

