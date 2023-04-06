function eI = calcEdgeInterfaceCircle(p,eI1,eI2,t)
% see calcEdgeInterfaceCircle_nb.m for a more detailed description

neI1 = size(eI1,2); neI2 = size(eI2,2);

% boundary between bodies
eI = zeros(10,(neI1+neI2)*2);

e1 = eI1(1,:); e2 = eI1(2,:);
% e3 and e4 swapped because the edges are equally orientated
e4 = eI2(1,:); e3 = eI2(2,:);

p1 = p(:,e1); p2 = p(:,e2);
p3 = p(:,e3); p4 = p(:,e4);

% calculate angle of points
p1theta = cart2pol(p1(1,:),p1(2,:));
p2theta = cart2pol(p2(1,:),p2(2,:));
p3theta = cart2pol(p3(1,:),p3(2,:));
p4theta = cart2pol(p4(1,:),p4(2,:));

% There is a problem if one of the two points of an edge has an angle close
% to -pi and the other close to pi. The following code resolves that.
for i = 1:2
    if i == 1
        pi_term = 2*pi;
    else
        pi_term = pi/2;
    end
    recalc = abs(p1theta - p2theta) > pi | abs(p3theta - p4theta) > pi;
    p1theta(recalc) = mod(p1theta(recalc)+pi_term,2*pi);
    p2theta(recalc) = mod(p2theta(recalc)+pi_term,2*pi);
    p3theta(recalc) = mod(p3theta(recalc)+pi_term,2*pi);
    p4theta(recalc) = mod(p4theta(recalc)+pi_term,2*pi);
end

% checking if edges is intersecting
k = 1; % index for eI
for i = 1:neI1
    
    for j = 1:neI2

        lambda4 = (p4theta(j) - p1theta(i))/(p2theta(i) - p1theta(i));
        lambda3 = (p3theta(j) - p1theta(i))/(p2theta(i) - p1theta(i));

        intersect = 0;
        if lambda3 > 0
            if lambda4 < 0
                intersect = 1;
                eI(1,k) = e1(i);
                pI1 = p1theta(i);
                if lambda3 > 1
                    eI(2,k) = e2(i);
                    pI2 = p2theta(i);
                else
                    eI(2,k) = e3(j);
                    pI2 = p3theta(j);
                end
            elseif lambda4 < 1
                intersect = 1;
                eI(1,k) = e4(j);
                pI1 = p4theta(j);
                if lambda3 > 1
                    eI(2,k) = e2(i);
                    pI2 = p2theta(i);
                else
                    eI(2,k) = e3(j);
                    pI2 = p3theta(j);
                end
            end
        end
        
        if intersect 
            lambda1 = (pI1 - p1theta(i))/(p2theta(i) - p1theta(i));
            lambda2 = (pI2 - p1theta(i))/(p2theta(i) - p1theta(i));
            lambda3 = (pI1 - p3theta(j))/(p4theta(j) - p3theta(j));
            lambda4 = (pI2 - p3theta(j))/(p4theta(j) - p3theta(j));
            eI(7:10,k) = [lambda1; lambda2; lambda3; lambda4];

            tInd1 = eI1(4,i); tInd2 = eI2(4,j);
            eI(3:4,k) = [tInd1; tInd2];
    
            [eI(5,k),~] = find(t(1:3,tInd1)==e1(i));
            [eI(6,k),~] = find(t(1:3,tInd2)==e3(j));
    
            % if the two new edge points are too close, do not consider it
            % as an element
            if abs(lambda1-lambda2) > 1e-15
                k = k+1;
            end
        end
    end
end

eI = eI(:,1:k-1);

