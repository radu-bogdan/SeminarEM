function [p1,e1,t1,u1]=refinemesh(p0,e0,t0,u0)
% [p,e,t,u]=refinemesh(p,e,t,u)
%    uniform refinement of mesh
%    if 4 arguments are given, then the function 
%    u given as np x 1 vector is interpolated 
%    to finer mesh as well

np0=size(p0,2); nt0=size(t0,2); ne0=size(e0,2);

if size(t0,1)<4, t0=[t0;1:nt0]; end
if size(e0,1)<3, e0=[e0;1:ne0]; end

% create new points at adge midpoints
[ue,te,be]=makeedgelist(t0,e0);
pn = 0.5*(p0(:,ue(1,:)) + p0(:,ue(2,:)));
p1 = [p0,pn];
  
% create elements
te = te + np0;
t1  = [t0(1,:), t0(2,:), t0(3,:), te(1,:); ...
       te(1,:), te(2,:), te(3,:), te(2,:); ...
       te(3,:), te(1,:), te(2,:), te(3,:); ...
       t0(4,:), t0(4,:), t0(4,:), t0(4,:)];

% create boundary elements
be = be + np0;
e1 = [e0(1,:),   be; ...
      be,        e0(2,:); ...
      e0(end,:), e0(end,:);];

if nargout>3 
   u1 = [u0;0.5*(u0(ue(1,:),:)+u0(ue(2,:),:))]; 
end


