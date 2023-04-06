function [ued,te,bed]=makeedgelist(t,e)
% [te,ue,bed] = makeedgelist
%    make list of all edges contained in mesh
%    input: 
%      t ... 3 x nt list of triangle point numbers
%      e ... 2 x ne list of boundary edge point numbers
%    output 
%      ued ... list of edge point numbers
%      te  ... list of triangle edge numbers
%      bed ... list of boundary edge numbers  

%% list of all edges
ed=[t(1,:),t(2,:),t(3,:);t(2,:),t(3,:),t(1,:)];
sed=sort(ed);
[ued,ie,je]=unique(sed','rows'); 
ued=ued';

%% triangle to edge list
ned=size(ed,2);
te=reshape(je,ned/3,3)';
  
%% find boundary edges
ed=sort(e(1:2,:));
[b,i,j]=unique([ed';ued'],'rows');
bed=j(1:size(ed,2))'; 

