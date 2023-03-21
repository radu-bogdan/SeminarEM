clear;
load -mat motor.mat
m = m/((10^7)/(4*pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate impressed current density
nt=size(t,2);
%Calculate indices
coilArrayInd  = strfind(regions_2d,'coil');
coilArray = find(not(cellfun('isempty',coilArrayInd))); %contains the t(4,:) values which corresponds to coils
jInd = zeros(nt,1);
for i = 1:1:nt
    var = char(regions_2d(t(4,i)));
    if ismember(t(4,i),coilArray)
        jInd(i) = sscanf(var,'coil%d');
        xd = t(4,i);
    else
        jInd(i) = 1;
    end
end
J3_vec = ismember(t(4,:),coilArray).*j3(int64(jInd));

save('J3_vec.mat','J3_vec')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate magnetiziation
magArrayInd  = strfind(regions_2d,'mag');
magArray = find(not(cellfun('isempty',magArrayInd)));   %contains the t(4,:) values which corresponds to magnets
mInd = zeros(nt,1);
for i = 1:1:nt
    var = char(regions_2d(t(4,i)));
    if ismember(t(4,i),magArray)
        mInd(i) = sscanf(var,'magnet%d');
    else
        mInd(i) = 1;
    end
end
M_vec = [ismember(t(4,:),magArray);ismember(t(4,:),magArray)].*m(:,int64(mInd));
save('M_vec.mat','M_vec')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create object with lists for every different material
material.coil = coilArray;
material.magnet = magArray;
airArrayInd  = strfind(regions_2d,'air');
airArray = find(not(cellfun('isempty',airArrayInd)));
material.air = airArray;
material.iron = [146,150];
material.non_conductive_shaft = [145];
material.all_non_conductive = [material.coil,material.magnet,material.air,material.non_conductive_shaft];

material.k1 = 49.4;
material.k2 = 1.46;
material.k3 = 520.6;
material.nu0 = (10^7)/(4*pi);
save('material.mat','material')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find boundary
edgeInd  = strfind(regions_1d,'stator_outer');
edgeArray = find(not(cellfun('isempty',edgeInd)));
isMem = ismember(e(3,:),edgeArray);
stator_outer_Indices = find(isMem);
e_stator_outer = e(:,stator_outer_Indices);
save('e_stator_outer.mat','e_stator_outer')

