function machine_geo = createBoschMachineGeo()
% Creates a decomposed geometry matrix (see "doc decsg" for more 
% information) which represents an electric machine from Bosch 
% (data obtained from Alessio). 


origin = [0,0];
%     #inner radius rotor
    r1 = 26.5*10^(-3);
%     #outer radius rotor
    r2 = 78.63225*10^(-3);
%     #sliding mesh rotor
    r4 = 78.8354999*10^(-3);
%     #sliding mesh stator
    r6 = 79.03874999*10^(-3);
%     #inner radius stator
    r7 = 79.242*10^(-3);
%     #outer radius stator
    r8 = 116*10^(-3);

%Points for magnet1 and air around magnet1
m1 = [69.23112999*10^(-3),7.535512*10^(-3)];
m2 = [74.828958945*10^(-3),10.830092744*10^(-3)];
m3 = [66.13621099700001*10^(-3),25.599935335*10^(-3)];
m4 = [60.53713*10^(-3),22.30748*10^(-3)];
a5 = [69.75636*10^(-3),5.749913*10^(-3)];
a6 = [75.06735*10^(-3),3.810523*10^(-3)];
a7 = [65.3506200*10^(-3),26.51379*10^(-3)];
a8 = [59.942145092*10^(-3),24.083661604*10^(-3)];
%Points for magnet2 and air around magnet2
m5 = [58.579985516*10^(-3), 27.032444757*10^(-3)];
m6 = [64.867251151*10^(-3),28.663475405*10^(-3)];
m7 = [60.570096319*10^(-3),45.254032279*10^(-3)];
m8 = [54.282213127*10^(-3),43.625389857*10^(-3)];
a1 = [53.39099766*10^(-3),45.259392713*10^(-3)];
a2 = [55.775078884*10^(-3),50.386185578*10^(-3)];
a3 = [59.41521771*10^(-3),25.355776837*10^(-3)];
a4 = [65.12210917100001*10^(-3),27.707477175*10^(-3)];
%Points for Stator Nut and air in the stator
s1 = [0.079143298919643,0.00395383359742138];
s2 = [80.143057128*10^(-3),4.0037794254*10^(-3)];
s3 = [80.387321219*10^(-3),2.965459706*10^(-3)];
s4 = [98.78501315600001*10^(-3),3.9007973292*10^(-3)];
s5 = [98.44904989600001*10^(-3),9.026606148400001*10^(-3)];
s6 = [80.086666706*10^(-3),7.5525611543*10^(-3)];
s7 = [79.980020247*10^(-3),6.4912415424*10^(-3)];
s8 = [0.078982295870193,0.00641026544483722];
% points for circle arcs on the inner side of the stator
s9 = [r7,0];
s10 = [cos(2*pi/48),sin(2*pi/48)]*r7;

% --------------- magnets --------------
tz = zeros(1,3); %three filling zero entries
magnets = [ % magnet 1
    2, m1, m2, 97, 99, tz; % start, end point
    2, m2, m3, 97, 146, tz;
    2, m3, m4, 97, 100 tz;
    2, m4, m1, 97, 146, tz;
    2, m1, a5, 99, 146, tz;
    2, a5, a6, 99, 146, tz;
    2, a6, m2, 99, 146, tz;
    2, m3, a7, 100, 146, tz;
    2, a7, a8, 100, 146, tz;
    2, a8, m4, 100, 146, tz;
    % magnet 2
    2, m5, m6, 98, 101, tz; 
    2, m6, m7, 98, 146, tz;
    2, m7, m8, 98, 102 tz;
    2, m8, m5, 98, 146, tz;
    2, m5, a3, 101, 146, tz;
    2, a3, a4, 101, 146, tz;
    2, a4, m6, 101, 146, tz;
    2, m7, a2, 102, 146, tz;
    2, a2, a1, 102, 146, tz;
    2, a1, m8, 102, 146, tz;
]; 

% adding rotated parts...
magnets_geo = magnets;
n_parts = 8;
for i = 1:(n_parts-1)
    theta = 2*pi/n_parts * i;
    Rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    magnets_rot = magnets;
    magnets_rot(:,[2,3]) = magnets_rot(:,[2,3])*Rot';
    magnets_rot(:,[4,5]) = magnets_rot(:,[4,5])*Rot';
    
    % adapting subdomain indices
    magnets_rot(:,6) = magnets_rot(:,6)+6*i;
    magnets_rot(magnets_rot(:,7)~=146,7) = magnets_rot(magnets_rot(:,7)~=146,7)+6*i;

    magnets_geo = [magnets_geo; magnets_rot];
end


% --------------- stator nuts --------------
statorNut = [ %coil6
    2, s2, s3, 1, 150, tz;
    2, s3, s4, 1, 150, tz;
    2, s4, s5, 1, 150, tz;
    2, s5, s6, 1, 150, tz;
    2, s6, s7, 1, 150, tz;
    2, s7, s2, 1, 2, tz;
    % air
    1, s9, s1, 149, 150, origin, r7;
    2, s1, s2, 2, 150, tz;
    2, s7, s8, 2, 150, tz;
    1, s1, s8, 149, 2, origin, r7;
    1, s8, s10,149, 150, origin, r7;
    ];

% adding rotated parts...
statorNut_geo = statorNut;
n_parts = 48;
for i = 1:(n_parts-1)
    theta = 2*pi/n_parts * i;
    Rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    statorNut_rot = statorNut;
    statorNut_rot(:,[2,3]) = statorNut_rot(:,[2,3])*Rot';
    statorNut_rot(:,[4,5]) = statorNut_rot(:,[4,5])*Rot';
    
    % adapting subdomain indices
    statorNut_rot(statorNut_rot(:,6)==1,6) = 1+2*i;
    statorNut_rot(statorNut_rot(:,6)==2,6) = 2+2*i;
    statorNut_rot(statorNut_rot(:,7)==2,7) = statorNut_rot(statorNut_rot(:,7)==2,7)+2*i;

    statorNut_geo = [statorNut_geo; statorNut_rot];
end

% --------------- additional circles --------------
circles_geo = [ %coil
    1, r1, 0, -r1, 0, 145, 146, origin, r1;
    1, -r1, 0, r1, 0, 145, 146, origin, r1;
    1, r2, 0, -r2, 0, 146, 147, origin, r2;
    1, -r2, 0, r2, 0, 146, 147, origin, r2;    
    1, r4, 0, -r4, 0, 147, 148, origin, r4;
    1, -r4, 0, r4, 0, 147, 148, origin, r4;
    1, r6, 0, -r6, 0, 148, 149, origin, r6;
    1, -r6, 0, r6, 0, 148, 149, origin, r6; 
    1, r8, 0, -r8, 0, 150, 0,   origin, r8;
    1, -r8, 0, r8, 0, 150, 0,   origin, r8;
    ];

machine_geo = [magnets_geo; statorNut_geo; circles_geo];
machine_geo(:,[3,4]) = machine_geo(:,[4,3]);
machine_geo = machine_geo';

end

