% close all
% % load('machine_geo.mat')
% machine_geo = createBoschMachineGeo();
% model = createpde;
% geometryFromEdges(model,machine_geo);
% figure;
% pdegplot(model,"FaceLabels","on")
% mesh = generateMesh(model,"GeometricOrder",'linear');
% figure;
% 
% pdeplot(mesh,"ElementLabels","on")
% 
% p = mesh.Nodes;
% t = mesh.Elements(:,1:100);
% figure;
% mypdemesh(p,t)


!python felismesh.py 2 1
load -mat motor.mat

size(t)