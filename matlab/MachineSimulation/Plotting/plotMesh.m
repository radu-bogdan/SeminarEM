function [h] = plotMesh(p,t,e,material)
    np=size(p,2); nt=size(t,2);
    h = pdeplot(p,e,t);
    h(2).Color = 'none';
    h(1).Color = 'none';
    IronPos = ismember(t(4,:),material.iron);
    MagnetPos = ismember(t(4,:),material.magnet);
    CoilPos = ismember(t(4,:),material.coil);
    
    str1 = '#FBCB00';
    color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;
    str2 = '#AEAEAE';
    color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;
    str3 = '#DADCDC';
    color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;
    
    patch('faces',t(1:3,IronPos)','vertices',p','facecolor',color2,'edgecolor',color2)
    patch('faces',t(1:3,MagnetPos)','vertices',p','facecolor','red','edgecolor','red')
    patch('faces',t(1:3,CoilPos)','vertices',p','facecolor','green','edgecolor','green')
    xlabel('x [m]') 
    ylabel('y [m]')
    %title(strcat('numper of points =  ',num2str(np),', number of triangles =  ',num2str(nt)))
end

