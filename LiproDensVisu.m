%%% Ce script permet la constitution d'une visualisation en densité dynamique à
%%% partir du fichier movers fourni par Lipro3 
cd /home/gascuel/Desktop
clear all;
close all;
N=2000;
rep=1;
v=-0.23;
pde=0.01; %pas d'espace
writerObj = VideoWriter(sprintf('LiproVisu_N%g_V%g_pde%g.avi',N,v,pde));
open(writerObj);

Movers=dlmread(sprintf("data_Movers_lipro_N_%g_V_%g0000.dat",N,v));
for i=1:N %remplissage de la condition initiale
    pos(i)=i-1;
end

h1=figure;

[row,col]=find(Movers(:,1)==-0.023); %on trouve les débuts de simu
for i=1:size(Movers,1)-1%row(2)-1 si on a plusieurs simu
    [row2,col2]=find(Movers(i,:)==0); %comme ce n'est pas une matrice NxN on a des zeros artefactuels
    for j=1:col2(1)-1 %on actualise les positions
        pos(j)=Movers(i,j);
    end
    hist(pos,N-round(pos(1)));
    figure(h1);
    axis([pos(1)-10 pos(col2(1))+10 -1 5]);
    F = getframe(gcf);
    writeVideo(writerObj,F);
    %pause(0.001);
end
close(writerObj);        