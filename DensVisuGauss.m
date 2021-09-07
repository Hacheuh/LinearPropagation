%%% Ce script permet la constitution d'une visualisation en densité dynamique
%%% à base de gaussiennes à partir du fichier movers fourni par Lipro3
cd /home/gascuel/Desktop/
clear all;
close all;
N=7825;
rep=1;
v=-3.0;
pde=0.001; %pas d'espace
sd=0.2; %ecart-type de la gaussienne
Movers=dlmread(sprintf("data_Movers_lipro_N_%g_V_%g.000000_rep_%g.dat",N,v,rep));

% writerObj = VideoWriter(sprintf('LiproVisu_N%g_V%g_pde%g_Gtest.avi',N,v,pde));
% open(writerObj);
%
% h1=figure;
% for i=1:N %remplissage de la condition initiale
%     pos(i)=i-1;
% end
% x=[pos(1)-1:pde:pos(N)+1]; % constitution de l'espace de visu
% 
% [row,col]=find(Movers(:,1)==-0.023); %on trouve les débuts de simu
% for i=1:size(Movers,1)-1%row(2)-1 si on a plusieurs simu
%     col2=find(Movers(i,:),1,'last'); %comme ce n'est pas une matrice NxN on a des zeros artefactuels, on cherche l'id du dernier élément avant un 0
%     for j=1:col2  %on actualise les positions
%         pos(j)=Movers(i,j);
%     end
%     sum=zeros(1,(N+1)/pde+1);%on initialise la somme pour ce pas de temps
%     for j=1:N
%         norm=normpdf(x,pos(j),sd);%on place les gaussiennes sur chaque indiv
%         for k=pos(j)-1:pde:pos(j)+1% à chaque gaussienne créée on edite la somme autour de la position
%             id=int32((k+1)/pde+1);
%             sum(id)=sum(id)+norm(id);
%         end
%     end
%     figure(h1);
%     plot(x,sum);
%     axis([pos(col2)-10 pos(col2)+10 -1 5]);% metrte pos(1) si on suit tous les movers
%     F = getframe(gcf);
%     writeVideo(writerObj,F);
% end
% close(writerObj);   


%%%%% Visualisation boule 1D
% %writerObj = VideoWriter(sprintf('LiproVisu_demo-.avi',N,v,pde));
% %open(writerObj);
% h1=figure;
% for i=1:N %remplissage de la condition initiale
%     pos(i)=i-1;
%     VectColor(i,1)= 1;VectColor(i,2)= 0;VectColor(i,3)= 0;
% end
% Z=zeros(size(pos,2),1);
% [row,col]=find(Movers(:,1)==-0.0131); %on trouve les débuts de simu
% for i=1:200%size(Movers,1)-1%row(2)-1 si on a plusieurs simu
%     col2=find(Movers(i,:),1,'last'); %comme ce n'est pas une matrice NxN on a des zeros artefactuels, on cherche l'id du dernier élément avant un 0
%     for j=1:col2  %on actualise les positions
%         pos(j)=Movers(i,j);
%         VectColor(j,1)= 0;VectColor(j,2)= 1;VectColor(j,3)= 0;
%     end
%     
%     figure(h1);
%     scatter(pos,Z,100,VectColor,'Filled');
%     axis([-1 10 -1 1]);% metrte pos(1) si on suit tous les movers
%     %F = getframe(gcf);
%     %writeVideo(writerObj,F);
% end
% %close(writerObj);   

% % %%%%% Visualisation boule 2D : Attention les fichiers movers utilisés
% %%%% ont été fait sur mesure
%writerObj = VideoWriter(sprintf('LiproVisu_demo2D_Disq.avi',N,v,pde));
%open(writerObj);
h1=figure;

for i=0:(size(Movers,1)/N)-1%size(Movers,1)-1%row(2)-1 si on a plusieurs simu
    
    for j=1:N  %on actualise les positions
        if Movers(i*N+j,3)==0
            VectColor(j,1)= 0;VectColor(j,2)= 1;VectColor(j,3)= 0;
        else
            VectColor(j,1)= 1;VectColor(j,2)= 0;VectColor(j,3)= 0;
        end
    end
    
    figure(h1);
    scatter(Movers(i*N+1:(i+1)*N,1),Movers(i*N+1:(i+1)*N,2),10,VectColor,'Filled');
    axis([-50 50 0 100]);% metrte pos(1) si on suit tous les movers
    %axis([-10 110 -10 110]);
    %F = getframe(gcf);
    %writeVideo(writerObj,F);
end

%close(writerObj);   


%%%%%%% Chemographe 1D
% h1=figure;
% hold on;
% for i=1:N %remplissage de la condition initiale
%     pos(i)=i-1;
%     VectColor(i,1)= 1;VectColor(i,2)= 0;VectColor(i,3)= 0;
% end
% 
% for i=1:1000%size(Movers,1)-1%row(2)-1 si on a plusieurs simu
%     Z=zeros(size(pos,2),1)+i*0.1/1.53;
%     col2=find(Movers(i,:),1,'last'); %comme ce n'est pas une matrice NxN on a des zeros artefactuels, on cherche l'id du dernier élément avant un 0
%     for j=1:col2  %on actualise les positions
%         pos(j)=Movers(i,j);
%         VectColor(j,1)= 0;VectColor(j,2)= 1;VectColor(j,3)= 0;
%     end
%     
%     figure(h1);
%     scatter(pos,Z,20,VectColor,'s','Filled');
%     axis([-10 250 -5 1000*0.1/1.53+10]);% metrte pos(1) si on suit tous les movers
% end
% saveName = ('ChemoTestSq20_T1000'); 
% saveas(gcf, saveName, 'png'); 

%%chemographe double sens
% h1=figure;
% hold on;
% for i=1:N %remplissage de la condition initiale
%     pos(i)=i-1;
%     VectColor(i,1)= 1;VectColor(i,2)= 0;VectColor(i,3)= 0;
% end
% Tmp=zeros(N,1);
% for i=1:1000%size(Movers,1)-1%row(2)-1 si on a plusieurs simu
%     Z=zeros(size(pos,2),1)+i*0.1/1.53;
%     col2=find(Movers(i,:),1,'last'); %comme ce n'est pas une matrice NxN on a des zeros artefactuels, on cherche l'id du dernier élément avant un 0
%     if isempty(col2)==0
%         Tmp(i)=Movers(i,col2);
%     end
%     for j=1:col2  %on actualise les positions
%         pos(round(Tmp(j),0))=Movers(i,j);
%         VectColor(round(Tmp(j),0),1)= 0;VectColor(round(Tmp(j),0),2)= 1;VectColor(round(Tmp(j),0),3)= 0;
%     end
%     
%     figure(h1);
%     scatter(pos,Z,3,VectColor,'s','Filled');
%     axis([500 1500 -1 1000*0.1/1.53+10]);% metrte pos(1) si on suit tous les movers
% end
% saveName = ('ChemoTestDoubleSq1_T1000'); 
% saveas(gcf, saveName, 'png'); 