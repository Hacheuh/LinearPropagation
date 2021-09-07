%%% Ce script permet la construction de toutes les analyses stats relatives
%%% aux quantité extraites de Lipro3
%cd '/home/gascuel/PHD/Lipro/Lipro3/Données-/V=-1.0'
%cd '/home/gascuel/Desktop/Supersonic'
cd  /tmp/wkm5/DoubleS
clear all;
close all;

N=2000;
rep=20;
v=1000;
seuil=20;

% Finalid=dlmread(sprintf("data_Finalid_lipro_N_%g_V_%g0000.dat",N,v));
% % FrontVague=dlmread(sprintf("data_FrontVague_lipro_N_%g_V_%g0000.dat",N,v));
% % nM=dlmread(sprintf("data_nM_lipro_N_%g_V_%g0000.dat",N,v));
FrontVague=dlmread(sprintf("data_FrontVague_lipro_N_%g_V_%g.000000_rep_%g_vc.dat",N,v,rep));
nM=dlmread(sprintf("data_nM_lipro_N_%g_V_%g.000000_rep_%g_vc.dat",N,v,rep));
FrontVagueL=dlmread(sprintf("data_FrontVagueLeft_lipro_N_%g_V_%g.000000_rep_%g_vc.dat",N,v,rep));
%Delta=dlmread(sprintf("data_Delta_lipro_N_%g_V_%g0000.dat",N,v));

% %Distributions des arrêt : Finalid
% figure;
% hist(Finalid,40);
% xlabel('Longueur de chaîne');
% ylabel("Nombre d'occurence");
% title("Histogramme de la longueur de chaîne de transmission d'information"); 
% saveName = (['fig_Finalid_v_', num2str(v)]); 
% saveas(gcf, saveName, 'jpg'); 

%Distributions des deltas : Deltas
% figure;
% hist(Delta(:,2),20);
% xlabel('Ecart entre le plus récent en mouvement et le prochain statique');
% ylabel("Nombre d'occurence")
% title('Histogramme des deltas')
% saveName = (['fig_Delta_v_', num2str(v)]); 
% saveas(gcf, saveName, 'jpg'); 

%Distribution des pentes : FrontVague, nM
figure;
hold on;
Reg=[];
[row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
for i=1:size(row,1)-1
    x=nM(row(i):row(i+1)-1,1);
    y=nM(row(i):row(i+1)-1,2);
    if max(y)>seuil
        reg=x\y;
        Reg=[Reg reg];
        yCalc=reg*x;
        %plot(x,yCalc);%plot des regressions linéaires
        plot(x,y);%plot des raw data
    end
end
xlabel('Time');
ylabel("Number of active individuals");
title('Progression of the state Active for several simulations (nM)');
% saveName = (['fig_nMCompil_v_', num2str(v)]); 
% saveas(gcf, saveName, 'png'); 

figure;
histogram(Reg,25);
xlabel("Vitesse de transmission de l'information");
ylabel("Nombre d'occurence");
title('Histogramme des regressions linéaires pour nM (Vitesses)');
% saveName = (['fig_nMHistLin_v_', num2str(v)]); 
% saveas(gcf, saveName, 'jpg'); 
m1=mean(Reg)
sd1=std(Reg)
[h,p]=kstest((Reg-mean(Reg))/std(Reg))

%FrontVague
figure;
hold on;
Reg=[];
[row,col]=find(FrontVague(:,1)==0); %on récupère toute les lignes de début de simu
for i=1:size(row,1)-1
    x=FrontVague(row(i):row(i+1)-1,1);
    y=FrontVague(row(i):row(i+1)-1,2);
    if max(y)>seuil
        reg=x\y;
        Reg=[Reg reg];
        yCalc=reg*x;
        plot(x,yCalc);
    end
end
xlabel('Temps');
ylabel("Indicatif de l'individus le plus proche du front mobile");
title('Compilation des courbes de progression des états pour chaque simulation FV');
% saveName = (['fig_FrontVCompil_v_', num2str(v)]); 
% saveas(gcf, saveName, 'jpg'); 

figure;
histogram(Reg,25);
xlabel("Vitesse de transmission de l'information");
ylabel("Nombre d'occurence");
title('Histogramme des regressions linéaires pour FrontV (Vitesses)');
% saveName = (['fig_FrontVHistLin_v_', num2str(v)]); 
% saveas(gcf, saveName, 'jpg'); 
m2=mean(Reg)
sd2=std(Reg)
[h,p]=kstest((Reg-mean(Reg))/std(Reg))


%FrontVagueLeft
figure;
hold on;
Reg=[];
[row,col]=find(FrontVagueL(:,1)==0); %on récupère toute les lignes de début de simu
for i=1:size(row,1)-1
    x=FrontVagueL(row(i):row(i+1)-1,1);
    y=FrontVague(row(i):row(i+1)-1,2);
    if max(y)>seuil
        reg=x\y;
        Reg=[Reg reg];
        yCalc=reg*x;
        plot(x,yCalc);
    end
end
xlabel('Temps');
ylabel("Indicatif de l'individus le plus proche du front mobile de droite");
title('Compilation des courbes de progression des états pour chaque simulation FV');
% saveName = (['fig_FrontVLCompil_v_', num2str(v)]); 
% saveas(gcf, saveName, 'jpg'); 


figure;
histogram(Reg,25);
xlabel("Vitesse de transmission de l'information");
ylabel("Nombre d'occurence");
title('Histogramme des regressions linéaires pour FrontVL (Vitesses)');
% saveName = (['fig_FrontVLHistLin_v_', num2str(v)]); 
% saveas(gcf, saveName, 'jpg'); 
m3=mean(Reg)
sd3=std(Reg)
[h,p]=kstest((Reg-mean(Reg))/std(Reg))

%%% Pour un fichier de taille rep=1
% plot([FrontVague(:,1) FrontVagueL(:,1) nM(:,1)], [FrontVague(:,2) FrontVagueL(:,2) nM(:,2)])
% xlabel('Temps');
% ylabel("ID front/ Population");
% title({"Superposition des 3 courbes de mesures:" "ID du front à droite, à gauche et population d'actif"});
% legend({"Front droit" "Front gauche" "Pop"}, "Location", "northwest")
% saveName = (['fig_Sup3mesuresLRP_v_', num2str(v)]); 
% saveas(gcf, saveName, 'png');  

%%%% Faire un fichier des instants d'evenements

% j=0;
% for i=1:size(nM,1)
%     if (i+1~=size(nM,1)+1)
%         if (nM(i,2)~=nM(i+1,2))
%             j=j+1;
%             Ev(j)=nM(i+1,1);
%         end
%     end
% end
% 
% %%% Calcul des différences 
% for i=1:size(Ev,2)
%     if i+1~=size(Ev,2)+1
%         if Ev(i+1)-Ev(i)>=0
%             DiffEve(i)=Ev(i+1)-Ev(i);
%         end
%     end
% end
% 
% histogram(DiffEve);
% 
% fid=fopen('diffeve.dat','w');
% for i=1:size(DiffEve,2)
%     fprintf(fid,"%g\t",DiffEve(i));
%     fprintf(fid,"\n");
% end
% fclose(fid);

