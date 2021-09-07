%%% Ce script permet la construction d'une figure compilant toutes les analyses stats relatives
%%% à la mesure de c extraites de Lipro3

% clear all;
% close all;
% vitesse=["0.0" "-0.2" "-0.4" "-0.6" "-0.8" "-1.0" "-1.2" "-1.4" "-1.6" "-2.0"];
% vitesse=["0.0" "0.2" "0.4" "0.6" "0.8" "1.0"];
% vitesse=["-1.0" "-0.8" "-0.6" "-0.4" "-0.2" "0.0" "0.2" "0.4" "0.6" "0.8" "1.0"]
vitesse=["-2.0" "-1.6" "-1.4" "-1.2" "-1.0" "-0.8" "-0.6" "-0.4" "-0.2" "0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0"];
% speed=[0.0 -0.2 -0.4 -0.6 -0.8 -1.0 -1.2 -1.4 -1.6 -2.0];
% speed=[0.0 0.2 0.4 0.6 0.8 1.0];
% speed=[-1.0 -0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 1.0]
speed=[-2.0 -1.6 -1.4 -1.2 -1.0 -0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];
N=2000;
rep=1000;
seuil=20;
graphique=[[] []];
j=1;
for v=vitesse
    v
    pause(0.0001)
%     eval(sprintf('cd /home/gascuel/PHD/Lipro/Lipro3/TriangleSup/V=%s', v));
    clear col row FrontVague nM Reg x y 
    if speed(j)>0 
        eval(sprintf('cd /home/gascuel/PHD/Lipro/Lipro3/Données+/V=%s', v));
    else
        eval(sprintf('cd /home/gascuel/PHD/Lipro/Lipro3/Données-/V=%s', v));
    end
    FrontVague=dlmread(sprintf("data_FrontVague_lipro_N_%g_V_%s00000_rep_%g.dat",N,v,rep));
    nM=dlmread(sprintf("data_nM_lipro_N_%g_V_%s00000_rep_%g.dat",N,v,rep)); 
    
    Reg=[];
    [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
    for i=1:size(row,1)-1
        x=nM(row(i):row(i+1)-1,1);
        y=nM(row(i):row(i+1)-1,2);
        if max(y)>seuil
            reg=x\y;
            Reg=[Reg reg];
        end
    end
    %analyse des propagations
    [row2,col2]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
    for i=2:size(row2,1)
        actifin(i-1)=nM(row2(i)-1,2)*100/N;
    end
    m2=mean(actifin);
    sd2=std(actifin);
    es2=sd2/sqrt(size(actifin,2));
    [row3,col3]=find(actifin'==100);
    if isempty(row3)==0
        s(j)=size(row3,1)/rep;
    else
        s(j)=0;
    end
%     [h,p]=kstest((Reg-mean(Reg))/std(Reg))
    graphique2(1,j)=m2;
    graphique2(2,j)=sd2;

    m1=mean(Reg);
    sd1=std(Reg);
    es1=sd1/sqrt(size(Reg,2));
    %[h,p]=kstest((Reg-mean(Reg))/std(Reg))
    graphique(1,j)=m1;
    graphique(2,j)=es1;
    j=j+1;
end

subplot(2,1,1);
errorbar(speed,graphique2(1,:),graphique2(2,:),'-x')

ylabel({"Nombre d'individus" "propagés en moyenne"});
title("Compilation des moyennes d'individus propagés");


subplot(2,1,2);
errorbar(speed,graphique(1,:),graphique(2,:))
xlabel('Vitesses de mouvement');
ylabel({"Vitesse moyenne" "de la vague d'activation"});
title('Compilation des moyennes de vitesses mesurées par ID');
% saveName = ('fig_cfonctionvtot_triangle'); 
% saveas(gcf, saveName, 'jpg'); 
% csvwrite('vitesse_moy_tri.dat',[speed' graphique']);

% % %%%%Regression linéaire
% % hold on;
% % b=[ones(length(speed(6:20)'),1) speed(6:20)']\graphique(1,6:20)';
% % ycalc=[ones(length(speed'),1) speed']*b;
% % plot(speed,ycalc)
% % saveName = ('fig_cfonctionvtot_résumé+fit'); 
% % legend('Data','Reg')
% % saveas(gcf, saveName, 'jpg');



