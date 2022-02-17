%%% Ce script permet la construction d'une figure compilant toutes les analyses stats relatives
%%% à la mesure de c extraites de Lipro3

clear all;
close all;
%vitesse=["0.0" "-0.2" "-0.4" "-0.6" "-0.8" "-1.0" "-1.2" "-1.4" "-1.6" "-2.0"];
% vitesse=["0.0" "0.2" "0.4" "0.6" "0.8" "1.0"];
% vitesse=["-1.0" "-0.8" "-0.6" "-0.4" "-0.2" "0.0" "0.2" "0.4" "0.6" "0.8" "1.0"]
%vitesse=["-2.0" "-1.6" "-1.4" "-1.2" "-1.0" "-0.8" "-0.6" "-0.4" "-0.2"]% "0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0"];
%vitesse=["0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0" "8.0" "50.0" "100.0" "300.0" "600.0" "900.0" "1000.0"];
vitesse=["0.0" "0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0" "3.0" "4.0" "5.0" "6.0" "7.0" "8.0" "9.0" "10.0" "11.0" "12.0" "15.0" "20.0" "50.0"];
%speed=[0.0 -0.2 -0.4 -0.6 -0.8 -1.0 -1.2 -1.4 -1.6 -2.0];
% speed=[0.0 0.2 0.4 0.6 0.8 1.0];
% speed=[-1.0 -0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 1.0]
%speed=[-2.0 -1.6 -1.4 -1.2 -1.0 -0.8 -0.6 -0.4 -0.2]% 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];
%speed=[0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 8.0 50.0 100.0 300.0 600 900 1000];
speed=[0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 15.0 20.0 50.0];
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
    if speed(j)>0 && speed(j)<2.5
        eval(sprintf('cd /home/gascuel/PHD/Lipro/Lipro3/Données+/V=%s', v));
        rep=1000;
        fin="";
    elseif speed(j)<=0 && speed(j)<2.5
        eval(sprintf('cd /home/gascuel/PHD/Lipro/Lipro3/Données-/V=%s', v));
        rep=1000;
        fin="";
    else
        cd /tmp/wkm5/DoubleS
        rep=20;
        fin="_vc";
    end
    %FrontVague=dlmread(sprintf("data_FrontVague_lipro_N_%g_V_%s00000_rep_%g%s.dat",N,v,rep,fin));
    nM=dlmread(sprintf("data_nM_lipro_N_%g_V_%s00000_rep_%g%s.dat",N,v,rep,fin)); 
%     
% %     Reg=[];
    [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
% %     for i=1:size(row,1)-1
% %         x=nM(row(i):row(i+1)-1,1);
% %         y=nM(row(i):row(i+1)-1,2);
% %         if max(y)>seuil
% %             reg=x\y;
% %             Reg=[Reg reg];
% %         end
% %     end
    %analyse des propagations
    [row2,col2]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
    clear actifin
    for i=2:size(row2,1)
        actifin(i-1)=nM(row2(i)-1,2)*100/N;
    end
    m2=mean(actifin);
    sd2=std(actifin);
    %h1=figure; hold on;
[row,col]=find(nM(:,1)==0);
for i=2:size(row,1)
    diff(i-1)=row(i)-row(i-1);
end
D=min(diff); %calcul du nomber de pas de temps min d'une simu
% 
% for i=1:100:D
%     %max(nM(find(nM(:,1)==0)+i)) %affichage du pas de temps
%     data=nM(row+i,2);
%     figure;
%     h=histogram(data,10)
%     hval=h.Values
%     close
%     figure(h1)
%     plot(hval)
%     pause(0.001)
% end
%     es2=sd2/sqrt(size(actifin,2));
%     clear row3
%     [row3,col3]=find(actifin'==100);
%     if isempty(row3)==0
%         s(j)=size(row3,1)/rep;
%     else
%         s(j)=0;
%     end
% %     [h,p]=kstest((Reg-mean(Reg))/std(Reg))
    graphique2(1,j)=m2;
    graphique2(2,j)=sd2;
% % % 
% % %     m1=mean(Reg);
% % %     sd1=std(Reg);
% % %     es1=sd1/sqrt(size(Reg,2));
% % %     %[h,p]=kstest((Reg-mean(Reg))/std(Reg))
% % %     graphique(1,j)=m1;
% % %     graphique(2,j)=es1;
% %     j=j+1;
% 
    %%% distributions des populations finales
    
    k=1;
    vec=[(row(2:size(row,1),1)-1)' size(nM,1)];
    for i=vec
        bp(k,j)=nM(i,2);
        k=k+1;
    end
    j=j+1;
end
% % 
% % %%% suite distrib finales
% %boxplot(bp,vitesse);
% X=[];
% Y=[];
% vect=2:6;
% for i=vect
%     X=[X speed(i)*ones(1,1000)];
%     Y=[Y bp(:,i)'];
% end    
% hist3([X',Y'],'CdataMode','auto','FaceColor','interp','Ctrs',{speed(vect) 50:100:1950});
% colorbar
% % view(2)%vue 2D
% az = 120;
% el = 20;
% view(az, el);%ou3D
%ylabel({"Number of active" "at the end of a simulation"});
%xlabel({"Speed of movement"});
%zlabel('Frequency');
%title({"Distributions of the number of active at the end of a simulation" "against the speed of movement"});
%savename = ('fig_distrib_nAfin_compil_hist3D_zoom');
% orient(gcf,'landscape')
% saveas(gcf, savename, 'pdf');print(savename,'-dpdf','-bestfit');
% 
% figure;hold on;
% X=[];
% for i=1:size(speed,2)
%     %h=histogram(bp(:,i),'BinEdges',(0:100:2000),'Visible','off');
%     h=histogram(bp(:,i),'BinEdges',(0:100:2000),'facealpha',0.3);
%     j=1;
%     for k=h.Values
%         hval(i,j)=k;
%         j=j+1;
%     end
%     X=[X (50:100:1950)'];
% end
% legend(vitesse)
% 
% figure;
% surf((50:100:1950)',speed',hval);
% xlabel('Number of departed');
% ylabel('speed');
% zlabel('Frequency');
% az = 75;
% el = 45;
% view(az, el);
% 
% % figure;
% % plot(X,hval');
% % legend(vitesse)
% % savename = ('fig_distrib_nAfin_compil_curv2D'); 
% % dlmwrite(savename,[speed;hval'],'delimiter',' ');
% % saveas(gcf, savename, 'png'); 

%%%suite number missed

%subplot(2,1,1);
figure
for i=1:size(vitesse,2)
    if graphique2(1,i)+graphique2(2,i)>100
        pos(i)=100-graphique2(1,i);
    else
        pos(i)=graphique2(2,i);
    end
    if graphique2(1,i)-graphique2(2,i)<0
        neg(i)=graphique2(1,i);
    else
        neg(i)=graphique2(2,i);
    end   
end
% errorbar(speed/7.39,graphique2(1,:),neg(:),pos(:),'-x')
% axis([-1/7.39 51/7.39 97 100.1]);
errorbar(speed,graphique2(1,:),neg(:),pos(:),'-x')
axis([-1 51 97 100.1]);

ylabel({"Nombre d'individus" "propagés en moyenne"});
title("Compilation des moyennes d'individus propagés");
savename="compil_finalid_v_0_50_zoom_urescaled"
% dlmwrite(savename,[speed graphique2(1,:) neg(:) pos(:)],'delimiter',' ');
% saveas(gcf, savename, 'pdf');

% 
% 
% subplot(2,1,2);
% errorbar(speed,graphique(1,:),graphique(2,:))
% xlabel('Vitesses de mouvement');
% ylabel({"Vitesse moyenne" "de la vague d'activation"});
% title('Compilation des moyennes de vitesses mesurées par ID');
% % saveName = ('fig_cfonctionvtot_triangle'); 
% % saveas(gcf, saveName, 'jpg'); 
% % csvwrite('vitesse_moy_tri.dat',[speed' graphique']);

% % %%%%Regression linéaire
% % hold on;
% % b=[ones(length(speed(6:20)'),1) speed(6:20)']\graphique(1,6:20)';
% % ycalc=[ones(length(speed'),1) speed']*b;
% % plot(speed,ycalc)
% % saveName = ('fig_cfonctionvtot_résumé+fit'); 
% % legend('Data','Reg')
% % saveas(gcf, saveName, 'jpg');

%% compilateur des fluctuations des missed
% on veut voir les distributions finales de missed

% cd /tmp/wkm5/DoubleS/data_nImiss
% 
% vect=[0 3 4 5 6 7 8 9 10 11 12 15 20 50];
% i=1;
% for v=vect
%     subplot(4,4,i)
%     nIm=dlmread(sprintf("fluct_20sim_%g.dat",v)); 
%     hist(nIm)
%     title(sprintf("v=%g",v))
%     axis([0 60 0 20])
%     i=i+1;
% end
% 
% saveName="compil_stat_nImiss_20sim"
% saveas(gcf,saveName,'pdf')

%% Compilateur liprostat 2D
% cd /tmp/wkm5/Lipro2D
% close all
% clear all
% N=10000;
% rep=10;
% seuil=1000;
% graphique=[[] []];
% j=1;
% 
% %vitesse =[-2:0.5:-0.5 0.2:0.2:2]
% % h1=figure
% % hold on
% % h2=figure
% % hold on
% figure;
% hold on;
% %vitesse=[-0.3 -0.1 0 0.1 0.3]
% vitesse=[-1:0.1:-0.3];
% vcol=["y" "c" "m" "r" "b" "k" "g" "y" "c"]
% %vitesse =[-100 -85 -75 -50 -45 -40 -35 -25 -10];
% 
% sup=10000;
% inf=0;
% for v=vitesse
%     v
% %     if v==3
% %         sup=6000;
% %     else
% %         sup=9000;
% %     end
%     pause(0.0001)
%     %FrontVague=dlmread(sprintf("data_FrontVague_lipro_N_%g_V_%.2f0000_rep_%g_vchc_F_alpha.dat",N,v,rep));
%     nM=dlmread(sprintf("data_nM_lipro_N_%g_V_%.2f0000_rep_%g_vchc_F_alpha.dat",N,v,rep));
%     
%     %nM
% %     figure;
% %     hold on;
%     Reg=[];
%     [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
%     size(row,1)
%     for i=1:size(row,1)-1
%         [row2,col2]=find(nM(row(i):row(i+1)-1,2)>=inf & nM(row(i):row(i+1)-1,2)<=sup);
%         if isempty(row2)
%             %Reg=[Reg -100];
%         else
%             x=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,1);
%             y=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,2);
%             if max(y)>seuil
%                 %reg=x\y;
%                 %Reg=[Reg reg];
%                 reg=[ones(length(x),1) x]\y;
%                 Reg=[Reg reg(2)/100];
%                 %yCalc=reg*x;
%                 %plot(x,yCalc);
%                 %figure(h1)
%                 plot(x,y,vcol(j));%plot des raw data
%             end
%         end
%     end
%     %%% on ajoute à la main la dernière simu
%     [row2,col2]=find(nM(row(i+1):size(nM,1),2)>=inf & nM(row(i+1):size(nM,1),2)<=sup);
%     if isempty(row2)
%         %Reg=[Reg -100];
%     else
%         x=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,1);
%         y=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,2);
%         if max(y)>seuil
%             %         reg=x\y; %lineaire
%             %         Reg=[Reg reg];
%             reg=[ones(length(x),1) x]\y; %correction affine
%             Reg=[Reg reg(2)/100];
%             %yCalc=reg*x;
%             %plot(x,yCalc);
%             %figure(h1)
%             plot(x,y,vcol(j));%plot des raw data
%         end
%     end
%     %     %Front
%     %     Reg=[];
%     %     %figure;
%     %     %hold on;
%     %     [row,col]=find(FrontVague(:,1)==0); %on récupère toute les lignes de début de simu
%     %     for i=1:size(row,1)-1
%     %         x=FrontVague(row(i):row(i+1)-1,1);
%     %         y=FrontVague(row(i):row(i+1)-1,2);
%     %         if max(y)>seuil
%     %             reg=x\y;
%     %             Reg=[Reg reg];
%     %             yCalc=reg*x;
%     %             %plot(x,yCalc);
%     %             %figure(h2)
%     %             %plot(x,y);%plot des raw data
%     %         end
%     %     end
%     % % on ajoute à la main la dernière simu
%     %     x=FrontVague(row(i+1):size(FrontVague,1),1);
%     %     y=FrontVague(row(i+1):size(FrontVague,1),2);
%     %     if max(y)>seuil
%     %         reg=x\y;
%     %         Reg=[Reg reg];
%     %         yCalc=reg*x;
%     %         %plot(x,yCalc);
%     %         %figure(h2)
%     %         %plot(x,y);%plot des raw data
%     %     end
%     
%     title(sprintf('v = %.1f',v));
%     
%     m1=mean(Reg);
%     sd1=std(Reg);
%     size(Reg,2)
%     es1=sd1/sqrt(size(Reg,2));
%     %[h,p]=kstest((Reg-mean(Reg))/std(Reg))
%     graphique(1,j)=m1;
%     graphique(2,j)=es1;
%     j=j+1;
% end
% 
% 
% figure;
% errorbar(vitesse,graphique(1,:),graphique(2,:),'x')
% xlabel('Vitesses de mouvement');
% ylabel({"Vitesse moyenne" "de la vague d'activation"});
% title('Compilation des moyennes de vitesses mesurées par ID');
% % saveName = ('fig_cfonctionvtot_2D_Front');
% % saveas(gcf, saveName, 'pdf');
% % dlmwrite('vitesse_moy_2D_Front.dat',[vitesse' graphique']);
% 
% %%%recherche point critique
% cd /tmp/wkm5/Lipro2D
% close all
% N=10000;
% rep=100;
% %vitesse =[-100 -85 -75 -50 -45 -40 -35 -25 -10 -5 -3 -1 0 2];
% %vitesse =[- 120 -110 -100 -90 -80 -70 -60 -55 -50 -45 -40 -35];
% vitesse=[-1.3:0.1:-1 -0.95 -0.9 -0.85 -0.8:0.1:-0.3 -0.1];
% % vitesse=[-2 -1.5 -1 -0.5]
% M=[];
% S=[];
% taille=10000; %correspond au nombre d'événements
% for v=vitesse
%     v
%     pause(0.001)
%     Mass=dlmread(sprintf("data_Mass_lipro_N_%g_V_%.2f0000_rep_%g_vchc_F_alpha.dat",N,v,rep));
%     size(Mass,1)
%     figure;
%     histogram(Mass,'BinEdges',[0:50:taille])
%     M=[M mean(Mass)];
%     S=[S std(Mass)];
%     title(sprintf('v = %.1f',v));
%     saveName = sprintf('hist_V_%.1f',v);
%     %     saveas(gcf, saveName, 'pdf');
% end
% for i=1:size(vitesse,2)
%     if M(i)+S(i)>taille
%         pos(i)=taille-M(1,i);
%     else
%         pos(i)=S(i);
%     end
%     if M(i)-S(i)<0
%         neg(i)=M(i);
%     else
%         neg(i)=S(i);
%     end
% end
% errorbar(vitesse,M,neg(:),pos(:),'-x');
% axis([vitesse(1)-0.1 vitesse(size(vitesse',1))+0.1 -10 taille+10])
% %print("Mass_criticality_paperversion", '-dpdf','-painters','-bestfit')
% 
% %%% Image finale
% cd /tmp/wkm5/Lipro2D
% %close all
% N=10000;
% rep=10;
% v=-1.2;
% simu=1;
% i=simu-1;
% Movers=dlmread(sprintf("data_Movers_lipro_N_%g_V_%.1f00000_rep_%g_vchc_F_alpha.dat",N,v,rep));
% 
% for i=0:2%rep-1
%     %B=abs(Movers(i*N+1,1)-Movers(i*N+10000,1))
%     %% plot de l'image finale
%     for j=1:N %on actualise les positions
%         if Movers(i*N+j,3)==0
%             VectColor(j,1)= 0;VectColor(j,2)= 1;VectColor(j,3)= 0;
%         else
%             VectColor(j,1)= 1;VectColor(j,2)= 0;VectColor(j,3)= 0;
%         end
%     end
%     figure
%     scatter(Movers(i*N+1:(i+1)*N,1),Movers(i*N+1:(i+1)*N,2),10,VectColor,'Filled');
%     
%     %% scan de la densité en largeur
%     figure;
%     H=histogram(Movers(i*N+1:(i+1)*N,1),"BinEdges",Movers(i*N+1,1)+0.5:1:Movers(i*N+10000,1));%scanne le long de l'axe x
%     hvalue=H.Values;
%     close
%     figure
%     X=Movers(i*N+1,1)+0.5+0.5:1:Movers(i*N+10000,1)-1+0.5;
%     plot(X/abs(v),hvalue*abs(v))
%     for j=1:size(hvalue')-1
%         Diff(j)=hvalue(j+1)-hvalue(j);
%     end
%     %plot(Diff);
%     %axis([-10 110 -10 110]);% mettre pos(1) si on suit tous les movers
%     %% analyse des fluctuations
% %     bin=10;
% %     H=hist3([Movers(i*N+1:(i+1)*N,1),Movers(i*N+1:(i+1)*N,2)],'CdataMode','auto','Ctrs',{Movers(i*N+1,1)+0.5:bin:Movers(i*N+10000,1) (0.5:bin:99.5)},'LineStyle','none');
% %     h1=surf(Movers(i*N+1,1)+0.5:bin:Movers(i*N+10000,1),(0.5:bin:99.5),H')
% %     set(h1,'LineStyle','none');
% %     view(2); colorbar; 
% %     Z=h1.ZData;
% end

%% calcul de la sensibilité du débit
% cd /tmp/wkm5/Lipro2D
% N=10000;
% k=1;
% rep=100;
% close all
% clear MSensi
% clear MMDebit
% filt1=10; %seuil de filtrage de la partie stationnaire
% filt2=50;
% vitesse=[-1.2:0.1:-1 -0.95 -0.9 -0.85 -0.8:0.1:-0.3 -0.1 -0.03 -0.01];
% for v=vitesse
%     if v>-0.15
%         rep=10;
%     end
%     Movers=dlmread(sprintf("data_Movers_lipro_N_%g_V_%.2f0000_rep_%g_vchc_F_alpha.dat",N,v,rep));
%     clear Mdebit
%     clear MDiff
%     for i=0:size(Movers,1)/10000-1
%         [row,col]=find(Movers(i*N+1:(i+1)*N,3)==0); %on filtre les actifs
%         figure;
%         H=histogram(Movers(i*N+row,1),"BinEdges",Movers(i*N+1,1)+filt1:1:Movers(i*N+10000,1)-filt2);%scanne le long de l'axe x, on filtre le début et la fin 
%         hvalue=H.Values;
%         close
%         %hvalue=hvalue(find(hvalue>0));
%         %figure
%         %X=Movers(i*N+1,1)+filt+0.5:1:Movers(i*N+10000,1)-filt-1+0.5;
%         %plot(X/abs(v),hvalue*abs(v))
%         clear Diff
%         for j=1:size(hvalue')-1
%             Diff(j)=(hvalue(j)*abs(v)-mean(hvalue*abs(v)))^2;
%         end
%         Mdebit(i+1)=mean(hvalue*abs(v));
%         MDiff(i+1)=mean(Diff);
%     end
%     %histogram(MDiff)
%     MSensi(k,1)=mean(MDiff);
%     MSensi(k,2)=std(MDiff);
%     MMDebit(k,1)=mean(Mdebit);
%     MMDebit(k,2)=std(Mdebit);
%     k=k+1;
% end
% errorbar(vitesse,MSensi(:,1),MSensi(:,2),'-x');
% figure
% errorbar(vitesse,MMDebit(:,1),MMDebit(:,2),'-x');

%% moyenneur de signal
% cd /tmp/wkm5/Lipro2D
% close all
% clear all
% N=10000;
% rep=10;
% seuil=5000;
% graphique=[[] []];
% j=1;
% 
% %vitesse =[-2:0.5:-0.5 0.2:0.2:2]
% 
% h1=figure;
% hold on;
% h2=figure;
% hold on;
% h3=figure;
% hold on;
% %vitesse=[-0.3 -0.1 0 0.1 0.3]
% vitesse=[-3 -1.3:0.1:-1 -0.95 -0.9 -0.85 -0.8:0.1:-0.3 -0.1 0 0.1 0.3 1 3];
% %vitesse=[-1.1 -0.9 -0.7 -0.6]
% %vitesse=[-0.95 -0.9 -0.85]
% vcol=["k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k"]
% %vitesse =[-120 -110 -100 -90 -80 -70 -60 -50];
% v=-1;
% sup=10000; %borne de controle des regressions
% inf=100;
% T=[0.1095    0.6164    1.7979    7.2458    8.6394    9.1170    8.8325    9.7980    9.8911]*1000;
% sample=100;
% Eqtemps=100000; % permet le prolongement des courbes après saturation valeur pour les alpha 21000
% 
% channel=4;
% 
% for v=vitesse
%     v
%     if(v>-3)
%         rep=100 ;
%     end
%     if(v>-0.2)
%         rep=10;
%     end
%     pause(0.0001)
%     nM=dlmread(sprintf("data_nM_lipro_N_%g_V_%.2f0000_rep_%g_vchc_F_alpha.dat",N,v,rep));
%     pdt=nM(2,1)-nM(1,1);
%     Reg=[];
%     [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
%     size(row,1) % nombre de simus
%     
%     Diff=0;
%     for i=1:size(row,1)-1
%         diff=row(i+1)-row(i);
%         Diff=max(Diff,diff);
%         if(Diff==diff)
%             lab=i;
%         end
%     end
%     diff=size(nM,1)-row(i+1);
%     Diff=max(Diff,diff);
%     if(Diff==diff)
%             lab=i+1;
%     end
%     MatTot=zeros(Eqtemps,size(row,1));
%     MatDeriv=zeros(Eqtemps,size(row,1));
%     MaxStock=zeros(1,size(row,1));
%     for i=1:size(row,1)-1
%         [row2,col2]=find(nM(row(i):row(i+1)-1,channel)>=inf & nM(row(i):row(i+1)-1,channel)<=sup);
%         if isempty(row2)
%             %Reg=[Reg -100];
%         else
%             MatTot(row2(1):row2(size(row2,1))-1,i)=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
%             MatTot(row2(size(row2,1)):size(MatTot,1),i)=nM(row(i)+row2(size(row2,1))-1,channel);
%             x=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,1);
%             y=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
%             if max(y)>seuil
% X=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,1);
%             Y=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
%                 %reg=x\y;
%                 %Reg=[Reg reg];
%                 reg=[ones(length(x),1) x]\y;
%                 Reg=[Reg reg(2)/100];
%                 %yCalc=reg*x;
%                 %plot(x,yCalc);
%                 clear Deriv;
%                 figure(h2)
% 
%                 for k=1:size(y)-sample
%                     Deriv(k)=(y(k+sample)-y(k))/(sample*pdt);
%                 end
%                 MatDeriv(1:size(Deriv',1),i)=Deriv;
%                 MatDeriv(size(Deriv',1)+1:size(MatDeriv,1),i)=Deriv(k);
%                 %plot(x(1:size(x)-sample),Deriv,vcol(j));
%                 %plot(x,y,vcol(j));%plot des raw data
%                 %pause(1)
%             end
%         end
%     end
% 
%     %%% on ajoute à la main la dernière simu
%     [row2,col2]=find(nM(row(i+1):size(nM,1),channel)>=inf & nM(row(i+1):size(nM,1),channel)<=sup);
%            
%     if isempty(row2)
%         %Reg=[Reg -100];
%     else
%         MatTot(row2(1):row2(size(row2,1))-1,i+1)=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,channel);
%         MatTot(row2(size(row2,1)):size(MatTot,1),i+1)=nM(row(i+1)+row2(size(row2,1))-1,channel);
%         x=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,1);
%         y=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,channel);
%         if max(y)>seuil
% 
%             
%             %         reg=x\y; %lineaire
%             %         Reg=[Reg reg];
%             reg=[ones(length(x),1) x]\y; %correction affine
%             Reg=[Reg reg(2)/100];
%             %yCalc=reg*x;
%             %plot(x,yCalc);
%             %figure(h1)
%             clear Deriv;   
%             for k=1:size(y)-sample
%                 Deriv(k)=(y(k+sample)-y(k))/(sample*pdt);
%             end
%             MatDeriv(1:size(Deriv',1),i+1)=Deriv;
%             MatDeriv(size(Deriv',1)+1:size(MatDeriv,1),i+1)=Deriv(k);
%             %plot(x(1:size(x)-sample),Deriv,vcol(j)); 
%             %plot(x,y,vcol(j));%plot des raw data
%             %pause(1)
%         end
%     end
%     MaxStock=max(MatDeriv);
%     
%     title(sprintf('v = %.1f',v));
%     
%     m1=mean(Reg);
%     sd1=std(Reg);
%     size(Reg,2)
%     es1=sd1/sqrt(size(Reg,2));
%     %[h,p]=kstest((Reg-mean(Reg))/std(Reg))
%     graphique(1,j)=m1;
%     graphique(2,j)=es1;
%     
%     
% %     figure(h1);
% %     [row3,col3]=find(MatTot(size(MatTot,1),:)>T(j)/3);
% %     clear neg;
% %     clear pos;
% %     M=mean(MatTot(:,col3)');
% %     S=std(MatTot(:,col3)');
% %     taille=10000;
% %     for i=1:size(MatTot,1)
% %         if M(i)+S(i)>taille
% %             pos(i)=taille-M(1,i);
% %         else
% %             pos(i)=S(i);
% %         end
% %         if M(i)-S(i)<0
% %             neg(i)=M(i);
% %         else
% %             neg(i)=S(i);
% %         end
% %     end
% %     
% %     pdt=nM(2,1)-nM(1,1); % on va prolonger les courbes après leur saturation
% %     errorbar(pdt:pdt:Eqtemps*pdt,M,neg,pos,'-x');% avec prlgmt
%     %plot(nM(row(lab):row(lab)+Diff-1,1),M,vcol(j));
%     %errorbar(nM(row(lab):row(lab)+Diff-1,1),M,neg,pos,'-x');%sans prlgmt
%     
%     figure(h1)
%     [row3,col3]=find(MatTot(size(MatTot,1),:)>0);%T(j)/3);
%     M=mean(MatTot(:,col3)');Stock1M(1:size(M'),j)=M;
%       plot(pdt:pdt:Eqtemps*pdt,(N-M)./(pdt:pdt:Eqtemps*pdt).^0.829,vcol(j));%rescaled
%       %plot(pdt:pdt:Eqtemps*pdt,(N-M),vcol(j));%pas resscaled
% %     plot(pdt:pdt:Eqtemps*pdt,M+pos,vcol(j));
% %     plot(pdt:pdt:Eqtemps*pdt,M-neg,vcol(j));
% %     axis([1 4000*pdt 100 9000])
% %     set(gca, 'YScale', 'log')
% %     set(gca, 'XScale', 'log')
% %     plot(nM(row(lab):row(lab)+Diff-1,1),M,vcol(j)); % on se limite à
% %     l'arrêt de la simulation
% %     plot(nM(row(lab):row(lab)+Diff-1,1),M+pos,vcol(j));
% %     plot(nM(row(lab):row(lab)+Diff-1,1),M-neg,vcol(j));
%     clear M
%     figure(h3)
%     M=mean(MatDeriv');
%     Stock2M(1:size(M'),j)=M;
%     MaxiM(j,1)=mean(MaxStock);
%     MaxiM(j,2)=std(MaxStock);
%     
%     plot(pdt:pdt:Eqtemps*pdt,M,vcol(j));
% 
%     j=j+1;
% 
% end
% figure(h1)
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% axis([40 160 2 140]);
% leg=legend(num2str(vitesse'),'Location','southwest');
% title(leg,'v');
% %print("comp_v_nI_mean_paperversion", '-dpdf','-painters','-bestfit')
% Stock1M=[pdt:pdt:Eqtemps*pdt;Stock1M'];
% Stock1M=Stock1M';
% save('IBMStock1M.dat','Stock1M','-ASCII') %stockage signaux moyens
% 
% figure
% MaxStockM= max(Stock2M);
% plot(vitesse,MaxStockM/100,"-x")
% hold on 
% errorbar(vitesse,MaxiM(:,1)/100, MaxiM(:,2)/100,"-x")
% axis([-3.25 3.25 -10 450]);
% dlmwrite('vitesse_moy_2D_Front_alpha_CMAX.dat',[vitesse;MaxStockM/100]);
% 
% cd /tmp/wkm5/Lipro2D
% A=dlmread("vitesse_moy_2D_Front_alpha.dat");
% A=A(3:13,1:2);
% A(:,2)=A(:,2)*100;
% 
% %B=[vitesse;MaxStockM/100]
% B=dlmread("vitesse_moy_2D_Front_alpha_CMAX.dat");
% C=[B';A]
% C(20,1)=-0.31;
% scatter(C(:,1),C(:,2),"x")
% 
% C=C([2:20],:);
% %C=C([1:11 13:19],:);
% C(:,1)=C(:,1)+1.00001;
% 
% C=log(C);
% coef=C(:,1)\C(:,2)

%% PROJET FDA
%moyenneur de signal
cd /tmp/wkm5/FDA-STOPS/FDA/
close all
clear all
N=2000;
rep=20;
seuil=0;
graphique=[[] []];
j=1;

h1=figure;
hold on;
h2=figure;
hold on;
h3=figure;
hold on;

vitesse=[-10 -4:1:4 10];

vcol=["k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k"];
v=1;
sup=2000; %borne de controle des regressions
inf=0;
T=[0.1095    0.6164    1.7979    7.2458    8.6394    9.1170    8.8325    9.7980    9.8911]*1000;
sample=100;
Eqtemps=721000; % permet le prolongement des courbes après saturation valeur pour les alpha 21000

channel=2;

for v=vitesse
    v

    pause(0.0001)
    if abs(v)>4 || v>0
        nM=dlmread(sprintf("data_nM_FDA_N_%g_V_%.2f0000_rep_%g_vchc.dat",N,v,rep));
    else
        nM=dlmread(sprintf("data_nM_FDA_N_%g_V_%.2f0000_rep_%g.dat",N,v,rep));
    end
    pdt=nM(2,1)-nM(1,1);
    Reg=[];
    [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
    size(row,1) % nombre de simus
    
    Diff=0;
    for i=1:size(row,1)-1
        diff=row(i+1)-row(i);
        Diff=max(Diff,diff);
        if(Diff==diff)
            lab=i;
        end
    end
    diff=size(nM,1)-row(i+1);
    Diff=max(Diff,diff);
    if(Diff==diff)
            lab=i+1;
    end
    MatTot=zeros(Eqtemps,size(row,1));
    MatDeriv=zeros(Eqtemps,size(row,1));
    MaxStock=zeros(1,size(row,1));
    for i=1:size(row,1)-1
        [row2,col2]=find(nM(row(i):row(i+1)-1,channel)>=inf & nM(row(i):row(i+1)-1,channel)<=sup);
        if isempty(row2)
            %Reg=[Reg -100];
        else
            MatTot(row2(1):row2(size(row2,1))-1,i)=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
            MatTot(row2(size(row2,1)):size(MatTot,1),i)=nM(row(i)+row2(size(row2,1))-1,channel);
            x=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,1);
            y=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
            if max(y)>seuil
                X=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,1);
                Y=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
                %reg=x\y;
                %Reg=[Reg reg];
                reg=[ones(length(x),1) x]\y;
                Reg=[Reg reg(2)];
                %yCalc=reg*x;
                %plot(x,yCalc);
                clear Deriv;
                figure(h2)

                for k=1:size(y)-sample
                    Deriv(k)=(y(k+sample)-y(k))/(sample*pdt);
                end
                MatDeriv(1:size(Deriv',1),i)=Deriv;
                MatDeriv(size(Deriv',1)+1:size(MatDeriv,1),i)=Deriv(k);
                %plot(x(1:size(x)-sample),Deriv,vcol(j));
                plot(x,y,vcol(j));%plot des raw data
                %pause(1)
            end
        end
    end

    %% on ajoute à la main la dernière simu
%     [row2,col2]=find(nM(row(i+1):size(nM,1),channel)>=inf & nM(row(i+1):size(nM,1),channel)<=sup);
%            
%     if isempty(row2)
%         %Reg=[Reg -100];
%     else
%         MatTot(row2(1):row2(size(row2,1))-1,i+1)=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,channel);
%         MatTot(row2(size(row2,1)):size(MatTot,1),i+1)=nM(row(i+1)+row2(size(row2,1))-1,channel);
%         x=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,1);
%         y=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,channel);
%         if max(y)>seuil
% 
%             
%             %reg=x\y; %lineaire
%             %Reg=[Reg reg];
%             reg=[ones(length(x),1) x]\y; %correction affine
%             Reg=[Reg reg(2)];
%             %yCalc=reg*x;
%             %plot(x,yCalc);
%             figure(h2)
%             clear Deriv;   
%             for k=1:size(y)-sample
%                 Deriv(k)=(y(k+sample)-y(k))/(sample*pdt);
%             end
%             MatDeriv(1:size(Deriv',1),i+1)=Deriv;
%             MatDeriv(size(Deriv',1)+1:size(MatDeriv,1),i+1)=Deriv(k);
%             %plot(x(1:size(x)-sample),Deriv,vcol(j)); 
%             plot(x,y,vcol(j));%plot des raw data
%             %pause(1)
%         end
%     end
%     MaxStock=max(MatDeriv);
    
    title(sprintf('v = %.1f',v));
    
    m1=mean(Reg);
    sd1=std(Reg);
    size(Reg,2)
    es1=sd1/sqrt(size(Reg,2));
    graphique(1,j)=m1;
    graphique(2,j)=es1;
    
    figure(h1)
     [row3,col3]=find(MatTot(size(MatTot,1),:)>0);%T(j)/3);
     M=mean(MatTot(:,col3)');Stock1M(1:size(M'),j)=M;
     plot(pdt:pdt:Eqtemps*pdt,M,vcol(j))
%     plot(pdt:pdt:Eqtemps*pdt,(N-M)./(pdt:pdt:Eqtemps*pdt).^0.829,vcol(j));%rescaled

    clear M
    figure(h3)
    M=mean(MatDeriv');
    Stock2M(1:size(M'),j)=M;
    MaxiM(j,1)=mean(MaxStock);
    MaxiM(j,2)=std(MaxStock);
    
    plot(pdt:pdt:Eqtemps*pdt,M,vcol(j));

    j=j+1;
end
figure(h1)
xlabel('Time')
ylabel('# Actives')
title('mean signal')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
%axis([40 160 2 140]);
leg=legend(num2str(vitesse'),'Location','southeast');
title(leg,'v');
%print("Actives_meansignal_FDA", '-dpdf','-painters','-bestfit')

figure(h2)
xlabel('Time')
ylabel('# Actives')
title('All signals')
%print("Actives_allsignal_FDA", '-dpdf','-painters','-bestfit')

figure(h3)
xlabel('Time')
ylabel('# Actives')
title('Derivative of nM')
%print("Actives_derivative_FDA", '-dpdf','-painters','-bestfit')

figure
for i=1:size(vitesse,2)
    pos(i)=graphique(2,i);
    if graphique(1,i)-graphique(2,i)<0
        neg(i)=graphique(1,i);
    else
        neg(i)=graphique(2,i);
    end   
end
errorbar(vitesse([1:9 11]),graphique(1,[1:9 11]),neg([1:9 11]),pos([1:9 11]),'-x')%[1:9 11] à remplacer par ":"
axis([-11 11 0 1])
xlabel('Speed')
ylabel('Activation speed')
title('Mean slope of signals')
Stock=[vitesse' graphique']
%Stock=Stock([1:9 11],:)
%save('cfunctionofv_FDA.dat','Stock1M','-ASCII') %stockage 
%print("cfunctionofv_FDA", '-dpdf','-painters','-bestfit')

%%% Image final
cd /tmp/wkm5/FDA-STOPS/FDA/
%close all
N=2000;
rep=1;
v=-7;
simu=1;
i=simu-1;
Movers=dlmread(sprintf("data_Movers_FDA_N_%g_V_%.2f0000_rep_%g.dat",N,v,rep));
N/(max(Movers(:,1))-min(Movers(:,1)))
for i=0:0%rep-1
    %B=abs(Movers(i*N+1,1)-Movers(i*N+10000,1))
    %% plot de l'image finale
    for j=1:N %on actualise les positions
        if Movers(i*N+j,2)==0
            VectColor(j,1)= 0;VectColor(j,2)= 1;VectColor(j,3)= 0;
        else
            VectColor(j,1)= 1;VectColor(j,2)= 0;VectColor(j,3)= 0;
        end
    end
    figure
    scatter(Movers(i*N+1:(i+1)*N,1),zeros(2000,1),10,VectColor,'Filled');
    
%     %% scan de la densité en largeur
%     figure;
%     H=histogram(Movers(i*N+1:(i+1)*N,1),"BinEdges",Movers(i*N+1,1)+0.5:1:Movers(i*N+10000,1));%scanne le long de l'axe x
%     hvalue=H.Values;
%     close
%     figure
%     X=Movers(i*N+1,1)+0.5+0.5:1:Movers(i*N+10000,1)-1+0.5;
%     plot(X/abs(v),hvalue*abs(v))
%     for j=1:size(hvalue')-1
%         Diff(j)=hvalue(j+1)-hvalue(j);
%     end
%     %plot(Diff);
%     %axis([-10 110 -10 110]);% mettre pos(1) si on suit tous les movers
%     %% analyse des fluctuations
% %     bin=10;
% %     H=hist3([Movers(i*N+1:(i+1)*N,1),Movers(i*N+1:(i+1)*N,2)],'CdataMode','auto','Ctrs',{Movers(i*N+1,1)+0.5:bin:Movers(i*N+10000,1) (0.5:bin:99.5)},'LineStyle','none');
% %     h1=surf(Movers(i*N+1,1)+0.5:bin:Movers(i*N+10000,1),(0.5:bin:99.5),H')
% %     set(h1,'LineStyle','none');
% %     view(2); colorbar; 
% %     Z=h1.ZData;
end


%% Projet STOP
%moyenneur de signal
cd /tmp/wkm5/Stops
close all
clear all
N=2000;
rep=20;
seuil=0;
graphique=[[] []];
j=1;

h1=figure;
hold on;
h2=figure;
hold on;
h3=figure;
hold on;

vitesse=[-50 -10 -5:1:-1];

vcol=["k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k" "r" "c" "m" "b" "g" "k"];
v=-1;
sup=2000; %borne de controle des regressions
inf=0;
T=[0.1095    0.6164    1.7979    7.2458    8.6394    9.1170    8.8325    9.7980    9.8911]*1000;
sample=100;
Eqtemps=721000; % permet le prolongement des courbes après saturation valeur pour les alpha 21000

channel=2;

for v=vitesse
    v

    pause(0.0001)
    nM=dlmread(sprintf("data_nM_STOP_N_%g_V_%.2f0000_rep_%g_vchc.dat",N,v,rep));
    pdt=nM(2,1)-nM(1,1);
    Reg=[];
    [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
    size(row,1) % nombre de simus
    
    Diff=0;
    for i=1:size(row,1)-1
        diff=row(i+1)-row(i);
        Diff=max(Diff,diff);
        if(Diff==diff)
            lab=i;
        end
    end
    diff=size(nM,1)-row(i+1);
    Diff=max(Diff,diff);
    if(Diff==diff)
            lab=i+1;
    end
    MatTot=zeros(Eqtemps,size(row,1));
    MatDeriv=zeros(Eqtemps,size(row,1));
    MaxStock=zeros(1,size(row,1));
    for i=1:size(row,1)-1
        [row2,col2]=find(nM(row(i):row(i+1)-1,channel)>=inf & nM(row(i):row(i+1)-1,channel)<=sup);
        if isempty(row2)
            %Reg=[Reg -100];
        else
            MatTot(row2(1):row2(size(row2,1))-1,i)=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
            MatTot(row2(size(row2,1)):size(MatTot,1),i)=nM(row(i)+row2(size(row2,1))-1,channel);
            x=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,1);
            y=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
            if max(y)>seuil
                X=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,1);
                Y=nM(row(i)+row2(1):row(i)+row2(size(row2,1))-1,channel);
                %reg=x\y;
                %Reg=[Reg reg];
                reg=[ones(length(x),1) x]\y;
                Reg=[Reg reg(2)];
                %yCalc=reg*x;
                %plot(x,yCalc);
                clear Deriv;
                figure(h2)

                for k=1:size(y)-sample
                    Deriv(k)=(y(k+sample)-y(k))/(sample*pdt);
                end
                MatDeriv(1:size(Deriv',1),i)=Deriv;
                MatDeriv(size(Deriv',1)+1:size(MatDeriv,1),i)=Deriv(k);
                %plot(x(1:size(x)-sample),Deriv,vcol(j));
                plot(x,y,vcol(j));%plot des raw data
                %pause(1)
            end
        end
    end

    % on ajoute à la main la dernière simu
%     [row2,col2]=find(nM(row(i+1):size(nM,1),channel)>=inf & nM(row(i+1):size(nM,1),channel)<=sup);
%            
%     if isempty(row2)
%         %Reg=[Reg -100];
%     else
%         MatTot(row2(1):row2(size(row2,1))-1,i+1)=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,channel);
%         MatTot(row2(size(row2,1)):size(MatTot,1),i+1)=nM(row(i+1)+row2(size(row2,1))-1,channel);
%         x=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,1);
%         y=nM(row(i+1)+row2(1):row(i+1)+row2(size(row2,1))-1,channel);
%         if max(y)>seuil
% 
%             
%             %reg=x\y; %lineaire
%             %Reg=[Reg reg];
%             reg=[ones(length(x),1) x]\y; %correction affine
%             Reg=[Reg reg(2)];
%             %yCalc=reg*x;
%             %plot(x,yCalc);
%             figure(h2)
%             clear Deriv;   
%             for k=1:size(y)-sample
%                 Deriv(k)=(y(k+sample)-y(k))/(sample*pdt);
%             end
%             MatDeriv(1:size(Deriv',1),i+1)=Deriv;
%             MatDeriv(size(Deriv',1)+1:size(MatDeriv,1),i+1)=Deriv(k);
%             %plot(x(1:size(x)-sample),Deriv,vcol(j)); 
%             plot(x,y,vcol(j));%plot des raw data
%             %pause(1)
%         end
%     end
%     MaxStock=max(MatDeriv);
    
    title(sprintf('v = %.1f',v));
    
    m1=mean(Reg);
    sd1=std(Reg);
    size(Reg,2)
    es1=sd1/sqrt(size(Reg,2));
    graphique(1,j)=m1;
    graphique(2,j)=es1;
    
    figure(h1)
     [row3,col3]=find(MatTot(size(MatTot,1),:)>0);%T(j)/3);
     M=mean(MatTot(:,col3)');Stock1M(1:size(M'),j)=M;
     plot(pdt:pdt:Eqtemps*pdt,M,vcol(j))
%     plot(pdt:pdt:Eqtemps*pdt,(N-M)./(pdt:pdt:Eqtemps*pdt).^0.829,vcol(j));%rescaled

    clear M
    figure(h3)
    M=mean(MatDeriv');
    Stock2M(1:size(M'),j)=M;
    MaxiM(j,1)=mean(MaxStock);
    MaxiM(j,2)=std(MaxStock);
    
    plot(pdt:pdt:Eqtemps*pdt,M,vcol(j));

    j=j+1;
end
figure(h1)
xlabel('Time')
ylabel('# Actives')
title('mean signal')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
%axis([40 160 2 140]);
leg=legend(num2str(vitesse'),'Location','southeast');
title(leg,'v');

figure(h2)
xlabel('Time')
ylabel('# Actives')
title('All signals')

figure(h3)
xlabel('Time')
ylabel('# Actives')
title('Derivative of nM')

figure
for i=1:size(vitesse,2)
    pos(i)=graphique(2,i);
    if graphique(1,i)-graphique(2,i)<0
        neg(i)=graphique(1,i);
    else
        neg(i)=graphique(2,i);
    end   
end
errorbar(abs(vitesse),graphique(1,:),neg(:),pos(:),'-x')
axis([0 50 0 70])
xlabel('Speed')
ylabel('Activation speed')
title('Mean slope of signals')

%Comparaison avec lipro classique
%vitmoy=dlmread("vitesse_moy.dat")
hold on 
errorbar(vitmoy(10:22,1),vitmoy(10:22,2),vitmoy(10:22,3),'-x')
legend({"lipro stop" "lipro classic"})

%%% Image final
cd /tmp/wkm5/Stops
%close all
N=2000;
rep=20;
v=-1;
simu=20;
i=simu-1;
Movers=dlmread(sprintf("data_Movers_STOP_N_%g_V_%.2f0000_rep_%g_vchc.dat",N,v,rep));

for i=0:0%rep-1
    %B=abs(Movers(i*N+1,1)-Movers(i*N+10000,1))
    %% plot de l'image finale
    for j=1:N %on actualise les positions
        if Movers(i*N+j,2)==0
            VectColor(j,1)= 1;VectColor(j,2)= 0;VectColor(j,3)= 0;
        else
            VectColor(j,1)= 0;VectColor(j,2)= 1;VectColor(j,3)= 0;
        end
    end
    figure
    scatter(Movers(i*N+1:(i+1)*N,1),zeros(2000,1),10,VectColor,'Filled');
    
%     %% scan de la densité en largeur
%     figure;
%     H=histogram(Movers(i*N+1:(i+1)*N,1),"BinEdges",Movers(i*N+1,1)+0.5:1:Movers(i*N+10000,1));%scanne le long de l'axe x
%     hvalue=H.Values;
%     close
%     figure
%     X=Movers(i*N+1,1)+0.5+0.5:1:Movers(i*N+10000,1)-1+0.5;
%     plot(X/abs(v),hvalue*abs(v))
%     for j=1:size(hvalue')-1
%         Diff(j)=hvalue(j+1)-hvalue(j);
%     end
%     %plot(Diff);
%     %axis([-10 110 -10 110]);% mettre pos(1) si on suit tous les movers
%     %% analyse des fluctuations
% %     bin=10;
% %     H=hist3([Movers(i*N+1:(i+1)*N,1),Movers(i*N+1:(i+1)*N,2)],'CdataMode','auto','Ctrs',{Movers(i*N+1,1)+0.5:bin:Movers(i*N+10000,1) (0.5:bin:99.5)},'LineStyle','none');
% %     h1=surf(Movers(i*N+1,1)+0.5:bin:Movers(i*N+10000,1),(0.5:bin:99.5),H')
% %     set(h1,'LineStyle','none');
% %     view(2); colorbar; 
% %     Z=h1.ZData;
end
axis([-1 2000 -1 1])

%%% nImiss
cd /tmp/wkm5/Stops
%close all
N=2000;
rep=20;
vitesse=[-50 -10 -5:1:-1];
v=-1;
simu=20;
i=simu-1;
j=1

for v=vitesse
   v
    pause(0.0001)
    nM=dlmread(sprintf("data_nM_STOP_N_%g_V_%.2f0000_rep_%g_vchc.dat",N,v,rep));

    [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu

    %analyse des propagations
    [row2,col2]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
    clear actifin
    for i=2:size(row2,1)
        actifin(i-1)=nM(row2(i)-1,2)*100/N; % à la différence du liproclassic, les "actifs" sont ceux qui ont fait l'arrêt. La definition étendue d'actifs sur les deux modèles est donc ceux qui ont fait un arrêt
    end
    m2=mean(actifin);
    sd2=std(actifin);
[row,col]=find(nM(:,1)==0);
for i=2:size(row,1)
    diff(i-1)=row(i)-row(i-1);
end
D=min(diff); %calcul du nomber de pas de temps min d'une simu

    graphique2(1,j)=m2;
    graphique2(2,j)=sd2;
    
    k=1;
    vec=[(row(2:size(row,1),1)-1)' size(nM,1)];
    for i=vec
        bp(k,j)=nM(i,2);
        k=k+1;
    end
    j=j+1;
end
% % 

figure
for i=1:size(vitesse,2)
    if graphique2(1,i)+graphique2(2,i)>100
        pos(i)=100-graphique2(1,i);
    else
        pos(i)=graphique2(2,i);
    end
    if graphique2(1,i)-graphique2(2,i)<0
        neg(i)=graphique2(1,i);
    else
        neg(i)=graphique2(2,i);
    end   
end

errorbar(abs(vitesse),graphique2(1,:),neg(:),pos(:),'-x')
axis([-1 51 97 100.1]);

ylabel({"Nombre d'individus" "propagés en moyenne"});
title("Compilation des moyennes d'individus propagés");
savename="compil_finalid_v_0_50_zoom_urescaled"


%% Lipro full model cohesion

cd /tmp/wkm5/FDA-STOPS/2-state/2particules/
cd /tmp/wkm5/FDA-STOPS/2-state/seeingAll
%cd /tmp/wkm5/FDA-STOPS/2-state/'T=10000 w=1'/
clear all;
close all;
cd /tmp/wkm10/CollabLuis/2part/'indep S'/asym/

N=2;
rep=500;
v=-1;
seuil=20;
perc=2;
sim=20;
T=1000;

%K=15301;
%Cohe=dlmread(sprintf("data_Cohesion_FULL_N_%g_V_%.2f0000_rep_%g_vc_a0.dat",N,v,rep));
Cohe=dlmread(sprintf("data_Cohesion_FULL_N_%g_param_1.dat",N));
Npdt=size(Cohe,1)/rep;
% %[row,col]=find(Cohe(:,1)==0);
% %Npdt=K;
% %rep=400
% 
% %rep simu individuelles
% % figure;hold on;
% % for i=1:rep
% %     plot(Cohe(1+(i-1)*Npdt:Npdt+(i-1)*Npdt,1),Cohe(1+(i-1)*Npdt:Npdt+(i-1)*Npdt,2));
% % end
% % title("all simulations")
% % ylabel('Cohesion Coefficient');
% % xlabel('Time');
% % print("Cohesion_allsim_w_01_a_1", '-dpdf','-painters','-bestfit')
% 
% 
% 
% figure;hold on;
% for i=1:rep
%     plot(Cohe(1+(i-1)*Npdt:Npdt+(i-1)*Npdt,1),Cohe(1+(i-1)*Npdt:Npdt+(i-1)*Npdt,2)-1);
% end
% title("delta (= Delta-Delta_0) function of time (all sim)")
% ylabel('delta');
% xlabel('Time');
% print("Cohesion_delta_2p_SIndep", '-dpdf','-painters','-bestfit')
% 
% 
% %moyenne des rep simus
% % figure;hold on;
% % CoheMoy=zeros(1,Npdt);
% % for i=1:rep
% %     for j=1:Npdt
% %         CoheMoy(j)=CoheMoy(j)+Cohe(j+(i-1)*Npdt,2);
% %     end
% % end
% % plot(Cohe(1:Npdt,1),CoheMoy(:)./rep);
% % title('Mean over all simulations')
% % ylabel('Cohesion Coefficient');
% % xlabel('Time');
% % print("Cohesion_meansim_w_01_a_0", '-dpdf','-painters','-bestfit')
% 
% %%%histogramme des distances à temps t
% % CoheHist=zeros(1,rep);
% % t=15301;%temps final 15301
% % for i=1:rep
% %      CoheHist(i)=Cohe(t+(i-1)*Npdt,2);
% % end
% % figure;
% % hist(CoheHist(:));
% % title("Histogram of final Cohesion coeff")
% % xlabel('Final Coheficient')
% % print("Cohesion_finCoheHist_w_01_a_1_THHV", '-dpdf','-painters','-bestfit')
% % 
% % figure;hold on;
CoheMoy=zeros(1,Npdt);
Cohe2Moy=zeros(1,Npdt);
for i=1:rep
    for j=1:Npdt
        CoheMoy(j)=CoheMoy(j)+Cohe(j+(i-1)*Npdt,2);
        Cohe2Moy(j)=Cohe2Moy(j)+Cohe(j+(i-1)*Npdt,3);
    end
end

% plot(Cohe(1:Npdt,1),CoheMoy(:)./rep);
% title('Mean over all simulations')
% ylabel('Cohesion Coefficient');
% xlabel('Time');
%print("Cohesion_meansim_w_01_a_0_2", '-dpdf','-painters','-bestfit')

% figure
% plot(Cohe(1:Npdt,1),Cohe2Moy(:)./rep);
% title('Cohesion coefficient squared over time')
% ylabel('Cohesion Coefficient squared');
% xlabel('Time');
%print("Cohesion2_meansim_w_01_a_0_2", '-dpdf','-painters','-bestfit')

figure
plot(Cohe(1:Npdt,1),Cohe2Moy(:)./rep-(CoheMoy(:)./rep).^2);
title('Variance of Delta over time (mean sim)')
ylabel('Variance of Delta');
xlabel('Time');
print("Cohesion3_Varmeansim_2p_a0_vdiff", '-dpdf','-painters','-bestfit')

%% COLLAB LUIS
cd /tmp/wkm10/CollabLuis/2part/asym/'extension 10P'/

Vmax=100;
Vmin=0;
%param=[Vmin+abs(Vmax/20):abs(Vmax/20):Vmax];
param=[Vmin+abs(Vmax-Vmin)/20:abs(Vmax-Vmin)/20:Vmax];
%figure('units','normalized','outerposition',[0 0 1 1])
for i=1:1
    A=dlmread(sprintf('data_Variance_Cohesion_N_2_param_%g.dat',i));
    %subplot(4,5,i)
    plot(A(:,1),A(:,2))
    %title(param(i))
    %title(i)
    ylabel('\sigma')
    xlabel('T')
end
orient landscape
print("variance_tendency_allRangeParam", '-dpdf','-painters','-bestfit') %range param
print("variance_tendency_asym_10P", '-dpdf','-painters','-bestfit') %sans range

%%% support simu numerical integration Luis
%cd /tmp/wkm5/FDA-STOPS/2-state/2particules
cd /tmp/wkm10/CollabLuis/2part/'indep S'/

N=2;
rep=500;
v=-1;
T=1000;

%Cohe=dlmread(sprintf("data_Cohesion_FULL_N_%g_V_%.2f0000_rep_%g_vchc_2_a1.dat",N,v,rep));
Cohe=dlmread(sprintf("data_Cohesion_FULL_N_%g_param_1.dat",N));

pas=100;
Npdt=size(Cohe,1)/rep;
Scan=[1:Npdt:size(Cohe,1)-Npdt];
hvalue=zeros((Npdt-1)/pas,20);
for i=1:pas:Npdt-1
    figure;
    H=histogram(abs(Cohe(Scan+(i-1),2)),[0:20]);%scanne le long de l'axe x
    hvalue((i-1)/pas+1,1:20)=H.Values;
    close
end
    figure
    X=Cohe(1:pas:Npdt-1,1);
    Y=[0:19];
%     h1=surf(X,Y',hvalue')
    mesh(X,Y',hvalue'/rep)
%     set(h1,'LineStyle','none');% colorbar
%view(2)%vue 2D
az = 225;
el = 20;
view(az, el);%ou3D
xlabel('Time')
ylabel('Delta')
%zlabel("Freq")
colorbar
axis([0 100 0 20])
%print("probability_SIndep", '-dpdf','-painters','-bestfit')

%%% extreme value stat
cd /tmp/wkm10/CollabLuis/2part/Slope/exp/extremeVS/

j=10
figure
hold on
for i=j:j+2
    A=dlmread(sprintf('data_Variance_Cohesion_N_2_param_%g.dat',i));
    %subplot(4,5,i)
    plot(A(:,1),A(:,2))
    M(i-j+1)=max(A(:,2))
    %title(param(i))
    %title(i)
    ylabel('\sigma')
    xlabel('T')
end
scatter([1000 5000 10000],M)
legend('T=1000','T=5000','T=10000','maxima')
orient landscape
print("extreme_val_stat_asym25", '-dpdf','-painters','-bestfit')

%%% switching behaviour
clear all;
close all;
cd /tmp/wkm10/CollabLuis/2part/asym/'extension 10P'

N=2; %pas adapté pour >2P
rep=500;
v=-1;
seuil=20;
perc=2;
sim=20;
T=1000;

%Cohe=dlmread(sprintf("data_Cohesion_FULL_N_%g_V_%.2f0000_rep_%g_vc_a0.dat",N,v,rep));
Cohe=dlmread(sprintf("data_Cohesion_FULL_N_%g_param_1.dat",N));
Npdt=size(Cohe,1)/rep;

SiCohe= [Cohe(:,1) sign(Cohe(:,2))];
    %find(SiCohe(:,2)==0)
    j=1;
    clear MatSwitch
for i=1:size(Cohe,1)-1
    if mod(i,Npdt)==0
        
    else
        if SiCohe(i,2)~=SiCohe(i+1,2)
        MatSwitch(j)=SiCohe(i,1);
        j=j+1;
        end
    end
end
h=histogram(MatSwitch,[0:5:1000])
hvalue=h.Values;
plot([5:5:1000],hvalue/rep)
xlabel('Time')
ylabel('# of switch')
print("numberOfSwitch_T1000_asym_10P", '-dpdf','-painters','-bestfit')

%%% first passage time
clear all;
close all;
cd /tmp/wkm10/CollabLuis/2part/'indep S'/asym/FirstPassageTime/

N=2;
rep=500;
v=-1;
seuil=20;
perc=2;
sim=20;
T=1000;

%Cohe=dlmread(sprintf("data_Cohesion_FULL_N_%g_V_%.2f0000_rep_%g_vc_a0.dat",N,v,rep));
FPT=dlmread(sprintf("data_FirstPT_FULL_N_%g_param_5.dat",N));
%Npdt=size(Cohe,1)/rep;

histogram(FPT,[0:100:5000])
print("hist_firstpassagetime_D20", '-dpdf','-painters','-bestfit')

%Compil
h1=figure
for i=1:5
    FPT=dlmread(sprintf("data_FirstPT_FULL_N_%g_param_%g.dat",N,i));
    figure;
    h=histogram(FPT,20)
    hvalue1=h.BinEdges(1:h.NumBins)+h.BinWidth/2;
    hvalue2=h.Values;
    close;
    figure(h1)
    hold on
    plot(hvalue1,hvalue2)
end
xlabel('First passage time')
ylabel('Frequence')
legend({'D=1','D=5','D=10','D=20','D=25'})
set(gca,'YScale','log')
set(gca,'XScale','log')
print("compil_firstpassagetime", '-dpdf','-painters','-bestfit')

%%% barycentre
clear all;
close all;
cd /tmp/wkm10/CollabLuis/2part/Slope/exp/Movers/

N=10;
rep=500;
v=-1;
seuil=20;
perc=2;
sim=20;
T=1000;

%Cohe=dlmread(sprintf("data_Cohesion_FULL_N_%g_V_%.2f0000_rep_%g_vc_a0.dat",N,v,rep));
barycentre=dlmread(sprintf("data_barycentre_FULL_N_%g_param_1.dat",N));
x=barycentre(:,1);
y=-(barycentre(:,2));
plot(x,y)
reg=x\y
xlabel('Time')
ylabel('pos barycentre')
print("dyn_barycentre_10P", '-dpdf','-painters','-bestfit')

%comparaison dyn bary
hold on
N=2
barycentre2=dlmread(sprintf("data_barycentre_FULL_N_%g_param_1.dat",N));
x2=barycentre2(:,1);
y2=-(barycentre2(:,2));
plot(x2,y2)
reg=x2\y2
print("dyn_barycentre_comp2P10P", '-dpdf','-painters','-bestfit')

%%% single sim cohesion
i=1
A=dlmread(sprintf('data_Cohesion_FULL_N_%g_param_%g.dat',N,i));
plot(A(:,1),A(:,2))
xlabel('Time')
ylabel('Delta')
print("dyn_delta_1rep_10P", '-dpdf','-painters','-bestfit')

%%% mesure de synchronicité
clear all;
close all;
cd /tmp/wkm10/CollabLuis/2part/'indep S'/exp/compRegLin/
%cd /tmp/wkm10/CollabLuis/2part/Slope/exp/compRegLin/

N=2;
rep=500;
v=-1;
seuil=20;
perc=2;
sim=20;
T=1000;

nM=load(sprintf('data_nM_FULL_N_%g_param_1.dat',N));

% figure;
% hold on;
% 
% for i=1:1
%     plot(nM((i-1)*size(nM,1)/rep+1:(i-1)*size(nM,1)/rep+size(nM,1)/rep,1),nM((i-1)*size(nM,1)/rep+1:(i-1)*size(nM,1)/rep+size(nM,1)/rep,2)/N*100)
% end

%%% mesure de temps de synchro
[row,col]=find(nM(:,2:4)==N);
TimePercTot=size(nM(row,:),1)/size(nM(:,:),1)*100
[row2,col2]=find(nM(:,2)==N);
TimePercStop=size(nM(row2,:),1)/size(nM(:,:),1)*100
[row3,col3]=find(nM(:,4)==N);
TimePercMove=size(nM(row3,:),1)/size(nM(:,:),1)*100
