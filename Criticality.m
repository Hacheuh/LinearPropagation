%% lipro
cd /tmp/wkm5/DoubleS
cd /home/gascuel/Desktop/Criticality

N=2000;
rep=1000;
b=100;
vit=(-1.22:-0.02:-1.38);
j=1;
for v=vit
    v
    nM=dlmread(sprintf("data_nM_lipro_N_%g_V_%.2f0000_rep_%g_perc_2_vc.dat",N,v,rep));
    [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
    
    k=1;
    vec=[(row(2:size(row,1),1)-1)' size(nM,1)];
    for i=vec
        bp(k,j)=nM(i,2);
        k=k+1;
    end
    
    buffer=bp(:,j);
%     figure
%     h=histogram(eval(sprintf('bp(:,%g)',j)),b,'facealpha',0.3);
    %save(sprintf('histo_lipro_v_%.2f.txt',v), 'buffer', '-ASCII')
    
    j=j+1;
end



% h1=figure
% hold on
% h2=figure
% ax=gca;
% ax.YScale='log';
% hold on
% for i=(2:1:4)
%     Y=load(sprintf('histo_v1p%g.txt',i));
%     figure(h1)
%     plot(Y)
%     figure(h2)
%     semilogy(Y)
% end
% for i=(22:2:38)
%     Y=load(sprintf('histo_v1p%g.txt',i));
%     figure(h1)
%     plot(Y)
%     figure(h2)
%     semilogy(Y)
% end

%% SIR
% cd /tmp/wkm5/SIR
% h1=figure;hold on;
% for i=0:0.2:2
% A=load(sprintf('data_Finalid_SIR_N_2000_beta_%.1f00000_rep_1000',round(i,1)));
% b=100;
% clear hval
% figure;
% h=histogram(A,'facealpha',0.3,'BinEdges',(0:10:1990));
% hval=h.Values;
% hval=transpose(hval);
% %save(sprintf('histo_v%.1f.txt',i), 'hval', '-ASCII')
% figure(h1)
% plot(hval)
% end


%%%Liprostat pour SIR
cd /tmp/wkm5/SIR
%clear all;
%close all;

N=2000;
rep=1000;
v=2;
beta=v/2;
seuil=0;

FrontVague=dlmread(sprintf("data_FrontVague_SIR_N_%g_beta_%.1f00000_rep_%g",N,beta,rep));
nM=dlmread(sprintf("data_nM_SIR_N_%g_beta_%.1f00000_rep_%g",N,beta,rep));

figure;
hold on;
% h2=figure;
% hold on;
% h3=figure;
% hold on;
%%%Attention le changement de figure prend du temps
Reg=[];
[row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
for i=1:size(row,1)-1
    x=nM(row(i):row(i+1)-1,1);
    y1=nM(row(i):row(i+1)-1,2);
    y2=nM(row(i):row(i+1)-1,3);
    y3=nM(row(i):row(i+1)-1,4);
    if max(y1)>seuil
        reg=x\y2;
        Reg=[Reg reg];
        yCalc=reg*x;
        %plot(x,yCalc);%plot des regressions linéaires
        %figure(h1)
        plot(x,y1);%plot des raw data
%         figure(h2)
         %plot(x,y2);%plot des raw data
%         figure(h3)
%         plot(x,y3);%plot des raw data
    end
end
xlabel('Time');
ylabel("Number of active individuals");
title('Progression of the state Active for several simulations (nM)');
% saveName = (['fig_nMCompil_v_', num2str(v)]); 
% saveas(gcf, saveName, 'png'); 

%%% Analyse du faisceau
clear diff
h1=figure; hold on;
[row,col]=find(nM(:,1)==0);
for i=2:size(row,1)
    diff(i-1)=row(i)-row(i-1);
end
diff(size(row,1))=size(nM,1)-row(i);%on ajoute la dernière simu qui ne peut pas être calculée autrement
Dmin=min(diff); %calcul du nombre de pas de temps min d'une simu
Dmax=max(diff); %si on veut utiliser le filtrage quand splitting fort
figure
hD=histogram(diff);
hDval=hD.Values;
hDedges=hD.BinEdges;
close
Dseuil=hDedges(max(find(hDval>10))+1); % pour atteindre un seuil précis

for i=1:10:Dseuil
    %max(nM(find(nM(:,1)==0)+i)) %affichage du pas de temps
    
    data=nM(row(find(diff>i))+i,2);
    figure;
    h=histogram(data,10);
    hval=h.Values;
    close
    figure(h1)
    plot(hval)
    pause(0.001)
end

%%% Analyse des temps finaux
cd /tmp/wkm5/SIR
%clear all;
%close all;

N=2000;
rep=1000;
v=2;
beta=v/2;
seuil=0;

% Finalid=dlmread(sprintf("data_Finalid_SIR_N_%g_beta_%g.000000_rep_%g",N,beta,rep));
Mass=dlmread(sprintf("data_Mass_SIR_N_%g_beta_%g.000000_rep_%g",N,beta,rep));

% figure
% h=histogram(Finalid,'BinEdges',(0:10:140))
figure
h2=histogram(Mass,'BinEdges',(0:10:140))

%figure;hold on;
b=100;

Beta=(0:0.2:2);
j=1;
for beta=Beta
    beta
    Mass=dlmread(sprintf("data_Mass_SIR_N_%g_beta_%.1f00000_rep_%g",N,beta,rep));
    nM=dlmread(sprintf("data_nM_SIR_N_%g_beta_%.1f00000_rep_%g",N,beta,rep));
    [row,col]=find(nM(:,1)==0); %on récupère toute les lignes de début de simu
    
    k=1;
    vec=[(row(2:size(row,1),1)-1)' size(nM,1)];
    for i=vec
        bp(k,j)=nM(i,2);
        k=k+1;
    end
    
    buffer=bp(:,j);
    figure
    h=histogram(eval(sprintf('bp(:,%g)',j)),b,'facealpha',0.3);
    save(sprintf('histo_SIR_beta_%.1f.txt',beta), 'buffer', '-ASCII')
    
    j=j+1;
end

X=[];
Y=[];
vect=1:5;
for i=vect
    X=[X Beta(i)*ones(1,1000)];
    Y=[Y bp(:,i)'];
end    
hist3([X',Y'],'CdataMode','auto','FaceColor','interp','Ctrs',{Beta(vect) 50:100:1950});
colorbar
az = 300;
el = 20;
view(az, el);
ylabel("Final active pop.");
xlabel("Beta");
zlabel("# of simulations");
saveas(gcf,'fig_distrib_nAfin_compil_hist3D_zoom_SIR','png')