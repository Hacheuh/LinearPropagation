clear all;
%clc;
close all;

%% SPLITTING
% figure; hold on;
% leg={};
% j=1;
% for v=-1%:1:1
%     x1=0:0.01:10;
%     ycalc1=exp(2/abs(v)*exp(-1/2)*(exp(-x1.*abs(v)/2)-1));
%     plot(x1,ycalc1);
%     leg=[leg num2str(v)];
%     j=j+1;     
% end
% xlabel('Temps');
% ylabel("p_0");
% axis([0 10 0 1]);
% %legend(leg);
% %title({"Probabilité de trouver la particule inactive au cours du temps" "pour plusieurs valeurs de vitesse"});
% saveName = ('fig_p0_2part_split'); 
% 
% %saveas(gcf, saveName, 'png'); 


%% OVERLAP

% figure; hold on;
% leg={};
% j=1;
% for v=1%:1:1
%     x1=0:0.01:abs(1/v);
%     x2=abs(1/v):0.01:10;
%     ycalc1=exp(-2/abs(v)*exp(-1/2)*(exp(x1.*abs(v)/2)-1));
%     ycalc2=exp(-2/abs(v)*exp(-1/2)*(exp(1/2)-1))*exp(2/abs(v)*exp(-1/2)*(exp(-((x2-abs(1/v))*abs(v))/2)-1));
%     plot([x1 x2],[ycalc1 ycalc2]);
%     temp_scat(1,j)= x1(int32(abs(round(1/v,2)/0.01)));
%     temp_scat(2,j)=ycalc1(int32(abs(round(1/v,2)/0.01)));
%     leg=[leg num2str(v)];
%     j=j+1;     
% end
% xlabel('Temps');
% ylabel("p_0");
% %legend(leg);
% %scatter(temp_scat(1,:),temp_scat(2,:),"ro");
% %title({"Probabilité de trouver la particule inactive au cours du temps" "pour plusieurs valeurs de vitesse"});
% saveName = ('fig_p0_2part_overlap');
% axis([0 10 0 1]);



% figure; hold on;
% leg={};
% j=1;
% for v=-1%:1:1
%     x1=0:0.01:10;
%     ycalc1=exp(2/abs(v)*exp(-1/2)*(exp(-x1.*abs(v)/2)-1));
%     plot(x1,ycalc1);
%     leg=[leg num2str(v)];
%     j=j+1;     
% end
% xlabel('Temps');
% ylabel("p_0");
% axis([0 10 0 1]);
% %legend(leg);
% %title({"Probabilité de trouver la particule inactive au cours du temps" "pour plusieurs valeurs de vitesse"});
% saveName = ('fig_p0_2part_split'); 

% %saveas(gcf, saveName, 'png'); 



%%%% dyn du Taux de transition par unité de temps 2 part
% v=0.05;
% x1=0:0.01:abs(1/v);
% x2=abs(1/v):0.01:100;
% x3=0:0.01:100;
% ycalc1=0.1*exp(-(1-abs(v).*x1)/2);
% ycalc2=0.1*exp(-(1+abs(v).*(x2-2/abs(v)))/2);
% ycalc3=0.1*exp(-(1+abs(v).*x3)/2);
% plot([x1 x2],[ycalc1 ycalc2],x3, ycalc3);

%%% dyn du Taux de transition par unité de temps n part dans le cas
%%% dépassement
% clear all
% v=1; 
% d=1;
% N=1; % nombre de particules inactives
% ycalc0=0.1;
% for i=1:N+1 %N+1 intervalles pour N particules inactives
%     %i
%     sumexp=0;
%     if i==N+1
%         %eval(sprintf('t%g=(i-1)*d/abs(v):0.01:2*i*abs(d/v);',i)) %% dernier intervalle de temps on ajoute le double de l'intervalle en guise de marge
%         eval(sprintf('t%g=(i-1)*d/abs(v):0.01:10;',i))
%     else
%         eval(sprintf('t%g=(i-1)*d/abs(v):0.01:i*abs(d/v);',i))
%     end
%     for j1=i:N
%         %j1
%         eval(sprintf('sumexp=sumexp+exp(-(j1*d-abs(v).*t%g)/2);',i)) %rapprochement
%     end
%     for j2=1:i-1
%         %j2
%         eval(sprintf('sumexp=sumexp+exp(-(j2*d+abs(v).*(t%g-2*j2*d/abs(v)))/2);',i)) %eloi
%     end
%     eval(sprintf('ycalc%g=0.1*sumexp;',i))
% end
% T=[];
% Y=[];
% for i=1:N+1
%     eval(sprintf('T=[T t%g];',i)) 
%     eval(sprintf('Y=[Y ycalc%g];',i))
% end
% plot(T,Y);
% saveName = ('fig_Gamma_2part_overlap'); 
% axis([0 10 0 0.11]);
% saveas(gcf, saveName, 'png'); 

%%%%% split
% v=1; 
% d=1;
% t=0:0.01:10;
% sumexp=exp(-(d+abs(v).*t)/2);
% ycalc=0.1*sumexp;
% plot(t,ycalc);
% saveName = ('fig_Gamma_2part_split'); 
% axis([0 10 0 0.11]);
% saveas(gcf, saveName, 'png'); 


%% dynamique du système étendu

%Na=1 N=K
figure; hold on;
d=2;
leg={};
j=1;
for v=1%:1:1
    for K=10:10
        sum=0;
        for l=0:(K-2)
            sum=sum+exp(l*1/2);
        end
        t=0:0.01:10;
        ycalc1=exp(2/abs(v)*exp(-(K-1)*1/2)*(sum)*(exp(-t.*abs(v)/2)-1));
        plot(t,ycalc1);
        leg=[leg num2str(v)];
        j=j+1;
    end
end
xlabel('Time');
ylabel("p");

axis([0 10 0 1]);
%legend(leg);
title({"Probability of existence for a 1A system"});
saveName = ('fig_p_1AKN'); 

saveas(gcf, saveName, 'png'); 

%% Probabilité d'existence des systèmes 2A
figure;
t=0;v=-1;d=2;
x=1:1:10;
G_tot=exp(-t*v/d)/(exp(1)^1/d-1);
y=exp(-(x-t*v)/d)/G_tot;
plot(x,y,'+')
title({"Probability of having the transition from 1 to x at t=0"});
xlabel('ID');
ylabel("p");
% saveName = ('fig_p_trans1A2A'); 
% saveas(gcf, saveName, 'png'); 