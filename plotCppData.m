function [legendCell, Ec_times_beta]=plotCppData(fileName,fig)
%close all; clc;
set(0,'DefaultLegendFontSize',12,'DefaultLegendFontSizeMode','manual');
pathName =strcat('Project/data/',fileName);
lengthAlphaVec = dlmread(pathName,' ',[2 3 2 3]);
lengthKvec = dlmread(pathName,' ',[2 4 2 4]);
Kvec=dlmread(pathName,' ',[0 0 0 lengthKvec-1]);
Ec_times_beta = dlmread(pathName,' ',[0 lengthKvec 0 lengthKvec]);
Lx = dlmread(pathName,' ',[0 lengthKvec+1 0 lengthKvec+1]);
alphaVec=dlmread(pathName,' ',[1 0 1 lengthAlphaVec-1]);
doubleAlphaVec =2*alphaVec;
%doubleAlphaVec =8*alphaVec; %??????
Nstat = dlmread(pathName,' ',[2 0 2 0]);
Nwarmup = dlmread(pathName,' ',[2 1 2 1])
Nprod = dlmread(pathName,' ',[2 2 2 2]);
alpha_lower = dlmread(pathName,' ',[2 5 2 5]);
alpha_upper = dlmread(pathName,' ',[2 6 2 6]);
extra_Nprod_factor = dlmread(pathName,' ',[2 7 2 7]);
factor = dlmread(pathName,' ',[2 8 2 8]);
lambda = dlmread(pathName,' ',[2 9 2 9]);
%mark={'.b-','*r-','^k-','dg-','.m-','*b-','^r-','.k:','.g:','*m-'};
marker={'.','*','^',};
orange = {1, 0.5, 0};
blue = [0 0.5 1];
%color={'b','r','k','g','m','c','y','o','orange','blue'};
%color = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
color = {[0 1 1], [1 0 1], [1 0 0], [0 1 0], [0 0 1], [0 0 0], [0.91 0.41 0.17], [255/255 153/255 255/255], [102/255 0/255 51/255], [153/255 255/255 0/255]};
if fig ~=0
     figure(fig);
end
   
for i=0:(length(Kvec)-1)
    K=Kvec(i+1);
    meanCk = (Lx-1)*dlmread(pathName,' ',[3+length(doubleAlphaVec)*i Nstat 3+length(doubleAlphaVec)*(i+1)-1 Nstat]);
    STDCk = (Lx-1)*dlmread(pathName,' ',[3+length(doubleAlphaVec)*i Nstat+1 3+length(doubleAlphaVec)*(i+1)-1 Nstat+1]);
    if K==0.1
        STDCk
    end
    %e=errorbar(-doubleAlphaVec*log(besseli(1,K)/besseli(0,K)),meanCk,STDCk);
    e=errorbar(doubleAlphaVec,meanCk,STDCk);
    %e.Marker = NONE;
    %e.MarkerSize = 10;
    color{i+1};
    e.Color = color{i+1};
    hold on;

    %e.CapSize = 15;
    %plot(alphaVec,meanCk,*);
    hold on;
end

set(0,'defaulttextinterpreter','latex');
xlabel('$2R/R_{Q}$', 'FontSize', 20);

legendCell=strcat('K=',strtrim(cellstr(num2str(Kvec(:)))));
%if fig ~= 0
    [h, ~,plots]=legend(legendCell,'Location','northwest');
%end
%[h, ~, plots] = legend('sin(x)', 'cos(x)');
for i=1:length(Kvec)
    K=Kvec(i);
    for j=1:length(doubleAlphaVec)
        meanCk = (Lx-1)*dlmread(pathName,' ',[3+length(doubleAlphaVec)*(i-1) Nstat 3+length(doubleAlphaVec)*(i)-1 Nstat]);
        alpha = doubleAlphaVec(j);
        if alpha < alpha_lower || alpha > alpha_upper
            plot(doubleAlphaVec(j),meanCk(j),'*','color',color{i});
            %plot(doubleAlphaVec(j)/(-1/K*log(besseli(0,K)/besseli(1,K))),meanCk(j),strcat('*',color{i}));
        else
            plot(doubleAlphaVec(j),meanCk(j),'.','color',color{i});
            %plot(doubleAlphaVec(j)/(log(besseli(1,K)/besseli(0,K))),meanCk(j),strcat('.',color{i}));
        end
    end 
end

if fig ~=0
    title(strcat('$\lambda=$',num2str(lambda),'$, L_x=$',num2str(Lx),'$, \beta=$',num2str(Ec_times_beta),'$/E_c$'));
    ylabel('$\tilde{C}_{k=0}^{-1}$','FontSize', 20);
end
set(h,'FontSize',20);

%set(gcf, 'PaperPosition', [0 0 15 10]); %Position plot at left hand corner with width 5 and height 5.
%set(gcf, 'PaperSize', [15 10]); %Set the paper to have width 5 and height 5.
%saveas(gcf, 'test', 'pdf') %Save figure
