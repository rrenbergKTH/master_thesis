function [legendCell, Ec_times_beta]=plotCppData(fileName,fig)

set(0,'DefaultLegendFontSize',12,'DefaultLegendFontSizeMode','manual','defaulttextinterpreter','latex');
pathName =strcat('data/',fileName);

% Read file-data
lengthAlphaVec = dlmread(pathName,' ',[2 3 2 3]);
lengthKvec = dlmread(pathName,' ',[2 4 2 4]);
Kvec=dlmread(pathName,' ',[0 0 0 lengthKvec-1]);
Ec_times_beta = dlmread(pathName,' ',[0 lengthKvec 0 lengthKvec]);
Lx = dlmread(pathName,' ',[0 lengthKvec+1 0 lengthKvec+1]);
alphaVec=dlmread(pathName,' ',[1 0 1 lengthAlphaVec-1]);
doubleAlphaVec =2*alphaVec;
Nstat = dlmread(pathName,' ',[2 0 2 0]);
Nwarmup = dlmread(pathName,' ',[2 1 2 1])
Nprod = dlmread(pathName,' ',[2 2 2 2]);
alpha_lower = dlmread(pathName,' ',[2 5 2 5]);
alpha_upper = dlmread(pathName,' ',[2 6 2 6]);
extra_Nprod_factor = dlmread(pathName,' ',[2 7 2 7]);
factor = dlmread(pathName,' ',[2 8 2 8]);
lambda = dlmread(pathName,' ',[2 9 2 9]);


color_vec = {[255/255 153/255 255/255] [0 1 1], [1 0 1], [1 0 0], [0 1 0], [0 0 1],[0 0 0],[0 0.5 1],[1 0.5 0],[153/255 255/255 0/255]};

key_set = [0.05 0.1, 0.2, 0.4, 0.8, 1.5, 2.5, 3,6,9];
color_map = containers.Map(key_set,color_vec)

if fig ~=0
     figure(fig);
end

%Plot error bar
for i=0:(length(Kvec)-1)
    K=Kvec(i+1);
    meanCk = (Lx-1)*dlmread(pathName,' ',[3+length(doubleAlphaVec)*i Nstat 3+length(doubleAlphaVec)*(i+1)-1 Nstat]);
    STDCk = (Lx-1)*dlmread(pathName,' ',[3+length(doubleAlphaVec)*i Nstat+1 3+length(doubleAlphaVec)*(i+1)-1 Nstat+1]);
    e=errorbar(doubleAlphaVec,meanCk,STDCk);
    e.Color = color_map(K);
    hold on;
end

xlabel('$2R/R_{Q}$', 'FontSize', 20);

legendCell=strcat('K=',strtrim(cellstr(num2str(Kvec(:)))));
[h, ~,~]=legend(legendCell,'Location','northwest');

%plot an extra marker
for i=1:length(Kvec)
    K=Kvec(i);
    for j=1:length(doubleAlphaVec)
        meanCk = (Lx-1)*dlmread(pathName,' ',[3+length(doubleAlphaVec)*(i-1) Nstat 3+length(doubleAlphaVec)*(i)-1 Nstat]);
        alpha = doubleAlphaVec(j);
        if alpha < alpha_lower || alpha > alpha_upper
            plot(doubleAlphaVec(j),meanCk(j),'*','color',color_map(K));
        else
            plot(doubleAlphaVec(j),meanCk(j),'.','color',color_map(K));
        end
    end 
end

if fig ~=0
    title(strcat('$\lambda=$',num2str(lambda),'$, L_x=$',num2str(Lx),'$, \beta=$',num2str(Ec_times_beta),'$/E_c$'));
    ylabel('$\tilde{C}_{k=0}^{-1}$','FontSize', 20);
end

set(h,'FontSize',20);
