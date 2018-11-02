function [legendCell, Ec_times_beta]=plotCppData(fileName, fig, color_type, fontsize_legend)

set(0,'DefaultLegendFontSize',fontsize_legend,'DefaultLegendFontSizeMode','manual');
pathName =strcat('data/',fileName);

% Read file
lengthAlphaVec = dlmread(pathName,' ',[2 3 2 3]);
lengthKvec = dlmread(pathName,' ',[2 4 2 4]);
Kvec=dlmread(pathName,' ',[0 0 0 lengthKvec-1]);
Ec_times_beta = dlmread(pathName,' ',[0 lengthKvec 0 lengthKvec]);
Lx = dlmread(pathName,' ',[0 lengthKvec+1 0 lengthKvec+1]);
alphaVec=dlmread(pathName,' ',[1 0 1 lengthAlphaVec-1]);
doubleAlphaVec =2*alphaVec;
Nstat = dlmread(pathName,' ',[2 0 2 0]);
Nwarmup = dlmread(pathName,' ',[2 1 2 1]);
Nprod = dlmread(pathName,' ',[2 2 2 2]);
alpha_lower = dlmread(pathName,' ',[2 5 2 5]);
alpha_upper = dlmread(pathName,' ',[2 6 2 6]);
extra_Nprod_factor = dlmread(pathName,' ',[2 7 2 7]);
factor = dlmread(pathName,' ',[2 8 2 8]);
lambda = dlmread(pathName,' ',[2 9 2 9]);

% Define colors
blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
yellow = [0.9290, 0.6940, 0.1250];
black = [0,0,0];
purple = [0.4940, 0.1840, 0.5560];
green = [0.4660, 0.6740, 0.1880];
lightblue = [0.3010, 0.7450, 0.9330];
darkred = [0.6350, 0.0780, 0.1840];
magenta = [0.75, 0, 0.75];
cyan = [0, 0.75, 0.75];

% Create color map
if color_type == 1
    key_set = [0.1, 0.2, 0.4, 0.8, 1.5,3,6,9];
    color = {blue, orange, magenta, yellow, black, green, darkred, cyan};
else
    key_set = [0.2 0.4 0.8 1.5 2.5 3.5];
    color = {blue, orange, magenta, yellow, black, green};
end

M = containers.Map(key_set,color);

if fig ~=0
     figure(fig);
end

% Plot errorbar
for i=0:(length(Kvec)-1)
    K=Kvec(i+1);
    meanCk = (Lx-1)*dlmread(pathName,' ',[3+length(doubleAlphaVec)*i Nstat 3+length(doubleAlphaVec)*(i+1)-1 Nstat]);
    STDCk = (Lx-1)*dlmread(pathName,' ',[3+length(doubleAlphaVec)*i Nstat+1 3+length(doubleAlphaVec)*(i+1)-1 Nstat+1]);
    errorbar(doubleAlphaVec,meanCk,STDCk,'color',M(K),'Marker','d','MarkerFaceColor','none','MarkerSize',10);
    hold on;
end

% Plot labels and legend
set(0,'defaulttextinterpreter','latex');
xlabel('$2R/R_{Q}$', 'FontSize', 20);

legendCell=strcat(' $K=',strtrim(cellstr(num2str(Kvec(:)))),'$ ');
[h, ~,plots]=legend(legendCell,'Location','northwest', 'interpreter','latex');

if fig ~=0
    title(strcat('$\lambda=$',num2str(lambda),'$, L_x=$',num2str(Lx),'$, \beta=$',num2str(Ec_times_beta),'$/E_c$'));
    ylabel('$\tilde{C}_{k=0}^{-1}$','FontSize', 20);
end
set(h,'FontSize',20);
