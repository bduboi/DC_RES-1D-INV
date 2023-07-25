%Import dataset
addpath(genpath('/Users/balthazar/Desktop/Master_Géophysique/M1/Stage_geophysique/Donnees'))
M=  readtable('VES_data.csv');
a = table2array(M(:,1));
m = table2array(M(:,2));
I= table2array(M(:,3));
R = table2array(M(:,4));
std= table2array(M(:,5));

AM = a-m;
BM=a+m;
K = (2.*pi.*1)./((2./AM) - (2./BM) );
pa = K.*R*10^-3;

std =(std.*pa ./ 100);

%define a starting model 
nl = 15;
rho = 50*ones(nl,1)';
thk = logspace(0.3,1.6,nl-1);
x0 = [rho thk]';
[x, obj_fn, lambda_history] = InversionDC1Dres(x0,nl,pa,a,std);
cumsum_thk = cumsum(x(nl+1:end));


%PLOT DU MODÈLE
figure(3);
plot_model(x(1:nl), thk, '.-');
hold on
plot_model(x0(1:nl), thk, '-');
hold off
legend('Modèle Final','Modèle de départ');
xlim([8 100]);
set(gcf, 'position', [1191 50 550 420]);
%PLOT DE L'ÉVOLUTION DE LA FONCTION OBJECTIF ET DU PARAMÈTRE DE
%RÉGULARISAITON
figure(4);
it_number = 0:(length(obj_fn) - 1);
yyaxis left
semilogy(obj_fn ./ obj_fn(1));
semilogy(it_number, obj_fn, '.-', 'MarkerSize', 12);
axl = gca();
xlabel("Nombre d'itérations");
ylabel('Fonction objectif relative');
%grid();
set(axl,'YGrid', 'on', 'YMinorGrid', 'on')
set(axl,'XGrid', 'on', 'XMinorGrid', 'on')

xlim([0 length(obj_fn(obj_fn>0))]);
yyaxis right
semilogy(it_number, lambda_history, '.-', 'MarkerSize', 12);
% grid();
ylabel('Paramètre de régularisation');
set(gcf, 'position', [1191 50 550 420]);
