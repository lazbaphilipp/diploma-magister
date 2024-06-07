clear;clc;
close all;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)
% 
c = physconst("Lightspeed");
% freq params
f = 3e9; % 3 GHz
w = 2*pi*f;
lam = freq2wavelen(f);
k = 2*pi/lam;
% array params
xmin = 0;xmax = 29;
dOk = lam/2;
% NOk = round((xmax-xmin)/dOk);
dxOk = dOk*(xmin:1:xmax)';


% dxOk = dxOk - max(dxOk)/2;
% aimAngles = [0];
% aimAmps   = [0];
aimAngles = [0 -20 30 50];
aimAmps   = [0 -1 -4 -2];

% % % % 

tmpar = [0.7 1.5 2.4 3.3 4.2 5.1 6.1 7.2 8.2 9.2 11.4 12.6 13.8 15]';
dxRareTmp = dOk*([-flipud(tmpar);0;tmpar]);
valuesTmp = distribution_former(dxRareTmp,f,aimAngles, aimAmps);
[rstTmp,theta] = aec_simulation(valuesTmp, dxRareTmp,f);
rstTmpT = mag2db(abs(rstTmp));rstTmpT = rstTmpT - max(rstTmpT);



% % % 
valuesOk = distribution_former(dxOk,f,aimAngles, aimAmps);
% valuesRare = distribution_former(dxOk,f,aimAngles, aimAmps);

[rstOk,theta] = aec_simulation(valuesOk, dxOk,f);
rstOkT = mag2db(abs(rstOk));rstOkT = rstOkT - max(rstOkT);
rstOk1 = rstOk.^2;
rstOk1T = mag2db(abs(rstOk1));rstOk1T = rstOk1T - max(rstOk1T);


dxRare = [-3*dOk; dxOk([1 3 5 9:20 26 28 30]); 32*dOk];
valuesRare = distribution_former(dxRare,f,aimAngles, aimAmps);
% valuesRare = valuesOk([1 3 5 9:20 26 28 30]);
rstRare = aec_simulation(valuesRare, dxRare,f);
rstRareT = mag2db(abs(rstRare));rstRareT = rstRareT - max(rstRareT);

valuesRare1 = valuesRare.*([[1,1, 3, 5, 6, 6.5, 7, 7.4, 7.7,8] fliplr([1,1, 3, 5, 6, 6.5, 7, 7.4, 7.7,8])]');
rstRare1 = aec_simulation(valuesRare1, dxRare,f);
rstRare1T = mag2db(abs(rstRare1));rstRare1T = rstRare1T - max(rstRare1T);


dxRare2 = dxOk([9:20]);
valuesRare2 = valuesOk([9:20]);
rstRare2 = aec_simulation(valuesRare2, dxRare2,f);
rstRare2T = mag2db(abs(rstRare2));rstRare2T = rstRare2T - max(rstRare2T);

valuesRare21 = valuesRare2.*([[5, 5, 5.5, 6, 6, 6.5] fliplr([5, 5, 5.5, 6, 6, 6.5])]');
rstRare21 = aec_simulation(valuesRare21, dxRare2,f);
rstRare21T = mag2db(abs(rstRare21));rstRare21T = rstRare21T - max(rstRare21T);

rstRare3 = rstRare21.*rstRare1;
rstRare3T = mag2db(abs(rstRare3));rstRare3T = rstRare3T - max(rstRare3T);

rstRare1_2=rstRare21.*rstRare;
rstRare1_2T = mag2db(abs(rstRare1_2));rstRare1_2T = rstRare1_2T - max(rstRare1_2T);
%%
% figure
% plot(theta, rstOkT)
% axis([-90 90 -40 1])
% title("ok")
figure; hold on
plot(theta, rstOkT,'b');
plot(theta, rstOk1T);
plot(theta, rstRareT)
plot(theta, rstRare1T)
plot(theta, rstRare2T)
plot(theta, rstRare21T)
plot(theta, rstRare3T)
plot(theta, rstRare1_2T,'r')
plot(theta, rstTmpT,'k')

lgdt = {
    'ok',
    'ok1',
    'rare',
    'rare1',
    'rare2',
    'rare21',
    'rare3',
    'rare1_2',
    'tmp1'
};
legend(lgdt)
axis([-90 90 -40 1])
xlabel("угол, град")
ylabel("мощность принятого сигнала, нормализованная, дБ")
% 
% figure
% plot(theta, rstRare2T,theta, rstRare3T,theta, rstRare4T)
% axis([-90 90 -40 1])
% title("rares")

%%
fig=figure(2);
fig.Position = [100 0 800 150];
hold on
grid on
plcs = zeros(size(dxOk));
plot(dxOk, plcs, 'r*')
plcs = zeros(size(dxRare))+0.1;
plot(dxRare, plcs, 'b*')
plcs = zeros(size(dxRare2))-0.1;
plot(dxRare2, plcs, 'k*')
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove x-axis

axis([-inf inf -0.15 0.15])