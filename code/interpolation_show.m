clear;clc;
% close all;
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
Nrare = round(length(dxOk)/2);
dxRare = linspace(0,xmax*dOk,Nrare)';


aimAngles = [30];
aimAmps   = [0];

% % % 
valuesOk = distribution_former(dxOk,f,aimAngles, aimAmps);
[rstOk,theta] = aec_simulation(valuesOk, dxOk,f);
rstOkT = pow2db(abs(rstOk));rstOkT = rstOkT - max(rstOkT);
% % % 
valuesRare = distribution_former(dxRare,f,aimAngles, aimAmps);
[rstRare,theta] = aec_simulation(valuesRare, dxRare,f);
rstRareT = pow2db(abs(rstRare));rstRareT = rstRareT - max(rstRareT);


figure; hold on
plot(theta, rstOkT);
plot(theta, rstRareT)

lgdt = {
    'ok',
    'rare'
};
legend(lgdt)

%%
figure
hold on
grid on
plcs = zeros(size(dxOk));
plot(dxOk, plcs, 'b*')
plcs = zeros(size(dxRare))+0.1;
plot(dxRare, plcs, 'r*')

ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove x-axis

axis([-inf inf -1.8 1.8])