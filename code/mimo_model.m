clear;clc;
close all;
format compact
format short
c = physconst("Lightspeed");

% sequence parameters
n_psseq = 7;
f_c = 14e9; lam = freq2wavelen(f_c);
Timp = 5e-9;Fimp = 1/Timp;

Fs = 2* Fimp;
Ts = 1/Fs;
R_max = 3e3;T_max = 2* R_max/c;
% % % % % %
T_seq = Timp* 127;
R_blind = T_seq* c/2

% % % % % % /
targets = [
    % Dist
    1.5e3 1e3;
    % Theta, deg
    20 -20;
    % Phi, deg
    25 -25;
    ];

targets_dist = targets(1, :);
targets_time = targets_dist* 2/c;
targets_theta = targets(2, :);
targets_phi = targets(3, :);

aimtable = table(targets_dist', targets_theta', targets_phi', 'VariableNames', {'R, m', 'theta, deg', 'phi, deg'});


d=lam* 0.49;
Nx=16;Ny=16;
tx_e = (0:(Nx-1))* d; %tx_e = tx_e-tx_e/2;
rx_e = (0:(Ny-1))* d; %rx_e = rx_e-rx_e/2;
rx_e=transpose(rx_e);

% % % find TX distributions
tx_distributions = zeros(length(tx_e), size(targets, 2));
for i = 1:size(targets, 2)
    phs_i = 2.* pi.* f_c.* tx_e.* sind(targets_theta(i))./ physconst("Lightspeed");
    tx_distributions(:, i) = exp(1j.* phs_i);
end

% generate M_sequence
ps_phases = zeros(length(tx_e), 2^n_psseq-1);
str = dec2bin(primpoly(n_psseq, 'nodisplay'));
myPoly = strlength(str)-find(str=='1');
initials = dec2bin([randperm(126, length(tx_e)-1), 127]);
initials=reshape(str2num(initials(:)), [length(tx_e), 7]);
for i = 1:length(tx_e)
    pnSequence = comm.PNSequence('Polynomial', myPoly, 'InitialConditions', initials(i, :), 'SamplesPerFrame', 2^n_psseq-1);
    ps_phases(i, :)=step(pnSequence)';
end
% upsample signal and cpnvert to
tx_sig_conv_phases = zeros(length(tx_e), round(Timp/Ts) * ((2^n_psseq-1)));
tx_phases = zeros(length(tx_e), round(Timp/Ts)* ((2^n_psseq-1)+round((T_max-T_seq)/Timp)));
for i = 1:length(tx_e)
    tx_sig_conv_phases(i, :) = fliplr(rectpulse(ps_phases(i, :), round(Timp/Ts))* pi);
    tx_phases(i, :) = rectpulse([ps_phases(i, :) zeros(1, round((T_max-T_seq)/Timp))], round(Timp/Ts))* pi;
end
tx_sig_conv = exp(1i.* tx_sig_conv_phases);

rx_coefs = rectpulse([zeros(1, size(ps_phases, 2)) ones(1, round((T_max-T_seq)/Timp))], round(Timp/Ts));
tx_coefs=1-rx_coefs;

t = 0:Ts:Ts* (length(tx_phases)-1); % system time
% find TX phases
tx_matrix = zeros(size(targets, 2), size(tx_phases, 2));
for i = 1:length(targets_theta)
    n_shift = ceil(targets_time(i)/Ts);
    if(n_shift > length(t))
        n_shift = length(t);
    end
    sig_tmp = exp(1i* tx_phases).* exp(1i* 2* pi* (f_c)* targets_time(i)).* tx_coefs;
    sig_tmp = circshift(sig_tmp, n_shift, 2).* rx_coefs.* tx_distributions(:, i);
    tx_matrix(i, :) = sum(sig_tmp, 1);
end

% % % find RX distributions
rx_distributions = zeros(size(targets, 2), length(rx_e));
for i = 1:size(targets, 2)
    phs_i = 2.* pi.* f_c.* tx_e.* sind(targets_phi(i))./ physconst("Lightspeed");
    rx_distributions(i, :) = exp(1j.* phs_i);
end

% % % create RX signal matrix
rx_matrix = zeros(length(tx_e), length(tx_e), size(tx_matrix, 2));
for i = 1:length(targets_theta)
    sig_tmp = tx_matrix(i, :);
    sig_tmp = permute(sig_tmp, [3 1 2]);
    sig_tmp = sig_tmp.* rx_distributions(i, :);
    sig_tmp = repmat(sig_tmp, length(tx_e), 1, 1);
    rx_matrix = rx_matrix + sig_tmp;
end

% % % convolve signal
for i = 1:length(tx_e)
    for j = 1:length(rx_e)
        sig_tmp = permute(rx_matrix(i, j, :), [1, 3, 2]);
        sig_tmp2 = conv(sig_tmp, tx_sig_conv(i, :), "same");
        rx_matrix(i, j, :) = permute(sig_tmp2, [3 2 1]);
    end
end

% % % find angles
rx_matrix = padarray(rx_matrix, [(128-size(rx_matrix, 1)), (128-size(rx_matrix, 2)), 0], 0, 'post');

res = zeros(size(rx_matrix));
for cori = 1:size(rx_matrix, 3)
    tmp = rx_matrix(:, :, cori);
    tmp = fft2(tmp);
    tmp = fftshift(fftshift(tmp, 1), 2);
    res(:, :, cori) = tmp;
end


% % %
resa = abs(res);
% % % calculate angles
sin_theta = linspace(sin(-pi/2), sin(pi/2), size(resa, 1));
theta_v = rad2deg(asin(sin_theta));
phi_v = theta_v;
[TH, PH] = meshgrid(theta_v, phi_v);

rsto = pow2db(resa);
rsto = rsto-max(rsto, [], 'all');
resa(rsto < -6) = 0;
rsto = pow2db(sum(resa, 3));
rsto = rsto-max(rsto, [], 'all');


figure;clf; hold on
contourf(TH, PH, rsto, [-3:1:0])
axis square
ax = gca;
ax.XTick = [-75:25:75];
grid on;
xlabel("\theta, deg")
ylabel("\phi, deg")
