%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTENDED SIMULATION: IRS vs. Relay in Urban Micro-Cell with SCM
%
% Requirements:
%   1. Multiple Ricean K-factors (Small, Average, High)
%   2. Spectral efficiency SE = log2(1 + SNR) as performance metric
%   3. Empirical CDFs of SE vs. K-factor, pathloss exponent, cluster count
%
% Based on:
%   - Bjornson et al., IEEE Wireless Commun. Lett., 2020 (IRS-relaying)
%   - Bjornson & Demir, mimobook, 2024 (SCM channel models)
%   - ETSI TR 125 996 (Spatial Channel Model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

rng(42);

%% ========================================================================
%  SECTION 1: System Parameters
%  ========================================================================

fc_GHz = 3;
fc = fc_GHz * 1e9;
lambda = 3e8 / fc;

B = 10e6;
noiseFiguredB = 10;
sigma2dBm = -174 + 10*log10(B) + noiseFiguredB;
sigma2 = db2pow(sigma2dBm);

antennaGainS = db2pow(5);
antennaGainR = db2pow(5);
antennaGainD = db2pow(0);

alpha_IRS = 1;

% Fixed transmit power for SE computation (20 dBm = 100 mW)
P_tx = db2pow(20);  % mW

% Geometry
d_SR = 80;
dv   = 10;
d1range = 40:2:100;

% IRS elements for SE analysis
N_irs = 200;

% Monte Carlo
nMonteCarlo = 1000;

%% ========================================================================
%  SECTION 2: Parameter Sweeps (Requirements 1 & 3)
%  ========================================================================

% --- REQUIREMENT 1: Multiple Ricean K-factors ---
% --- Sweep 1: K-factor ---
% UPDATED: Wider spread for significant visual separation
%   0 dB  = equal LOS and scatter (near-Rayleigh, worst case)
%   10 dB = moderate LOS dominance (typical 3GPP UMi)
%   25 dB = very strong LOS (open rooftop, rural, mmWave beam)
K_factors_dB = [0, 10, 25];
K_factor_labels = {'Small (K=0 dB)', 'Average (K=10 dB)', 'High (K=25 dB)'};
nK = length(K_factors_dB);

% --- REQUIREMENT 3a: Multiple pathloss exponents ---
% 3GPP UMi-LOS:  PL = 28 + 22*log10(d) -> exponent ~2.2
% 3GPP UMi-NLOS: PL = 22.7 + 36.7*log10(d) -> exponent ~3.67
% We also test intermediate and harsh exponents
%   0.5 = extreme waveguide/tunnel (guided propagation)
%   1.0 = strong waveguide effect (indoor corridor)
%   2.0 = free-space (Friis equation, theoretical baseline)
%   3.0 = suburban / moderate urban
%   4.5 = dense urban (severe shadowing)
PL_exponents = [0.5, 1.0, 2.0, 3.0, 4.5];
PL_labels = {'$\alpha_{\mathrm{PL}} = 0.5$', ...
             '$\alpha_{\mathrm{PL}} = 1.0$', ...
             '$\alpha_{\mathrm{PL}} = 2.0$', ...
             '$\alpha_{\mathrm{PL}} = 3.0$', ...
             '$\alpha_{\mathrm{PL}} = 4.5$'};
nPL = length(PL_exponents);

% --- REQUIREMENT 3b: Multiple cluster counts ---
% UPDATED: Option C — wider spacing for clear curve separation
%   1   = single dominant scatterer (mmWave, tunnel)
%   6   = moderate scattering (3GPP UMi default, sub-6 GHz)
%   20  = rich scattering (dense indoor)
%   50  = very rich scattering (indoor reflective environment)
%   100 = extreme scattering (reverberation-like, near-deterministic)
cluster_counts = [1, 6, 20, 50, 100];
cluster_labels = {'$L = 1$', '$L = 6$', '$L = 20$', '$L = 50$', '$L = 100$'};
nL = length(cluster_counts);

% Fixed SCM parameters (defaults)
M_subpaths = 20;
AS_deg     = 35;
DS_us      = 0.251;

%% ========================================================================
%  SECTION 3: Pathloss Functions
%  ========================================================================

% Generic pathloss with variable exponent (reference distance d0=1m)
% PL(d) = (lambda/(4*pi))^2 * (d/d0)^(-alpha)
% We normalise to match 3GPP at reference distances
gen_pathloss_LOS = @(d, alpha_pl) db2pow(-28 - 20*log10(fc_GHz) - ...
    alpha_pl * 10 * log10(d));
gen_pathloss_NLOS = @(d, alpha_pl) db2pow(-22.7 - 26*log10(fc_GHz) - ...
    alpha_pl * 10 * log10(d));

% Default 3GPP pathloss (alpha_LOS = 2.2, alpha_NLOS = 3.67)
pathloss_LOS  = @(d) db2pow(-28 - 20*log10(fc_GHz) - 22*log10(d));
pathloss_NLOS = @(d) db2pow(-22.7 - 26*log10(fc_GHz) - 36.7*log10(d));

%% ========================================================================
%  SECTION 4: Main Simulation
%  ========================================================================

fprintf('============================================================\n');
fprintf(' EXTENDED SIMULATION: IRS vs. Relay with SCM\n');
fprintf(' Requirements: K-sweep, SE metric, CDFs\n');
fprintf('============================================================\n');
fprintf('Monte Carlo: %d | Sub-paths/cluster: %d\n', nMonteCarlo, M_subpaths);
fprintf('K-factors: [%s] dB\n', num2str(K_factors_dB));
fprintf('PL exponents: [%s]\n', num2str(PL_exponents));
fprintf('Cluster counts: [%s]\n\n', num2str(cluster_counts));

% -------------------------------------------------------------------------
% 4A: SPECTRAL EFFICIENCY vs. DISTANCE for each K-factor
%     (Fixed: L=6, alpha_NLOS=3.67, N=200)
% -------------------------------------------------------------------------
fprintf('--- Sweep 1: K-factor sweep (L=6, alpha=3.67) ---\n');

% Storage: SE samples for CDF [nMonteCarlo x nD1 x nK x 3schemes]
% Schemes: 1=SISO, 2=DF, 3=IRS
nD1 = length(d1range);
SE_vs_K = zeros(nMonteCarlo, nD1, nK, 3);

for ik = 1:nK
    K_dB = K_factors_dB(ik);
    fprintf('  K = %d dB (%s)...\n', K_dB, K_factor_labels{ik});

    for kd = 1:nD1
        d1 = d1range(kd);
        d_SD = sqrt(d1^2 + dv^2);
        d_RD = sqrt((d1 - d_SR)^2 + dv^2);

        for mc = 1:nMonteCarlo
            % SISO: S->D (NLOS, no Ricean)
            h_SD = gen_scm_scalar(d_SD, false, fc, lambda, 6, M_subpaths, ...
                AS_deg, DS_us, 0, pathloss_LOS, pathloss_NLOS, ...
                antennaGainS, antennaGainD);

            % S->R (LOS with variable K)
            h_SR = gen_scm_scalar(d_SR, true, fc, lambda, 6, M_subpaths, ...
                AS_deg, DS_us, K_dB, pathloss_LOS, pathloss_NLOS, ...
                antennaGainS, antennaGainR);

            % R->D (LOS with variable K)
            h_RD = gen_scm_scalar(d_RD, true, fc, lambda, 6, M_subpaths, ...
                AS_deg, DS_us, K_dB, pathloss_LOS, pathloss_NLOS, ...
                antennaGainR, antennaGainD);

            % --- REQUIREMENT 2: Spectral Efficiency = log2(1 + SNR) ---

            % SISO SE
            SNR_siso = P_tx * abs(h_SD)^2 / sigma2;
            SE_vs_K(mc, kd, ik, 1) = log2(1 + SNR_siso);

            % DF Relay SE (half-duplex: factor 1/2)
            SNR_hop1 = P_tx * abs(h_SR)^2 / sigma2;
            SNR_hop2 = P_tx * abs(h_RD)^2 / sigma2;
            SNR_df = min(SNR_hop1, SNR_hop2);
            SE_vs_K(mc, kd, ik, 2) = 0.5 * log2(1 + SNR_df);

            % IRS SE (with N elements, optimal phase)
            h_irs_total = gen_irs_channel(N_irs, h_SD, h_SR, h_RD, ...
                alpha_IRS, lambda, AS_deg);
            SNR_irs = P_tx * abs(h_irs_total)^2 / sigma2;
            SE_vs_K(mc, kd, ik, 3) = log2(1 + SNR_irs);
        end
    end
end
fprintf('  Done.\n\n');

% -------------------------------------------------------------------------
% 4B: SPECTRAL EFFICIENCY CDFs for different pathloss exponents
%     (Fixed: K=10dB, L=6, d1=80m, N=200)
% -------------------------------------------------------------------------
fprintf('--- Sweep 2: Pathloss exponent sweep (K=10dB, L=6, d1=80m) ---\n');

d1_fixed = 80;
d_SD_fixed = sqrt(d1_fixed^2 + dv^2);
d_RD_fixed = sqrt((d1_fixed - d_SR)^2 + dv^2);

SE_vs_PL = zeros(nMonteCarlo, nPL, 3);

for ip = 1:nPL
    alpha_pl = PL_exponents(ip);
    fprintf('  alpha = %.1f (%s)...\n', alpha_pl, PL_labels{ip});

    % Create pathloss functions with this exponent
    pl_los_var  = @(d) gen_pathloss_LOS(d, alpha_pl);
    pl_nlos_var = @(d) gen_pathloss_NLOS(d, alpha_pl);

    for mc = 1:nMonteCarlo
        h_SD = gen_scm_scalar(d_SD_fixed, false, fc, lambda, 6, ...
            M_subpaths, AS_deg, DS_us, 0, pl_los_var, pl_nlos_var, ...
            antennaGainS, antennaGainD);

        h_SR = gen_scm_scalar(d_SR, true, fc, lambda, 6, M_subpaths, ...
            AS_deg, DS_us, 9, pl_los_var, pl_nlos_var, ...
            antennaGainS, antennaGainR);

        h_RD = gen_scm_scalar(d_RD_fixed, true, fc, lambda, 6, ...
            M_subpaths, AS_deg, DS_us, 9, pl_los_var, pl_nlos_var, ...
            antennaGainR, antennaGainD);

        % SISO
        SE_vs_PL(mc, ip, 1) = log2(1 + P_tx * abs(h_SD)^2 / sigma2);

        % DF
        SE_vs_PL(mc, ip, 2) = 0.5 * log2(1 + min(P_tx*abs(h_SR)^2, ...
            P_tx*abs(h_RD)^2) / sigma2);

        % IRS
        h_total = gen_irs_channel(N_irs, h_SD, h_SR, h_RD, ...
            alpha_IRS, lambda, AS_deg);
        SE_vs_PL(mc, ip, 3) = log2(1 + P_tx * abs(h_total)^2 / sigma2);
    end
end
fprintf('  Done.\n\n');

% -------------------------------------------------------------------------
% 4C: SPECTRAL EFFICIENCY CDFs for different cluster counts
%     (Fixed: K=10dB, alpha=3.67, d1=80m, N=200)
% -------------------------------------------------------------------------
fprintf('--- Sweep 3: Cluster count sweep (K=10dB, d1=80m) ---\n');

SE_vs_L = zeros(nMonteCarlo, nL, 3);

for il = 1:nL
    L = cluster_counts(il);
    fprintf('  L = %d (%s)...\n', L, cluster_labels{il});

    for mc = 1:nMonteCarlo
        h_SD = gen_scm_scalar(d_SD_fixed, false, fc, lambda, L, ...
            M_subpaths, AS_deg, DS_us, 0, pathloss_LOS, pathloss_NLOS, ...
            antennaGainS, antennaGainD);

        h_SR = gen_scm_scalar(d_SR, true, fc, lambda, L, M_subpaths, ...
            AS_deg, DS_us, 9, pathloss_LOS, pathloss_NLOS, ...
            antennaGainS, antennaGainR);

        h_RD = gen_scm_scalar(d_RD_fixed, true, fc, lambda, L, ...
            M_subpaths, AS_deg, DS_us, 9, pathloss_LOS, pathloss_NLOS, ...
            antennaGainR, antennaGainD);

        % SISO
        SE_vs_L(mc, il, 1) = log2(1 + P_tx * abs(h_SD)^2 / sigma2);

        % DF
        SE_vs_L(mc, il, 2) = 0.5 * log2(1 + min(P_tx*abs(h_SR)^2, ...
            P_tx*abs(h_RD)^2) / sigma2);

        % IRS
        h_total = gen_irs_channel(N_irs, h_SD, h_SR, h_RD, ...
            alpha_IRS, lambda, AS_deg);
        SE_vs_L(mc, il, 3) = log2(1 + P_tx * abs(h_total)^2 / sigma2);
    end
end
fprintf('  Done.\n\n');

% -------------------------------------------------------------------------
% 4D: Nmin vs. K-factor (extends original analysis)
% -------------------------------------------------------------------------
fprintf('--- Sweep 4: Nmin vs K-factor ---\n');

Nrange_search = 25:25:800;
Nmin_vs_K = zeros(nD1, nK);

for ik = 1:nK
    K_dB = K_factors_dB(ik);
    fprintf('  K = %d dB...\n', K_dB);

    for kd = 1:nD1
        d1 = d1range(kd);
        d_SD = sqrt(d1^2 + dv^2);
        d_RD = sqrt((d1 - d_SR)^2 + dv^2);

        % Collect median DF power
        P_df_samples = zeros(nMonteCarlo, 1);
        P_irs_samples = zeros(nMonteCarlo, length(Nrange_search));

        for mc = 1:nMonteCarlo
            h_SD = gen_scm_scalar(d_SD, false, fc, lambda, 6, ...
                M_subpaths, AS_deg, DS_us, 0, pathloss_LOS, ...
                pathloss_NLOS, antennaGainS, antennaGainD);
            h_SR = gen_scm_scalar(d_SR, true, fc, lambda, 6, ...
                M_subpaths, AS_deg, DS_us, K_dB, pathloss_LOS, ...
                pathloss_NLOS, antennaGainS, antennaGainR);
            h_RD = gen_scm_scalar(d_RD, true, fc, lambda, 6, ...
                M_subpaths, AS_deg, DS_us, K_dB, pathloss_LOS, ...
                pathloss_NLOS, antennaGainR, antennaGainD);

            beta_SD = abs(h_SD)^2;
            beta_SR = abs(h_SR)^2;
            beta_RD = abs(h_RD)^2;

            SINR_DF = 2^(2*4) - 1;
            if beta_SR >= beta_SD
                P_df_samples(mc) = SINR_DF * sigma2 * ...
                    (beta_SR + beta_RD - beta_SD) / (2*beta_RD*beta_SR);
            else
                P_df_samples(mc) = SINR_DF * sigma2 / beta_SD;
            end

            SINR_tgt = 2^4 - 1;
            for in = 1:length(Nrange_search)
                h_total = gen_irs_channel(Nrange_search(in), h_SD, ...
                    h_SR, h_RD, alpha_IRS, lambda, AS_deg);
                P_irs_samples(mc, in) = SINR_tgt * sigma2 / abs(h_total)^2;
            end
        end

        P_df_med = median(P_df_samples);
        P_irs_med = median(P_irs_samples, 1);

        idx = find(P_irs_med < P_df_med, 1, 'first');
        if ~isempty(idx) && idx > 1
            N_lo = Nrange_search(idx-1); N_hi = Nrange_search(idx);
            P_lo = P_irs_med(idx-1); P_hi = P_irs_med(idx);
            frac = (P_df_med - P_lo) / (P_hi - P_lo);
            Nmin_vs_K(kd, ik) = ceil(N_lo + frac*(N_hi - N_lo));
        elseif ~isempty(idx) && idx == 1
            Nmin_vs_K(kd, ik) = Nrange_search(1);
        else
            Nmin_vs_K(kd, ik) = NaN;
        end
    end
end
fprintf('  Done.\n\n');

%% ========================================================================
%  SECTION 5: FIGURES
%  ========================================================================

scheme_names = {'SISO', 'DF Relay', 'IRS (N=200)'};
scheme_lines = {'--', '-.', '-'};
scheme_colors_base = [0 0 0; 0 0.4470 0.7410; 0.8500 0.3250 0.0980];

% =========================================================================
% FIGURE 1: Mean SE vs. Distance for each K-factor (3 subplots)
% =========================================================================
figure('Position', [50 50 1500 450], 'Name', 'Fig1: Mean SE vs Distance');

for ik = 1:nK
    subplot(1, nK, ik);
    hold on; box on; grid on;

    for s = 1:3
        SE_mean = squeeze(mean(SE_vs_K(:, :, ik, s), 1));
        plot(d1range, SE_mean, scheme_lines{s}, ...
            'Color', scheme_colors_base(s,:), 'LineWidth', 2);
    end

    xlabel('Distance $d_1$ [m]', 'Interpreter', 'Latex');
    ylabel('Mean Spectral Efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
    title(K_factor_labels{ik}, 'Interpreter', 'Latex');
    legend(scheme_names, 'Location', 'NorthEast', 'Interpreter', 'Latex');
    set(gca, 'fontsize', 12);
end

% --- Unify vertical-axis limits across (a)-(c) ---
all_ylims = zeros(nK, 2);
for ik = 1:nK
    subplot(1, nK, ik);
    all_ylims(ik, :) = ylim;
end
common_ylim = [min(all_ylims(:,1)), max(all_ylims(:,2))];

for ik = 1:nK
    subplot(1, nK, ik);
    ylim(common_ylim);
end

sgtitle('Mean Spectral Efficiency vs.\ Distance for Different Ricean $K$-Factors', ...
    'Interpreter', 'Latex', 'FontSize', 14);

exportgraphics(gcf, 'figures/fig5_SE_vs_distance_Kfactors.pdf', ...
    'ContentType', 'vector');

% =========================================================================
% FIGURE 2: CDF of SE at d1=80m for each K-factor (3 panels × 3 schemes)
% =========================================================================
figure('Position', [50 550 1500 450], 'Name', 'Fig2: CDF of SE vs K-factor');

% Find index for d1 = 80 m
[~, idx_d80] = min(abs(d1range - 80));

K_colors = [0.8500 0.3250 0.0980;   % Small K - orange
            0.0    0.4470 0.7410;    % Average K - blue
            0.4660 0.6740 0.1880];   % High K - green

for s = 1:3
    subplot(1, 3, s);
    hold on; box on; grid on;

    for ik = 1:nK
        se_data = SE_vs_K(:, idx_d80, ik, s);
        [f, x] = ecdf(se_data);
        plot(x, f, '-', 'Color', K_colors(ik,:), 'LineWidth', 2);
    end

    xlabel('Spectral Efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
    ylabel('CDF', 'Interpreter', 'Latex');
    title(scheme_names{s}, 'Interpreter', 'Latex');
    legend(K_factor_labels, 'Location', 'SouthEast', ...
        'Interpreter', 'Latex', 'FontSize', 9);
    set(gca, 'fontsize', 12);
end

% --- UNIFIED X-AXIS LIMITS FOR CDF ---
all_xlims_fig6 = zeros(3, 2);
for s = 1:3
    subplot(1, 3, s);
    all_xlims_fig6(s, :) = xlim;
end
common_xlim_fig6 = [min(all_xlims_fig6(:,1)), max(all_xlims_fig6(:,2))];
for s = 1:3
    subplot(1, 3, s);
    xlim(common_xlim_fig6);
end

sgtitle('Empirical CDF of Spectral Efficiency at $d_1=80$~m for Different $K$-Factors', ...
    'Interpreter', 'Latex', 'FontSize', 14);

exportgraphics(gcf, 'figures/fig6_CDF_SE_Kfactors.pdf', ...
    'ContentType', 'vector');

% =========================================================================
% FIGURE 3: CDF of SE at d1=80m for each pathloss exponent
% =========================================================================
figure('Position', [50 50 1500 450], 'Name', 'Fig3: CDF of SE vs PL exponent');

% Colors for 5 pathloss values (0.5 to 4.5)
PL_colors = [0.00 0.80 0.60;    % teal/green  (0.5)
            0.00 0.45 0.74;    % blue        (1.0)
            0.93 0.69 0.13;    % yellow      (2.0)
            0.85 0.33 0.10;    % orange      (3.0)
            0.64 0.08 0.18];   % dark red    (4.5)

for s = 1:3
    subplot(1, 3, s);
    hold on; box on; grid on;

    for ip = 1:nPL
        se_data = SE_vs_PL(:, ip, s);
        [f, x] = ecdf(se_data);
        plot(x, f, '-', 'Color', PL_colors(ip,:), 'LineWidth', 2);
    end

    xlabel('Spectral Efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
    ylabel('CDF', 'Interpreter', 'Latex');
    title(scheme_names{s}, 'Interpreter', 'Latex');
    legend(PL_labels, 'Location', 'SouthEast', ...
        'Interpreter', 'Latex', 'FontSize', 9);
    set(gca, 'fontsize', 12);
end

% --- UNIFIED X-AXIS LIMITS FOR CDF ---
all_xlims_fig7 = zeros(3, 2);
for s = 1:3
    subplot(1, 3, s);
    all_xlims_fig7(s, :) = xlim;
end
common_xlim_fig7 = [min(all_xlims_fig7(:,1)), max(all_xlims_fig7(:,2))];
for s = 1:3
    subplot(1, 3, s);
    xlim(common_xlim_fig7);
end

sgtitle('Empirical CDF of Spectral Efficiency at $d_1=80$~m for Different Pathloss Exponents', ...
    'Interpreter', 'Latex', 'FontSize', 14);

exportgraphics(gcf, 'figures/fig7_CDF_SE_pathloss.pdf', ...
    'ContentType', 'vector');

% =========================================================================
% FIGURE 4: CDF of SE at d1=80m for each cluster count
% =========================================================================
figure('Position', [50 550 1500 450], 'Name', 'Fig4: CDF of SE vs Clusters');

% Colors for 5 cluster count values (L=1 to L=100)
L_colors = [0.6350 0.0780 0.1840;   % L=1   dark red
            0.8500 0.3250 0.0980;    % L=6   orange
            0.9290 0.6940 0.1250;    % L=20  yellow
            0.0    0.4470 0.7410;    % L=50  blue
            0.4940 0.1840 0.5560];   % L=100 purple

for s = 1:3
    subplot(1, 3, s);
    hold on; box on; grid on;

    for il = 1:nL
        se_data = SE_vs_L(:, il, s);
        [f, x] = ecdf(se_data);
        plot(x, f, '-', 'Color', L_colors(il,:), 'LineWidth', 2);
    end

    xlabel('Spectral Efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
    ylabel('CDF', 'Interpreter', 'Latex');
    title(scheme_names{s}, 'Interpreter', 'Latex');
    legend(cluster_labels, 'Location', 'SouthEast', ...
        'Interpreter', 'Latex', 'FontSize', 9);
    set(gca, 'fontsize', 12);
end

% --- UNIFIED X-AXIS LIMITS ---
all_xlims_fig8 = zeros(3, 2);
for s = 1:3
    subplot(1, 3, s);
    all_xlims_fig8(s, :) = xlim;
end
common_xlim_fig8 = [min(all_xlims_fig8(:,1)), max(all_xlims_fig8(:,2))];
for s = 1:3
    subplot(1, 3, s);
    xlim(common_xlim_fig8);
end

sgtitle('Empirical CDF of Spectral Efficiency at $d_1=80$~m for Different Cluster Counts', ...
    'Interpreter', 'Latex', 'FontSize', 14);

exportgraphics(gcf, 'figures/fig8_CDF_SE_clusters.pdf', ...
    'ContentType', 'vector');

% =========================================================================
% FIGURE 5: Nmin vs. Distance for each K-factor
% =========================================================================
figure('Position', [100 100 800 500], 'Name', 'Fig5: Nmin vs K-factor');
hold on; box on; grid on;

K_line_styles = {'-o', '-s', '-^'};
for ik = 1:nK
    valid = ~isnan(Nmin_vs_K(:, ik));
    plot(d1range(valid), Nmin_vs_K(valid, ik), K_line_styles{ik}, ...
        'Color', K_colors(ik,:), 'LineWidth', 2, 'MarkerSize', 4);
end

xlabel('Distance $d_1$ [m]', 'Interpreter', 'Latex');
ylabel('$N_{\min}$', 'Interpreter', 'Latex');
title('$N_{\min}$ vs.\ Distance for Different Ricean $K$-Factors (SCM)', ...
    'Interpreter', 'Latex');
legend(K_factor_labels, 'Interpreter', 'Latex', 'Location', 'North');
set(gca, 'fontsize', 14);

exportgraphics(gcf, 'figures/fig9_Nmin_vs_Kfactor.pdf', ...
    'ContentType', 'vector');

% =========================================================================
% FIGURE 6: Combined SE comparison — bar chart at d1=80m
% =========================================================================
figure('Position', [100 650 900 500], 'Name', 'Fig6: SE Bar Chart');

% Compute median SE at d1=80m for each K and scheme
SE_bar_data = zeros(nK, 3);
SE_bar_10pct = zeros(nK, 3);  % 10th percentile (outage)
SE_bar_90pct = zeros(nK, 3);  % 90th percentile

for ik = 1:nK
    for s = 1:3
        se_data = SE_vs_K(:, idx_d80, ik, s);
        SE_bar_data(ik, s) = median(se_data);
        SE_bar_10pct(ik, s) = prctile(se_data, 10);
        SE_bar_90pct(ik, s) = prctile(se_data, 90);
    end
end

b = bar(SE_bar_data, 'grouped');
b(1).FaceColor = scheme_colors_base(1,:);
b(2).FaceColor = scheme_colors_base(2,:);
b(3).FaceColor = scheme_colors_base(3,:);

hold on; box on; grid on;

% Add error bars for 10th-90th percentile range
ngroups = nK;
nbars = 3;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for s = 1:nbars
    x_pos = (1:ngroups) - groupwidth/2 + (2*s-1)*groupwidth/(2*nbars);
    errorbar(x_pos, SE_bar_data(:,s), ...
        SE_bar_data(:,s) - SE_bar_10pct(:,s), ...
        SE_bar_90pct(:,s) - SE_bar_data(:,s), ...
        'k.', 'LineWidth', 1.5);
end

set(gca, 'XTickLabel', {'K=0 dB\newline(Small)', ...
    'K=10 dB\newline(Average)', 'K=25 dB\newline(High)'});
xlabel('Ricean $K$-Factor', 'Interpreter', 'Latex');
ylabel('Spectral Efficiency [bit/s/Hz]', 'Interpreter', 'Latex');
title('Median SE at $d_1=80$~m (error bars: 10th--90th percentile)', ...
    'Interpreter', 'Latex');
legend(scheme_names, 'Location', 'NorthWest', 'Interpreter', 'Latex');
set(gca, 'fontsize', 14);

exportgraphics(gcf, 'figures/fig10_SE_bar_chart.pdf', ...
    'ContentType', 'vector');

% =========================================================================
% FIGURE 7: Heatmap — Median SE gain of IRS over DF vs. (K, d1)
% =========================================================================
figure('Position', [100 100 900 400], 'Name', 'Fig7: SE Gain Heatmap');

SE_gain_IRS_over_DF = zeros(nD1, nK);
for ik = 1:nK
    for kd = 1:nD1
        med_irs = median(SE_vs_K(:, kd, ik, 3));
        med_df  = median(SE_vs_K(:, kd, ik, 2));
        SE_gain_IRS_over_DF(kd, ik) = med_irs - med_df;
    end
end

imagesc(K_factors_dB, d1range, SE_gain_IRS_over_DF);
colorbar;
colormap(jet);
xlabel('Ricean $K$-Factor [dB]', 'Interpreter', 'Latex');
ylabel('Distance $d_1$ [m]', 'Interpreter', 'Latex');
title('Median SE Gain: IRS (N=200) over DF Relay [bit/s/Hz]', ...
    'Interpreter', 'Latex');
set(gca, 'fontsize', 14, 'YDir', 'normal');
set(gca, 'XTick', K_factors_dB);

exportgraphics(gcf, 'figures/fig11_SE_gain_heatmap.pdf', ...
    'ContentType', 'vector');

%% ========================================================================
%  SECTION 6: PRINT SUMMARY TABLES
%  ========================================================================

fprintf('\n============================================================\n');
fprintf(' SUMMARY: Median Spectral Efficiency at d1=80m, N_IRS=200\n');
fprintf('============================================================\n');
fprintf('%-20s | %-12s | %-12s | %-12s\n', ...
    'K-factor', 'SISO', 'DF Relay', 'IRS(N=200)');
fprintf('--------------------+--------------+--------------+--------------\n');
for ik = 1:nK
    fprintf('%-20s | %10.3f   | %10.3f   | %10.3f\n', ...
        K_factor_labels{ik}, ...
        SE_bar_data(ik,1), SE_bar_data(ik,2), SE_bar_data(ik,3));
end

fprintf('\n%-20s | %-12s | %-12s | %-12s\n', ...
    'K-factor', 'SISO(10%%)', 'DF(10%%)', 'IRS(10%%)');
fprintf('--------------------+--------------+--------------+--------------\n');
for ik = 1:nK
    fprintf('%-20s | %10.3f   | %10.3f   | %10.3f\n', ...
        K_factor_labels{ik}, ...
        SE_bar_10pct(ik,1), SE_bar_10pct(ik,2), SE_bar_10pct(ik,3));
end

fprintf('\n============================================================\n');
fprintf(' Nmin at d1=80m for each K-factor\n');
fprintf('============================================================\n');
for ik = 1:nK
    fprintf('  %s: Nmin = %d\n', K_factor_labels{ik}, ...
        Nmin_vs_K(idx_d80, ik));
end

fprintf('\n============================================================\n');
fprintf(' IRS SE Gain over DF at d1=80m\n');
fprintf('============================================================\n');
for ik = 1:nK
    fprintf('  %s: +%.3f bit/s/Hz\n', K_factor_labels{ik}, ...
        SE_gain_IRS_over_DF(idx_d80, ik));
end

fprintf('\nSimulation complete. Figures saved to figures/ directory.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function h = gen_scm_scalar(d, isLOS, fc, lambda, L, M, AS_deg, DS_us, ...
    K_dB, pl_los_fn, pl_nlos_fn, antGainTx, antGainRx)

    offsets_half = [0.0447 0.1413 0.2492 0.3715 0.5129 ...
                    0.6797 0.8844 1.1481 1.5195 2.1551];
    sp_offsets = [-fliplr(offsets_half), offsets_half];

    % If M ~= 20, resample offsets
    if M ~= 20
        sp_offsets = linspace(sp_offsets(1), sp_offsets(end), M);
    end

    % Cluster delays
    tau_raw = -DS_us * log(rand(L,1));
    tau = sort(tau_raw - min(tau_raw));

    % Cluster powers: ETSI Step 6
    cluster_pow = 10.^(-tau / DS_us + 0.2 * randn(L,1));
    cluster_pow = cluster_pow / sum(cluster_pow);

    % Cluster AoD
    cluster_AoD = AS_deg * randn(L,1);

    % Sub-path angles
    c_AS = 2;
    sp_angles = cluster_AoD + c_AS * sp_offsets;

    % Random phases
    phases = 2*pi*rand(L, M);

    % Sum sub-paths
    h = 0;
    for l = 1:L
        for m = 1:M
            h = h + sqrt(cluster_pow(l)/M) * exp(1i*phases(l,m));
        end
    end

    % Ricean LOS
    if isLOS
        K_lin = db2pow(K_dB);
        h = h / sqrt(1 + K_lin) + sqrt(K_lin/(1+K_lin)) * exp(1i*2*pi*rand());
    end

    % Pathloss
    if isLOS
        PL = pl_los_fn(d);
    else
        PL = pl_nlos_fn(d);
    end
    h = h * sqrt(PL * antGainTx * antGainRx);
end


function h_total = gen_irs_channel(N, h_SD, h_SR, h_RD, alpha, lambda, AS_deg)

    d_elem = lambda / 2;
    theta_SI = AS_deg * randn() * pi/180;
    theta_ID = AS_deg * randn() * pi/180;
    dphi_SI = 2*pi * d_elem * sin(theta_SI) / lambda;
    dphi_ID = 2*pi * d_elem * sin(theta_ID) / lambda;

    abs_SR = abs(h_SR);
    abs_RD = abs(h_RD);

    h_irs_sum = N * alpha * abs_SR * abs_RD;
    h_total = abs(h_SD) + h_irs_sum;
end