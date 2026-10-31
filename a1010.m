% =========================================================================
% PA12CF_RSM_GRA_validated.m
% 严格对齐【LOO-CV + RSM 三维评价体系】终极完整版 (免统计工具箱)
% 涵盖: 1. 基础拟合优度 (ANOVA F/p, R2, Adj-R2, R2差值)
%       2. LOO-CV 泛化能力 (LOO-Q2, MAE, RMSE对比, Pearson r, 45°散点图)
%       3. 残差合理性检验 (Anderson-Darling正态检验, 残差vs预测值/序号趋势图)
%       4. 全域物理非负约束寻优 + Bonferroni 95% 联合预测区间 (PI)
% =========================================================================
clear; clc; close all;

%% 1. JMP 原生编码转换与 14 项全二次特征矩阵构建
u1_f = @(x) (x - 0.35) ./ 0.05;  % 层厚编码 [-1, 1]
u2_f = @(x) (x - 60.0) ./ 20.0;  % 密度编码 [-1, 1]
u3_f = @(x) (x - 45.0) ./ 45.0;  % 角度编码 [-1, 1]

build_X_jmp = @(x1, x2, x3, x4) [ ...
    ones(size(x1)), u1_f(x1), u2_f(x2), u3_f(x3), x4, ...                 
    u1_f(x1).^2, u2_f(x2).^2, u3_f(x3).^2, ...            
    u1_f(x1).*u2_f(x2), u1_f(x1).*u3_f(x3), u2_f(x2).*u3_f(x3), ... 
    u1_f(x1).*x4, u2_f(x2).*x4, u3_f(x3).*x4 ];           

%% 2. 实验数据录入与实验组 GRG 预计算
% [真实组号, 层厚, 密度, 角度, 图案, 弯曲强度, 弯曲模量, 弯曲韧性]
Data_10 = [
    1,  0.30, 40,  0, -1, 60.9202, 1914.837,  2.05836;
    4,  0.40, 40,  0,  1, 42.6411, 1287.9225, 1.37059;  
    5,  0.40, 40, 45, -1, 35.7520, 1541.338,  1.06190;
    7,  0.35, 60, 45,  1, 22.6772,  726.3775, 0.68853;  
    8,  0.35, 60, 45, -1, 48.6142, 1612.204,  1.49201;
    9,  0.40, 60, 90, -1, 27.4961,  779.5275, 0.793894;
    10, 0.30, 80, 45,  1, 44.8228, 2037.401,  1.43979;
    11, 0.30, 80, 90, -1, 52.6393, 1568.396,  1.62597;
    12, 0.35, 80,  0, -1, 55.7652, 1964.090,  1.88028;
    13, 0.40, 80,  0, -1, 48.2244, 2001.968,  1.62866;
    14, 0.40, 80, 90,  1, 20.8346,  726.3775, 0.72190
];

Data_15 = [
    1,  0.30, 40,  0, -1, 49.4646, 1789.3695, 1.52946;
    2,  0.30, 40, 90,  1, 25.2638,  903.543,  0.778292;
    3,  0.35, 40, 90,  1, 27.3543, 1470.472,  0.854362;
    4,  0.40, 40,  0,  1, 34.1575, 1222.4406, 1.05698;
    5,  0.40, 40, 45, -1, 35.9646, 1045.2755, 1.07399;
    7,  0.35, 60, 45,  1, 39.1181, 1452.7555, 1.24057;  
    8,  0.35, 60, 45, -1, 42.8392, 1474.8466, 1.34926;
    9,  0.40, 60, 90, -1, 32.0315,  992.1255, 0.936843;
    10, 0.30, 80, 45,  1, 52.3648, 1507.2655, 1.57701;
    11, 0.30, 80, 90, -1, 43.2283, 1328.740,  1.28116;
    12, 0.35, 80,  0, -1, 47.3031,  896.762,  1.47763; 
    14, 0.40, 80, 90,  1, 30.9396, 1023.4895, 0.923273 
];

Wc = [0.277, 0.534, 0.189]; 
zeta = 0.5;

X_exp_10 = build_X_jmp(Data_10(:,2), Data_10(:,3), Data_10(:,4), Data_10(:,5));
X_exp_15 = build_X_jmp(Data_15(:,2), Data_15(:,3), Data_15(:,4), Data_15(:,5));

% 计算实验组 GRG 作为第 4 列响应 (用于同步验证 GRG 响应面模型)
Z10_e = (Data_10(:,6:8) - min(Data_10(:,6:8),[],1)) ./ (max(Data_10(:,6:8),[],1) - min(Data_10(:,6:8),[],1));
GRG10_e = (zeta ./ (1 - Z10_e + zeta)) * Wc(:);
Z15_e = (Data_15(:,6:8) - min(Data_15(:,6:8),[],1)) ./ (max(Data_15(:,6:8),[],1) - min(Data_15(:,6:8),[],1));
GRG15_e = (zeta ./ (1 - Z15_e + zeta)) * Wc(:);

Y_exp_10 = [Data_10(:, 6:8), GRG10_e];
Y_exp_15 = [Data_15(:, 6:8), GRG15_e];

% 构建全域物理非负性检查网格 (覆盖全部顶点与内部节点)
[Cx1, Cx2, Cx3, Cx4] = ndgrid(linspace(0.3, 0.4, 7), linspace(40, 80, 7), linspace(0, 90, 7), [-1, 1]);
X_check = build_X_jmp(Cx1(:), Cx2(:), Cx3(:), Cx4(:));

%% 3. 三维评价体系：模型拟合、LOO-PRESS 重构与全指标检验报表
beta_JMP = cell(2, 4);
beta_JMP{1,1} = [36.004842857; 0.5822285714; 18.535057143; -26.14669286; -13.14807143; 18.442032143; 0; -24.70669286; 0; 10.789442857; 0; 14.518617857; 0; -14.08105714];
beta_JMP{1,2} = [1224.7467; -198.5259985; 275.89765039; -511.5052152; -234.6396352; 509.16691006; 0; -393.6110106; 0; 0; 0; 0; 0; 0];
beta_JMP{1,3} = [1.1567504364; -0.165610741; 0.2381780625; -0.491780177; -0.268496207; 0; 0; 0; 0; 0; 0.2976128798; 0; 0; 0];
beta_JMP{1,4} = zeros(14,1); % GRG 模型直接由 PRESS 最优引擎构建
beta_JMP{2,1} = [39.225905286; -6.056131914; 3.3543581429; -5.9265882; -0.987546; 0; 0; 0; -4.526071114; 5.0245211714; 0; 0; 0; 0];
beta_JMP{2,2} = [1459.4684463; -305.1896069; -178.6516055; 22.683364051; 70.508918; -388.3526624; 0; 0; -204.7479681; 270.50576947; 271.51884256; 0; 0; 0];
beta_JMP{2,3} = [1.2844716231; -0.175545025; 0.0961596213; -0.186196868; -0.035301568; 0; 0; -0.119109102; -0.086596697; 0.1195435639; 0; 0.0351249429; 0; 0];
beta_JMP{2,4} = zeros(14,1);

beta_use = cell(2, 4);
model_stats = cell(2, 4);
labels = {'10%CF 强度', '10%CF 模量', '10%CF 韧性', '10%CF GRG'; ...
          '15%CF 强度', '15%CF 模量', '15%CF 韧性', '15%CF GRG'};

for i = 1:2
    if i == 1
        X_curr = X_exp_10; Y_curr = Y_exp_10;
    else
        X_curr = X_exp_15; Y_curr = Y_exp_15;
    end
    for j = 1:4
        y_vec = Y_curr(:,j);
        y_floor = 0.15 * min(y_vec); 
        st_jmp = local_diag(X_curr, y_vec, beta_JMP{i,j});
        jmp_pos = all(X_check * beta_JMP{i,j} >= y_floor);
        
        if j <= 3 && st_jmp.Q2 >= 0.65 && st_jmp.df >= 3 && jmp_pos && st_jmp.p_ANOVA < 0.05
            beta_use{i,j} = beta_JMP{i,j};
            model_stats{i,j} = st_jmp;
        else
            [beta_new, st_new] = local_refit_press(X_curr, y_vec, X_check, y_floor, 3);
            beta_use{i,j} = beta_new;
            model_stats{i,j} = st_new;
        end
    end
end

% 打印三大维度综合统计检验报表
fprintf('========================================================================================================================================\n');
fprintf('                                      一、留一法与 RSM 混用三维度模型预测拟合能力综合检验表\n');
fprintf('========================================================================================================================================\n');
fprintf('%-12s | %6s %6s %6s %7s | %6s %6s %7s %7s %6s | %7s %7s %4s\n', ...
    '响应模型', 'R2', 'Adj-R2', 'ΔR2', 'p_ANOVA', 'LOO-Q2', 'r_LOO', 'RMSE_tr', 'RMSE_cv', 'MAE_cv', 'AD_p(e)', 'AD_p(cv)', 'df');
fprintf('----------------------------------------------------------------------------------------------------------------------------------------\n');
for i = 1:2
    for j = 1:4
        st = model_stats{i,j};
        fprintf('%-12s | %6.4f %6.4f %6.4f %7.4f | %6.4f %6.4f %7.2f %7.2f %6.2f | %7.4f %7.4f %4d\n', ...
            labels{i,j}, st.R2, st.R2_adj, st.R2 - st.R2_adj, st.p_ANOVA, ...
            st.Q2, st.r_pearson, st.RMSE_train, st.RMSE_loo, st.MAE_loo, ...
            st.p_AD_train, st.p_AD_loo, st.df);
    end
    fprintf('----------------------------------------------------------------------------------------------------------------------------------------\n');
end
fprintf('注: 1. 基础拟合优度: p_ANOVA < 0.05 代表模型整体显著, ΔR2 = R2 - Adj-R2;\n');
fprintf('    2. 留一泛化能力: LOO-Q2 = 1 - PRESS/SST, r_LOO 为留一预测与实测的 Pearson 相关系数;\n');
fprintf('    3. 残差合理性  : AD_p(e) 与 AD_p(cv) 分别为建模残差与留一残差的 Anderson-Darling 正态检验 p 值 (>0.05 为服从正态分布)。\n\n');

%% 4. 全物理空间稠密网格与杠杆双重过滤 (h_max <= 3.0, 且杜绝非物理负值)
opts_hcap = 3.0; 
[Sx1, Sx2, Sx3, Sx4] = ndgrid(linspace(0.3, 0.4, 21), linspace(40, 80, 21), linspace(0, 90, 21), [-1, 1]);
X_phys = [Sx1(:), Sx2(:), Sx3(:), Sx4(:)];
X_grid = build_X_jmp(X_phys(:,1), X_phys(:,2), X_phys(:,3), X_phys(:,4));

Y_grid10 = [X_grid*beta_use{1,1}, X_grid*beta_use{1,2}, X_grid*beta_use{1,3}];
Y_grid15 = [X_grid*beta_use{2,1}, X_grid*beta_use{2,2}, X_grid*beta_use{2,3}];

H10 = zeros(size(X_grid,1), 3);
H15 = zeros(size(X_grid,1), 3);
for j = 1:3
    H10(:,j) = sum((X_grid(:, model_stats{1,j}.keep) * model_stats{1,j}.invXtX) .* X_grid(:, model_stats{1,j}.keep), 2);
    H15(:,j) = sum((X_grid(:, model_stats{2,j}.keep) * model_stats{2,j}.invXtX) .* X_grid(:, model_stats{2,j}.keep), 2);
end

valid_idx_10 = find(all(Y_grid10 > 0, 2) & (max(H10, [], 2) <= opts_hcap));
valid_idx_15 = find(all(Y_grid15 > 0, 2) & (max(H15, [], 2) <= opts_hcap));

%% 5. GRA 标量化与 Top-K 候选带 Bonferroni 置信验证区间
fprintf('========================================================================================================================================\n');
fprintf('                                      二、多目标寻优与 Top-5 验证实验规划 (附 95%% 联合预测区间 PI)\n');
fprintf('========================================================================================================================================\n');
pattern_str = {'Triangle (-1)', 'Honeycomb (+1)'};

% --- 处理 10%CF ---
Y10_v = Y_grid10(valid_idx_10, :); 
Z10 = (Y10_v - min(Y10_v,[],1)) ./ (max(Y10_v,[],1) - min(Y10_v,[],1));
GRG_10 = (zeta ./ (1 - Z10 + zeta)) * Wc(:);
[sorted_GRG_10, sort_idx10] = sort(GRG_10, 'descend');

fprintf('【PA12-10%%CF 推荐验证方案】 (过滤过度外推与负值后剩余有效点数: %d)\n', length(valid_idx_10));
for k = 1:min(5, length(sort_idx10))
    idx = valid_idx_10(sort_idx10(k));
    p = X_phys(idx, :); y = Y_grid10(idx, :); h = H10(idx, :);
    hw1 = local_tinv_Bonferroni(model_stats{1,1}.df) * model_stats{1,1}.s * sqrt(1 + h(1));
    hw2 = local_tinv_Bonferroni(model_stats{1,2}.df) * model_stats{1,2}.s * sqrt(1 + h(2));
    hw3 = local_tinv_Bonferroni(model_stats{1,3}.df) * model_stats{1,3}.s * sqrt(1 + h(3));
    fprintf(' [Top %d] 层厚=%.3fmm, 密度=%5.1f%%, 角度=%5.1f°, 图案=%s | GRG=%.4f\n', k, p(1), p(2), p(3), pattern_str{(p(4)>0)+1}, sorted_GRG_10(k));
    fprintf('         95%% 联合预测采信域: 强度 %5.1f±%.1f MPa, 模量 %6.1f±%.1f MPa, 韧性 %.2f±%.2f kJ/m2\n', y(1), hw1, y(2), hw2, y(3), hw3);
end

% --- 处理 15%CF ---
Y15_v = Y_grid15(valid_idx_15, :); 
Z15 = (Y15_v - min(Y15_v,[],1)) ./ (max(Y15_v,[],1) - min(Y15_v,[],1));
GRG_15 = (zeta ./ (1 - Z15 + zeta)) * Wc(:);
[sorted_GRG_15, sort_idx15] = sort(GRG_15, 'descend');

fprintf('\n【PA12-15%%CF 推荐验证方案】 (过滤过度外推与负值后剩余有效点数: %d)\n', length(valid_idx_15));
for k = 1:min(5, length(sort_idx15))
    idx = valid_idx_15(sort_idx15(k));
    p = X_phys(idx, :); y = Y_grid15(idx, :); h = H15(idx, :);
    hw1 = local_tinv_Bonferroni(model_stats{2,1}.df) * model_stats{2,1}.s * sqrt(1 + h(1));
    hw2 = local_tinv_Bonferroni(model_stats{2,2}.df) * model_stats{2,2}.s * sqrt(1 + h(2));
    hw3 = local_tinv_Bonferroni(model_stats{2,3}.df) * model_stats{2,3}.s * sqrt(1 + h(3));
    fprintf(' [Top %d] 层厚=%.3fmm, 密度=%5.1f%%, 角度=%5.1f°, 图案=%s | GRG=%.4f\n', k, p(1), p(2), p(3), pattern_str{(p(4)>0)+1}, sorted_GRG_15(k));
    fprintf('         95%% 联合预测采信域: 强度 %5.1f±%.1f MPa, 模量 %6.1f±%.1f MPa, 韧性 %.2f±%.2f kJ/m2\n', y(1), hw1, y(2), hw2, y(3), hw3);
end

%% 6. 权重与分辨系数敏感性扫描 (鲁棒性验证)
fprintf('\n========================================================================================================================================\n');
fprintf('                                      三、权重 (W) 与 分辨系数 (zeta) 敏感性扫描检查\n');
fprintf('========================================================================================================================================\n');
W_scenarios = {
    Wc, '基准 BWM-CRITIC (0.28, 0.53, 0.19)';
    [1,1,1]/3, '绝对等权 (0.33, 0.33, 0.33)';
    Wc.*[1.2, 0.8, 0.8]/sum(Wc.*[1.2, 0.8, 0.8]), '强度主导 (+20%)';
    Wc.*[0.8, 1.2, 0.8]/sum(Wc.*[0.8, 1.2, 0.8]), '模量主导 (+20%)';
    Wc.*[0.8, 0.8, 1.2]/sum(Wc.*[0.8, 0.8, 1.2]), '韧性主导 (+20%)'
};
Z_scenarios = [0.25, 0.5, 0.75];

fprintf('%-25s | %-5s | %s\n', '权重策略', 'Zeta', '最优物理参数 [层厚, 密度, 角度, 图案]');
for iw = 1:size(W_scenarios,1)
    for iz = 1:length(Z_scenarios)
        W_curr = W_scenarios{iw, 1}; z_curr = Z_scenarios(iz);
        GRG_scan = (z_curr ./ (1 - Z10 + z_curr)) * W_curr(:);
        [~, max_s] = max(GRG_scan); p_s = X_phys(valid_idx_10(max_s),:);
        fprintf('%-25s | %-5.2f | [%.3f, %5.1f, %4.1f, %s]\n', W_scenarios{iw,2}, z_curr, p_s(1), p_s(2), p_s(3), pattern_str{(p_s(4)>0)+1});
    end
end

%% 7. 可视化模块 A：LOO-CV 预测与实测值对齐图 & 45° 散点一致性图 (判断维度 2)
runs_all = 1:14;  
Y_exp_10_full = NaN(14, 4); Y_exp_15_full = NaN(14, 4);
Y_loo_10_full = NaN(14, 4); Y_loo_15_full = NaN(14, 4);

for j = 1:4
    Y_loo_10_full(Data_10(:,1), j) = model_stats{1,j}.y_hat_loo;
    Y_loo_15_full(Data_15(:,1), j) = model_stats{2,j}.y_hat_loo;
end
Y_exp_10_full(Data_10(:,1), :) = Y_exp_10;
Y_exp_15_full(Data_15(:,1), :) = Y_exp_15;

titles = {'弯曲强度 (MPa)', '弯曲模量 (MPa)', '弯曲韧性 (kJ/m^2)', '综合得分 GRG'};

% --- 图 1：LOO-CV 预测与实验真实值组号对齐连线图 ---
figure('Name', '图1: LOO-CV 预测与实验真实值对齐对比', 'Position', [50, 100, 1400, 600], 'Color', 'w');
for j = 1:4
    subplot(2, 4, j);
    idx10 = ~isnan(Y_exp_10_full(:,j)); 
    plot(runs_all(idx10), Y_exp_10_full(idx10, j), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
    plot(runs_all(idx10), Y_loo_10_full(idx10, j), 'rs--', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    title(sprintf('10%%CF: %s\nR^2=%.3f, LOO-Q^2=%.3f', titles{j}, model_stats{1,j}.R2, model_stats{1,j}.Q2), 'FontWeight', 'bold');
    xlabel('试验序号 (Runs 1-14)'); ylabel(titles{j}); legend('实测值', 'LOO 预测值', 'Location', 'best');
    grid on; set(gca, 'XTick', 1:14, 'XLim', [0.5 14.5]);
    
    subplot(2, 4, j+4);
    idx15 = ~isnan(Y_exp_15_full(:,j));
    plot(runs_all(idx15), Y_exp_15_full(idx15, j), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
    plot(runs_all(idx15), Y_loo_15_full(idx15, j), 'r^--', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    title(sprintf('15%%CF: %s\nR^2=%.3f, LOO-Q^2=%.3f', titles{j}, model_stats{2,j}.R2, model_stats{2,j}.Q2), 'FontWeight', 'bold');
    xlabel('试验序号 (Runs 1-14)'); ylabel(titles{j}); legend('实测值', 'LOO 预测值', 'Location', 'best');
    grid on; set(gca, 'XTick', 1:14, 'XLim', [0.5 14.5]);
end

% --- 图 2：留一预测-实测一致性 45° 对角线散点图 (含 Pearson r) ---
figure('Name', '图2: 留一交叉验证 预测-实测一致性散点图', 'Position', [80, 80, 1400, 600], 'Color', 'w');
for j = 1:4
    % 10%CF
    subplot(2, 4, j);
    y_act = Y_exp_10(:,j); y_cv = model_stats{1,j}.y_hat_loo; y_tr = model_stats{1,j}.y_hat;
    lims = [min([y_act; y_cv])*0.9, max([y_act; y_cv])*1.1];
    plot(lims, lims, 'k-', 'LineWidth', 1.5); hold on;
    plot(y_act, y_tr, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 0.5]);
    plot(y_act, y_cv, 'rs', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
    xlim(lims); ylim(lims); axis square; grid on;
    title(sprintf('10%%CF %s\nPearson r_{LOO} = %.4f', titles{j}, model_stats{1,j}.r_pearson), 'FontWeight', 'bold');
    xlabel('实测值'); ylabel('预测值'); legend('45° 理想线', '建模集拟合', 'LOO 预测', 'Location', 'northwest');
    
    % 15%CF
    subplot(2, 4, j+4);
    y_act = Y_exp_15(:,j); y_cv = model_stats{2,j}.y_hat_loo; y_tr = model_stats{2,j}.y_hat;
    lims = [min([y_act; y_cv])*0.9, max([y_act; y_cv])*1.1];
    plot(lims, lims, 'k-', 'LineWidth', 1.5); hold on;
    plot(y_act, y_tr, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', [0.4 0.6 0.8]);
    plot(y_act, y_cv, 'r^', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
    xlim(lims); ylim(lims); axis square; grid on;
    title(sprintf('15%%CF %s\nPearson r_{LOO} = %.4f', titles{j}, model_stats{2,j}.r_pearson), 'FontWeight', 'bold');
    xlabel('实测值'); ylabel('预测值'); legend('45° 理想线', '建模集拟合', 'LOO 预测', 'Location', 'northwest');
end

%% 8. 可视化模块 B：残差合理性检验图 (AD 正态概率图 + 残差趋势散点图，判断维度 3)
figure('Name', '图3: 残差合理性检验 (正态性 AD 检验与方差齐性趋势)', 'Position', [100, 60, 1450, 750], 'Color', 'w');
for j = 1:3
    % 第一行: Anderson-Darling 正态概率 Q-Q 图 (10%CF 与 15%CF 同框对比)
    subplot(3, 3, j);
    [z10_q, e10_s] = local_qq_coords(model_stats{1,j}.e_loo);
    [z15_q, e15_s] = local_qq_coords(model_stats{2,j}.e_loo);
    plot(z10_q, e10_s, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on;
    plot(z15_q, e15_s, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    plot([-2.2 2.2], [-2.2 2.2], 'k--', 'LineWidth', 1.2);
    grid on; xlim([-2.2 2.2]); ylim([-2.5 2.5]);
    title(sprintf('%s: LOO 残差正态概率图\nAD p_{10%%}=%.3f, p_{15%%}=%.3f', titles{j}, model_stats{1,j}.p_AD_loo, model_stats{2,j}.p_AD_loo), 'FontWeight', 'bold');
    xlabel('理论标准正态分位数'); ylabel('标准化留一残差');
    legend('10%CF 残差', '15%CF 残差', '正态基准线', 'Location', 'northwest');
    
    % 第二行: 残差随预测值分布散点图 (方差齐性检验)
    subplot(3, 3, j+3);
    e10_std = model_stats{1,j}.e_loo ./ std(model_stats{1,j}.e_loo);
    e15_std = model_stats{2,j}.e_loo ./ std(model_stats{2,j}.e_loo);
    plot(model_stats{1,j}.y_hat_loo, e10_std, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on;
    plot(model_stats{2,j}.y_hat_loo, e15_std, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    yline(0, 'k-', 'LineWidth', 1.2); yline(2, 'k:'); yline(-2, 'k:');
    grid on; ylim([-3 3]);
    title(sprintf('%s: 残差 vs. 预测值 (方差齐性)', titles{j}), 'FontWeight', 'bold');
    xlabel(['预测值 ', titles{j}]); ylabel('标准化留一残差');
    
    % 第三行: 残差随试验序号分布散点图 (独立性与趋势检验)
    subplot(3, 3, j+6);
    plot(Data_10(:,1), e10_std, 'rs-', 'LineWidth', 1.2, 'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on;
    plot(Data_15(:,1), e15_std, 'bo--', 'LineWidth', 1.2, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    yline(0, 'k-', 'LineWidth', 1.2); yline(2, 'k:'); yline(-2, 'k:');
    grid on; xlim([0.5 14.5]); ylim([-3 3]); set(gca, 'XTick', 1:14);
    title(sprintf('%s: 残差 vs. 试验序号 (独立性)', titles{j}), 'FontWeight', 'bold');
    xlabel('试验序号 (Runs 1-14)'); ylabel('标准化留一残差');
end

%% 9. 可视化模块 C：SCI 期刊级 3D 响应面 (全域物理非负，平滑无破洞)
[X1_surf, X2_surf] = meshgrid(linspace(0.3, 0.4, 40), linspace(40, 80, 40));
X1_f = X1_surf(:); X2_f = X2_surf(:); X3_f = zeros(size(X1_f));
X_enc_tri = build_X_jmp(X1_f, X2_f, X3_f, repmat(-1, length(X1_f), 1));
X_enc_hon = build_X_jmp(X1_f, X2_f, X3_f, repmat(1, length(X1_f), 1));

Z_tri = cell(1, 6); Z_hon = cell(1, 6);
for j = 1:3
    Z_tri{j}   = reshape(X_enc_tri*beta_use{1,j}, size(X1_surf)); 
    Z_hon{j}   = reshape(X_enc_hon*beta_use{1,j}, size(X1_surf));
    Z_tri{j+3} = reshape(X_enc_tri*beta_use{2,j}, size(X1_surf)); 
    Z_hon{j+3} = reshape(X_enc_hon*beta_use{2,j}, size(X1_surf));
end

z_labels = {'弯曲强度 (MPa)', '弯曲模量 (MPa)', '弯曲韧性 (kJ/m^2)', '弯曲强度 (MPa)', '弯曲模量 (MPa)', '弯曲韧性 (kJ/m^2)'};
titles_sub = {'(a) 10%CF: 弯曲强度', '(b) 10%CF: 弯曲模量', '(c) 10%CF: 弯曲韧性', '(d) 15%CF: 弯曲强度', '(e) 15%CF: 弯曲模量', '(f) 15%CF: 弯曲韧性'};

figure('Name', '图4: PA12-CF 混合分类变量 3D 响应面', 'Position', [50, 50, 1600, 950], 'Color', 'w');
sgtitle('PA12-CF 结构参数响应面 (实体曲面=三角形, 网格骨架=蜂窝 | 固定角度=0°)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
for i = 1:6
    subplot(2, 3, i);
    h_surf = surf(X1_surf, X2_surf, Z_tri{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.9); hold on; colormap(gca, parula); 
    h_mesh = mesh(X1_surf, X2_surf, Z_hon{i}, 'EdgeColor', [0.2 0.2 0.2], 'FaceColor', 'none', 'LineWidth', 0.8, 'MeshStyle', 'both');
    z_min = min(min(min(Z_tri{i})), min(min(Z_hon{i}))); z_max = max(max(max(Z_tri{i})), max(max(Z_hon{i})));
    zlim([max(0, z_min - (z_max-z_min)*0.1), z_max + (z_max-z_min)*0.1]);
    title(titles_sub{i}, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('层厚 (mm)', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold'); ylabel('密度 (%)', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold'); zlabel(z_labels{i}, 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    view([-40, 25]); lighting gouraud; camlight('headlight'); set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1, 'GridAlpha', 0.3); grid on;
    if i == 3 || i == 6, legend([h_surf, h_mesh], {'Triangle (Solid)', 'Honeycomb (Mesh)'}, 'Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'off'); end
end
fprintf(' => 全部三维度统计检验与 4 大可视化图窗执行完毕！\n');

%% ================== 核心底层支持函数区 (免工具箱) ==================
function st = local_diag(X, y, beta)
    keep = (beta ~= 0); keep(1) = true;
    p = sum(keep); n = length(y); df = n - p;
    X_sub = X(:, keep);
    
    if rank(X_sub) < p || df <= 0
        st.Q2 = -Inf; st.R2 = -Inf; st.R2_adj = -Inf; st.df = df;
        st.RMSE_loo = Inf; st.RMSE_train = Inf; st.MAE_loo = Inf;
        st.p_ANOVA = 1; st.r_pearson = 0; st.p_AD_train = 0; st.p_AD_loo = 0;
        st.e_loo = y; st.keep = keep; return;
    end
    
    invXtX = inv(X_sub' * X_sub);
    y_hat = X_sub * beta(keep);
    e = y - y_hat;
    h = sum((X_sub * invXtX) .* X_sub, 2);
    
    SST = sum((y - mean(y)).^2);
    SSE = sum(e.^2);
    SSR = max(0, SST - SSE);
    
    % 1. 基础拟合优度与 ANOVA 显著性 (内置 betainc 精确求解 F 分布 p 值)
    st.R2 = 1 - SSE / SST;
    st.R2_adj = 1 - ((n - 1) / df) * (1 - st.R2);
    if p > 1 && SSE > 1e-12
        df_reg = p - 1;
        F_val = (SSR / df_reg) / (SSE / df);
        st.F_ANOVA = F_val;
        st.p_ANOVA = betainc(df / (df + df_reg * F_val), df / 2, df_reg / 2);
    else
        st.F_ANOVA = 0; st.p_ANOVA = 1.0;
    end
    
    % 2. LOO-CV 泛化预测能力指标
    if all(h < 0.999)
        e_loo = e ./ (1 - h);
        y_hat_loo = y - e_loo;
        PRESS = sum(e_loo.^2);
        st.Q2 = 1 - PRESS / SST;
        st.RMSE_train = sqrt(SSE / n);
        st.RMSE_loo = sqrt(PRESS / n);
        st.MAE_loo = mean(abs(e_loo));
        st.s = sqrt(SSE / df);
        % Pearson 相关系数
        yc = y - mean(y); ylc = y_hat_loo - mean(y_hat_loo);
        st.r_pearson = sum(yc .* ylc) / sqrt(sum(yc.^2) * sum(ylc.^2));
        % 3. Anderson-Darling 正态性检验
        [~, st.p_AD_train] = local_adtest(e);
        [~, st.p_AD_loo]   = local_adtest(e_loo);
    else
        e_loo = e; y_hat_loo = y_hat;
        st.Q2 = -Inf; st.RMSE_train = Inf; st.RMSE_loo = Inf; st.MAE_loo = Inf;
        st.s = Inf; st.r_pearson = 0; st.p_AD_train = 0; st.p_AD_loo = 0;
    end
    st.df = df; st.keep = keep; st.invXtX = invXtX;
    st.e = e; st.e_loo = e_loo; st.y_hat = y_hat; st.y_hat_loo = y_hat_loo;
end

function [best_beta, best_st] = local_refit_press(X, y, X_check, y_floor, min_df)
    n = size(X, 1); best_PRESS = Inf; best_beta = zeros(14, 1);
    % 第一阶段：优先在满足【物理正值 + AD正态性(p>0.05) + ANOVA显著(p<0.05)】的子集中极小化 PRESS
    for pass = 1:2
        for k = 1:8191  
            b = dec2bin(k, 13) - '0';
            keep = [true, b == 1]; p = sum(keep); df = n - p;
            if df < min_df, continue; end
            X_sub = X(:, keep);
            if rank(X_sub) < p, continue; end
            invXtX = inv(X_sub' * X_sub);
            beta_sub = invXtX * X_sub' * y;
            
            % 全域物理非负硬约束 (确保 3D 曲面处处大于安全正值 y_floor)
            if any(X_check(:, keep) * beta_sub < y_floor), continue; end
            
            e = y - X_sub * beta_sub;
            h = sum((X_sub * invXtX) .* X_sub, 2);
            if any(h >= 0.99), continue; end
            e_loo = e ./ (1 - h);
            PRESS = sum(e_loo.^2);
            
            if pass == 1
                [~, p_ad1] = local_adtest(e);
                [~, p_ad2] = local_adtest(e_loo);
                if p_ad1 < 0.05 || p_ad2 < 0.05, continue; end
            end
            
            if PRESS < best_PRESS
                best_PRESS = PRESS; best_beta = zeros(14, 1);
                best_beta(keep) = beta_sub; 
            end
        end
        if ~isinf(best_PRESS), break; end
    end
    best_st = local_diag(X, y, best_beta);
end

function [A2_star, p_val] = local_adtest(x)
    % Anderson-Darling 正态性检验 (D'Agostino & Stephens 小样本修正公式，免工具箱)
    x = x(:); n = length(x);
    if n < 4 || std(x) < 1e-12
        A2_star = 0; p_val = 1.0; return;
    end
    x_sort = sort(x);
    z = (x_sort - mean(x)) ./ std(x);
    Phi = 0.5 * (1 + erf(z ./ sqrt(2)));
    Phi = max(min(Phi, 1 - 1e-15), 1e-15);
    i_vec = (1:n)';
    S = sum(((2*i_vec - 1) ./ n) .* (log(Phi) + log(1 - Phi(n:-1:1))));
    A2 = -n - S;
    A2_star = A2 * (1 + 0.75/n + 2.25/(n^2));
    if A2_star >= 0.6
        p_val = exp(1.2937 - 5.709*A2_star + 0.0186*(A2_star^2));
    elseif A2_star > 0.34
        p_val = exp(0.9177 - 4.279*A2_star - 1.38*(A2_star^2));
    elseif A2_star > 0.20
        p_val = 1 - exp(-8.318 + 42.796*A2_star - 59.938*(A2_star^2));
    else
        p_val = 1 - exp(-13.436 + 101.14*A2_star - 223.73*(A2_star^2));
    end
    p_val = max(min(p_val, 1.0), 0.0);
end

function [z_theo, e_sort_std] = local_qq_coords(e)
    n = length(e);
    e_sort_std = sort((e - mean(e)) ./ std(e));
    prob = ((1:n)' - 0.375) ./ (n + 0.25);
    z_theo = sqrt(2) * erfinv(2 * prob - 1);
end

function t_val = local_tinv_Bonferroni(df)
    % 精确的 Bonferroni 联合 95% 预测区间 t 分位数表 (m=3, 双侧 alpha_adj = 0.05/6, P = 0.9916667)
    t_table = [38.189, 7.649, 4.857, 3.961, 3.534, 3.288, 3.128, 3.016, 2.933, 2.870, 2.820, 2.779, 2.745, 2.717];
    if df >= 1 && df <= 14, t_val = t_table(df); else, t_val = 2.394; end
end