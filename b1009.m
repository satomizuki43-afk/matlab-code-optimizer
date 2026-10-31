% =========================================================================
% PA12-CF 夹芯结构多目标轻量化优化算法 (纯净防报错版)
% -------------------------------------------------------------------------
% 核心升级: 
% 1. 模型底座: LOO-PRESS 预测极小化 + 强遗传原则 (防过拟合)
% 2. 算法引擎: 修复版 NSGA-II (加入非负约束 Y>0 与 杠杆红线 h<=3.0 双重惩罚)
% 3. 决策机制: 纳什均衡博弈 (Nash Bargaining Solution) 规避非凸边界陷阱
% 4. 可视化  : 3D 曲面截断屏蔽策略 (NaN Masking) 掏空反物理预测负值
% =========================================================================
clear; clc; close all;

%% 1. JMP 原生编码转换与设计矩阵构建
u1_f = @(x) (x - 0.35) ./ 0.05;  % X1: 层厚 (0.3~0.4)
u2_f = @(x) (x - 60.0) ./ 20.0;  % X2: 密度 (40~80)
u3_f = @(x) (x - 45.0) ./ 45.0;  % X3: 角度 (0~90)
build_X_jmp = @(x1, x2, x3, x4) [ ...
    ones(size(x1)), u1_f(x1), u2_f(x2), u3_f(x3), x4, ...
    u1_f(x1).^2, u2_f(x2).^2, u3_f(x3).^2, ...
    u1_f(x1).*u2_f(x2), u1_f(x1).*u3_f(x3), u2_f(x2).*u3_f(x3), ...
    u1_f(x1).*x4, u2_f(x2).*x4, u3_f(x3).*x4 ];

%% 2. 实验数据先行录入
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
X_exp_10 = build_X_jmp(Data_10(:,2), Data_10(:,3), Data_10(:,4), Data_10(:,5));
X_exp_15 = build_X_jmp(Data_15(:,2), Data_15(:,3), Data_15(:,4), Data_15(:,5));
Y_exp_10 = Data_10(:, 6:8); Y_exp_15 = Data_15(:, 6:8);

%% 3. [后验环节] 执行 LOO-PRESS 动态重拟合
beta_JMP = cell(2, 3);
beta_JMP{1,1} = [36.004842857; 0.5822285714; 18.535057143; -26.14669286; -13.14807143; 18.442032143; 0; -24.70669286; 0; 10.789442857; 0; 14.518617857; 0; -14.08105714];
beta_JMP{1,2} = [1224.7467; -198.5259985; 275.89765039; -511.5052152; -234.6396352; 509.16691006; 0; -393.6110106; 0; 0; 0; 0; 0; 0];
beta_JMP{1,3} = [1.1567504364; -0.165610741; 0.2381780625; -0.491780177; -0.268496207; 0; 0; 0; 0; 0; 0.2976128798; 0; 0; 0];
beta_JMP{2,1} = [39.225905286; -6.056131914; 3.3543581429; -5.9265882; -0.987546; 0; 0; 0; -4.526071114; 5.0245211714; 0; 0; 0; 0];
beta_JMP{2,2} = [1459.4684463; -305.1896069; -178.6516055; 22.683364051; 70.508918; -388.3526624; 0; 0; -204.7479681; 270.50576947; 271.51884256; 0; 0; 0];
beta_JMP{2,3} = [1.2844716231; -0.175545025; 0.0961596213; -0.186196868; -0.035301568; 0; 0; -0.119109102; -0.086596697; 0.1195435639; 0; 0.0351249429; 0; 0];

fprintf('====================================================================\n');
fprintf('  [后验启动] 执行严格遵守强遗传原则的交叉验证重修剪...\n');
fprintf('====================================================================\n');
beta_use = cell(2, 3); model_stats = cell(2, 3);
for i = 1:2
    if i == 1
        X_curr = X_exp_10; Y_curr = Y_exp_10; 
    else
        X_curr = X_exp_15; Y_curr = Y_exp_15; 
    end
    for j = 1:3
        st_jmp = local_diag(X_curr, Y_curr(:,j), beta_JMP{i,j});
        if st_jmp.Q2 >= 0.5 && st_jmp.df >= 3 && local_check_heredity(st_jmp.keep, 'strong')
            beta_use{i,j} = beta_JMP{i,j}; model_stats{i,j} = st_jmp;
        else
            [beta_use{i,j}, model_stats{i,j}] = local_refit_press_heredity(X_curr, Y_curr(:,j), 3);
        end
    end
end

%% 4. MCDM 评价机制参数设定
Y_min_global = [17.5295,  726.3775, 0.68853]; 
Y_max_global = [60.9202, 2442.6917, 5.66375];
Wc = [0.277, 0.534, 0.189]; 
zeta = 0.5;

%% 5. 标准 NSGA-II 算法执行: GRG (力学) vs Density (轻量化)
fprintf('\n====================================================================\n');
fprintf('  [寻优进化] 标准 NSGA-II 引擎 (内置 h<=3.0 与 Y>0 惩罚墙) 求解中...\n');
fprintf('====================================================================\n');
rng(42); 
knee_params = zeros(2, 4); knee_preds = zeros(2, 3); 
knee_scores = zeros(2, 2); knee_h = zeros(2, 3);
PF_data = cell(2,1);

for row = 1:2
    X_pf_comb = []; Y_obj_comb = []; Y_pred_comb = [];
    
    for pattern = [-1, 1]
        fit_fcn = @(u) physical_surrogate_fitness(u, pattern, row, beta_use, ...
            model_stats, build_X_jmp, Y_min_global, Y_max_global, Wc, zeta);
            
        [x_nsga, obj_nsga] = nsga2_engine_standard(fit_fcn, 3, [0.30, 40, 0], [0.40, 80, 90], 300, 200);
        
        valid_nsga = obj_nsga(:,2) < 1000; % 剔除被越界惩罚的高分个体
        x_nsga = x_nsga(valid_nsga, :);
        obj_nsga = obj_nsga(valid_nsga, :);
        
        if isempty(x_nsga), continue; end
        
        X_mat = build_X_jmp(x_nsga(:,1), x_nsga(:,2), x_nsga(:,3), repmat(pattern, size(x_nsga,1), 1));
        preds = [X_mat*beta_use{row,1}, X_mat*beta_use{row,2}, X_mat*beta_use{row,3}];
        
        X_pf_comb = [X_pf_comb; [x_nsga, repmat(pattern, size(x_nsga,1), 1)]];
        Y_obj_comb = [Y_obj_comb; obj_nsga];
        Y_pred_comb = [Y_pred_comb; preds];
    end
    
    [fronts, ~] = fast_non_dominated_sort(Y_obj_comb);
    is_pf = false(size(Y_obj_comb, 1), 1);
    if ~isempty(fronts{1}), is_pf(fronts{1}) = true; end
    
    X_pf = X_pf_comb(is_pf, :);
    Y_obj_pf = Y_obj_comb(is_pf, :);
    Y_pred_pf = Y_pred_comb(is_pf, :);
    
    PF_GRG = -Y_obj_pf(:, 1);  
    PF_Density = Y_obj_pf(:, 2);  
    
    % [决策升级]：智能定位 Knee Point (基于纳什均衡博弈 Nash Bargaining)
    U_GRG = (PF_GRG - min(PF_GRG)) ./ (max(PF_GRG) - min(PF_GRG) + 1e-6);
    U_Den = (max(PF_Density) - PF_Density) ./ (max(PF_Density) - min(PF_Density) + 1e-6);
    Nash_Product = U_GRG .* U_Den;
    [~, knee_idx] = max(Nash_Product);
    
    knee_params(row, :) = X_pf(knee_idx, :);
    knee_preds(row, :) = Y_pred_pf(knee_idx, :);
    knee_scores(row, :) = [PF_GRG(knee_idx), PF_Density(knee_idx)];
    PF_data{row} = [PF_Density, PF_GRG, X_pf];
    
    X_k = build_X_jmp(knee_params(row,1), knee_params(row,2), knee_params(row,3), knee_params(row,4));
    for j = 1:3
        knee_h(row, j) = sum((X_k(model_stats{row,j}.keep) * model_stats{row,j}.invXtX) .* X_k(model_stats{row,j}.keep), 2);
    end
end

%% 6. 打印工程推荐报告
pattern_str = {'Triangle (-1)', 'Honeycomb (+1)'};
sys_name = {'PA12-10%CF', 'PA12-15%CF'};
for i = 1:2
    fprintf(' -> 体系: %s (Nash 均衡决策点)\n', sys_name{i});
    p_str = pattern_str{(knee_params(i,4)>0)+1};
    
    hw1 = local_tinv_Bonferroni(model_stats{i,1}.df) * model_stats{i,1}.s * sqrt(1 + knee_h(i,1));
    hw2 = local_tinv_Bonferroni(model_stats{i,2}.df) * model_stats{i,2}.s * sqrt(1 + knee_h(i,2));
    hw3 = local_tinv_Bonferroni(model_stats{i,3}.df) * model_stats{i,3}.s * sqrt(1 + knee_h(i,3));
    
    % 安全换行输出，避免意外的行尾报错
    fprintf(['   [最优工艺]: 层厚=%.3fmm | 密度=%5.1f%% | ', ...
             '角度=%5.1f度 | 图案=%s\n'], ...
             knee_params(i,1), knee_params(i,2), knee_params(i,3), p_str);
             
    fprintf(['   [95%% 预测区间]: 强度 %5.1f +/- %4.1f | ', ...
             '模量 %6.1f +/- %5.1f | 韧性 %.2f +/- %.2f\n'], ...
             knee_preds(i,1), hw1, knee_preds(i,2), hw2, knee_preds(i,3), hw3);
             
    fprintf('--------------------------------------------------------------------\n');
end

%% 7. 高精度 Pareto 曲线可视化
figure('Name', 'Pareto Optimization for Lightweight Design', 'Color', 'w', 'Position', [100, 150, 1100, 500]);
color_10 = [0.12, 0.52, 0.90]; color_15 = [0.90, 0.35, 0.15];
for row = 1:2
    subplot(1, 2, row); hold on; grid on; box on;
    data = PF_data{row};
    [sorted_Den, sort_idx] = sort(data(:,1)); sorted_GRG = data(sort_idx, 2);
    
    plot(sorted_Den, sorted_GRG, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
    c_use = color_10; if row == 2, c_use = color_15; end
    scatter(sorted_Den, sorted_GRG, 48, c_use, 'filled', 'MarkerEdgeColor', 'k');
    
    k_den = knee_scores(row, 2); k_grg = knee_scores(row, 1);
    plot(k_den, k_grg, 'p', 'MarkerSize', 18, 'MarkerFaceColor', [1 0.84 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    
    str_knee = sprintf('Nash Point\nDen: %.1f%%\nGRG: %.3f', knee_params(row,2), k_grg);
    text(k_den+0.5, k_grg, str_knee, 'FontSize', 10, 'FontWeight', 'bold');
    
    xlabel('物理填充密度 (%)', 'FontWeight', 'bold'); 
    ylabel('综合力学 GRG', 'FontWeight', 'bold');
    title(sprintf('%s: 强度 vs 轻量化 博弈', sys_name{row}));
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
end

%% 8. 三维画图：SCI 期刊级 3D 响应面 (包含 NaN 掏空策略)
[X1_surf, X2_surf] = meshgrid(linspace(0.3, 0.4, 40), linspace(40, 80, 40));
X1_f = X1_surf(:); X2_f = X2_surf(:); X3_f = zeros(size(X1_f));
X_enc_tri = build_X_jmp(X1_f, X2_f, X3_f, repmat(-1, length(X1_f), 1));
X_enc_hon = build_X_jmp(X1_f, X2_f, X3_f, repmat(1, length(X1_f), 1));
Z_tri = cell(1, 6); Z_hon = cell(1, 6);
for j = 1:3
    Z_tri{j} = reshape(X_enc_tri*beta_use{1,j}, size(X1_surf)); 
    Z_hon{j} = reshape(X_enc_hon*beta_use{1,j}, size(X1_surf));
    Z_tri{j+3} = reshape(X_enc_tri*beta_use{2,j}, size(X1_surf)); 
    Z_hon{j+3} = reshape(X_enc_hon*beta_use{2,j}, size(X1_surf));
end

% [视觉屏蔽策略]：掏空所有非物理的负值预测
for i = 1:6
    Z_tri{i}(Z_tri{i} <= 0) = NaN;
    Z_hon{i}(Z_hon{i} <= 0) = NaN;
end

z_labels = {'弯曲强度 (MPa)', '弯曲模量 (MPa)', '弯曲韧性 (kJ/m^2)', ...
            '弯曲强度 (MPa)', '弯曲模量 (MPa)', '弯曲韧性 (kJ/m^2)'};
titles_sub = {'(a) 10%CF 强度', '(b) 10%CF 模量', '(c) 10%CF 韧性', ...
              '(d) 15%CF 强度', '(e) 15%CF 模量', '(f) 15%CF 韧性'};
figure('Name', 'PA12-CF 混合分类变量响应面图', 'Position', [50, 50, 1600, 950], 'Color', 'w');
sgtitle('PA12-CF 结构参数响应面 (截断屏蔽策略)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

for i = 1:6
    subplot(2, 3, i);
    h_surf = surf(X1_surf, X2_surf, Z_tri{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.9); hold on; 
    colormap(gca, parula); 
    h_mesh = mesh(X1_surf, X2_surf, Z_hon{i}, 'EdgeColor', [0.2 0.2 0.2], 'FaceColor', 'none', 'LineWidth', 0.8, 'MeshStyle', 'both');
    
    title(titles_sub{i}, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('层厚 (mm)', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold'); 
    ylabel('密度 (%)', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold'); 
    zlabel(z_labels{i}, 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
    
    view([-40, 25]); lighting gouraud; camlight('headlight'); 
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1, 'GridAlpha', 0.3); grid on;
end

%% ========================================================================
%% 附属底层函数库
% =========================================================================

% --- 适应度函数：内置双重惩罚防线 (h<=3.0 & Y>0) ---
function objs = physical_surrogate_fitness(X, pattern, row, beta_use, model_stats, build_X_jmp, Y_min, Y_max, Wc, zeta)
    N = size(X, 1);
    X_mat = build_X_jmp(X(:,1), X(:,2), X(:,3), repmat(pattern, N, 1));
    Y_pred = [X_mat*beta_use{row,1}, X_mat*beta_use{row,2}, X_mat*beta_use{row,3}];
    
    % [绝杀防线 1]：物理非负约束
    neg_violation = sum(Y_pred <= 0, 2); 
    
    Z = (Y_pred - Y_min) ./ (Y_max - Y_min);
    GRG = sum((zeta ./ (1 - Z + zeta)) .* Wc, 2);
    
    % [绝杀防线 2]：DOE 外推红线
    h_max = zeros(N, 1);
    for j = 1:3
        st = model_stats{row,j};
        h_curr = sum((X_mat(:, st.keep) * st.invXtX) .* X_mat(:, st.keep), 2);
        h_max = max(h_max, h_curr);
    end
    h_violation = max(0, h_max - 3.0);
    
    penalty = (h_violation + neg_violation) * 10000; 
    objs = [-GRG + penalty, X(:,2) + penalty];
end

% --- 标准 NSGA-II 引擎 ---
function [X_pf, Y_pf] = nsga2_engine_standard(fit_fcn, nvars, lb, ub, popSize, maxGen)
    pop = repmat(lb, popSize, 1) + rand(popSize, nvars) .* repmat((ub - lb), popSize, 1);
    obj = fit_fcn(pop);
    for gen = 1:maxGen
        [fronts, rank] = fast_non_dominated_sort(obj);
        cd = crowding_distance(obj, fronts);
        parent_idx = zeros(popSize, 1);
        for i = 1:popSize
            t1 = randi(popSize); t2 = randi(popSize);
            if rank(t1) < rank(t2)
                parent_idx(i) = t1;
            elseif rank(t1) > rank(t2)
                parent_idx(i) = t2;
            else
                if cd(t1) >= cd(t2), parent_idx(i) = t1; else, parent_idx(i) = t2; end
            end
        end
        parents = pop(parent_idx, :); 
        offspring = zeros(popSize, nvars);
        for i = 1:2:popSize
            if rand < 0.9 
                u = rand(1, nvars); beta = zeros(1, nvars);
                beta(u <= 0.5) = (2 * u(u <= 0.5)).^(1/21); 
                beta(u > 0.5) = (1 ./ (2 * (1 - u(u > 0.5)))).^(1/21);
                offspring(i, :)   = 0.5 * ((1 + beta).*parents(i,:) + (1 - beta).*parents(i+1,:));
                offspring(i+1, :) = 0.5 * ((1 - beta).*parents(i,:) + (1 + beta).*parents(i+1,:));
            else
                offspring(i, :) = parents(i, :); 
                offspring(i+1, :) = parents(i+1, :); 
            end
        end
        for i = 1:popSize
            for j = 1:nvars
                if rand < (1 / nvars)
                    r = rand; 
                    if r < 0.5, delta = (2 * r)^(1/21) - 1; else, delta = 1 - (2 * (1 - r))^(1/21); end
                    offspring(i, j) = offspring(i, j) + delta * (ub(j) - lb(j));
                end
            end
        end
        offspring = max(min(offspring, repmat(ub, popSize, 1)), repmat(lb, popSize, 1));
        off_obj = fit_fcn(offspring);
        combined_pop = [pop; offspring]; 
        combined_obj = [obj; off_obj];
        [comb_fronts, ~] = fast_non_dominated_sort(combined_obj);
        next_pop = zeros(popSize, nvars); 
        next_obj = zeros(popSize, size(obj,2));
        fill_count = 0; current_f = 1;
        while fill_count + length(comb_fronts{current_f}) <= popSize
            idx = comb_fronts{current_f}; 
            next_pop(fill_count+1 : fill_count+length(idx), :) = combined_pop(idx, :);
            next_obj(fill_count+1 : fill_count+length(idx), :) = combined_obj(idx, :); 
            fill_count = fill_count + length(idx); 
            current_f = current_f + 1;
        end
        if fill_count < popSize
            idx = comb_fronts{current_f}; 
            last_cd = crowding_distance(combined_obj, {idx});
            [~, sorted_idx] = sort(last_cd(idx), 'descend');
            needed = popSize - fill_count; 
            final_idx = idx(sorted_idx(1:needed));
            next_pop(fill_count+1 : end, :) = combined_pop(final_idx, :); 
            next_obj(fill_count+1 : end, :) = combined_obj(final_idx, :);
        end
        pop = next_pop; obj = next_obj;
    end
    [fronts, ~] = fast_non_dominated_sort(obj); 
    X_pf = pop(fronts{1}, :); 
    Y_pf = obj(fronts{1}, :);
end

function [fronts, rank] = fast_non_dominated_sort(Y)
    N = size(Y, 1); fronts = {[]}; rank = zeros(N, 1); S = cell(N, 1); n = zeros(N, 1);
    for i = 1:N
        for j = 1:N
            if i==j, continue; end
            if all(Y(i,:) <= Y(j,:)) && any(Y(i,:) < Y(j,:)), S{i} = [S{i}, j];
            elseif all(Y(j,:) <= Y(i,:)) && any(Y(j,:) < Y(i,:)), n(i) = n(i) + 1; end
        end
        if n(i) == 0, rank(i) = 1; fronts{1} = [fronts{1}, i]; end
    end
    i = 1;
    while ~isempty(fronts{i})
        Q = []; 
        for p = fronts{i}
            for q = S{p}
                n(q) = n(q) - 1; 
                if n(q) == 0, rank(q) = i + 1; Q = [Q, q]; end
            end
        end
        i = i + 1; fronts{i} = Q;
    end
    fronts(end) = [];
end

function cd = crowding_distance(Y, fronts)
    N = size(Y, 1); cd = zeros(N, 1); M = size(Y, 2);
    for f = 1:length(fronts)
        front = fronts{f}; l = length(front);
        if l <= 2, cd(front) = Inf; continue; end
        for m = 1:M
            [~, idx] = sort(Y(front, m)); 
            cd(front(idx(1))) = Inf; 
            cd(front(idx(end))) = Inf;
            f_min = Y(front(idx(1)), m); f_max = Y(front(idx(end)), m);
            if f_max == f_min, continue; end
            for i = 2:l-1
                cd(front(idx(i))) = cd(front(idx(i))) + (Y(front(idx(i+1)), m) - Y(front(idx(i-1)), m)) / (f_max - f_min); 
            end
        end
    end
end

function st = local_diag(X, y, beta)
    keep = (beta ~= 0); keep(1) = true; p = sum(keep); n = length(y); df = n - p; X_sub = X(:, keep);
    if rank(X_sub) < p
        st.Q2 = -Inf; st.df = df; st.RMSE_loo = Inf; st.e_loo = y; return; 
    end
    invXtX = inv(X_sub' * X_sub); 
    e = y - X_sub * beta(keep); 
    h = sum((X_sub * invXtX) .* X_sub, 2);
    if df > 0 && all(h < 0.999)
        e_loo = e ./ (1 - h); PRESS = sum(e_loo.^2); SST = sum((y - mean(y)).^2); 
        st.Q2 = 1 - PRESS / SST; st.RMSE_loo = sqrt(PRESS / n); 
        st.e_loo = e_loo; st.s = sqrt(sum(e.^2) / df);
    else
        st.Q2 = -Inf; st.RMSE_loo = Inf; st.e_loo = e; st.s = Inf; 
    end
    st.df = df; st.keep = keep; st.invXtX = invXtX;
end

function [best_beta, best_st] = local_refit_press_heredity(X, y, min_df)
    n = size(X, 1); best_PRESS = Inf; best_beta = zeros(14, 1);
    for k = 0:8191
        b = dec2bin(k, 13) - '0'; keep = [true, b == 1]; p = sum(keep); df = n - p;
        if df >= min_df && local_check_heredity(keep, 'strong')
            X_sub = X(:, keep); 
            if rank(X_sub) < p, continue; end
            invXtX = inv(X_sub' * X_sub); 
            beta_sub = invXtX * X_sub' * y; 
            e = y - X_sub * beta_sub; 
            h = sum((X_sub * invXtX) .* X_sub, 2);
            if any(h >= 0.99), continue; end
            PRESS = sum((e ./ (1 - h)).^2);
            if PRESS < best_PRESS
                best_PRESS = PRESS; 
                best_beta = zeros(14, 1); 
                best_beta(keep) = beta_sub; 
            end
        end
    end
    best_st = local_diag(X, y, best_beta);
end

function is_valid = local_check_heredity(keep, mode)
    is_valid = true; x1 = keep(2); x2 = keep(3); x3 = keep(4); x4 = keep(5);
    if keep(6) && ~x1, is_valid = false; return; end
    if keep(7) && ~x2, is_valid = false; return; end
    if keep(8) && ~x3, is_valid = false; return; end 
    if strcmpi(mode, 'strong')
        if keep(9) && (~x1 || ~x2), is_valid = false; return; end
        if keep(10) && (~x1 || ~x3), is_valid = false; return; end
        if keep(11) && (~x2 || ~x3), is_valid = false; return; end 
        if keep(12) && (~x1 || ~x4), is_valid = false; return; end
        if keep(13) && (~x2 || ~x4), is_valid = false; return; end
        if keep(14) && (~x3 || ~x4), is_valid = false; return; end 
    end
end

function t_val = local_tinv_Bonferroni(df)
    t_table = [38.189, 6.577, 3.864, 3.037, 2.661, 2.446, 2.308, 2.213, 2.143, 2.089, 2.046, 2.011, 1.982, 1.957];
    if df >= 1 && df <= 14, t_val = t_table(df); else, t_val = 2.394; end
end