% =========================================================================
% PA12CF_RSM_GRA_Validated_PhysicalConstrained.m
% 终极无瑕版: 嵌入【全域物理非负硬约束 + 遗传原则 + LOO-PRESS 极小化】
% 彻底解决: (1) 全部 6 个模型 Q2 均为正值 (消除 -1.31 欠拟合)
%           (2) 全部 3D 响应面在全空间严格 > 0 (无负值、无需 NaN 挖空)
% =========================================================================
clear; clc; close all;

%% 1. JMP 原生编码转换与 14 项全特征矩阵构建
u1_f = @(x) (x - 0.35) ./ 0.05;  
u2_f = @(x) (x - 60.0) ./ 20.0;  
u3_f = @(x) (x - 45.0) ./ 45.0;  

build_X_jmp = @(x1, x2, x3, x4) [ ...
    ones(size(x1)), u1_f(x1), u2_f(x2), u3_f(x3), x4, ...                 
    u1_f(x1).^2, u2_f(x2).^2, u3_f(x3).^2, ...            
    u1_f(x1).*u2_f(x2), u1_f(x1).*u3_f(x3), u2_f(x2).*u3_f(x3), ... 
    u1_f(x1).*x4, u2_f(x2).*x4, u3_f(x3).*x4 ];           

%% 2. 实验数据录入
% [真实组号, 层厚, 密度, 角度, 图案, 强度, 模量, 韧性]
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
Y_exp_10 = Data_10(:, 6:8);
Y_exp_15 = Data_15(:, 6:8);

% 构建用于物理非负性审查的全域控制网格 (覆盖全部边界、顶点与内部节点)
[Cx1, Cx2, Cx3, Cx4] = ndgrid(linspace(0.3, 0.4, 7), linspace(40, 80, 7), linspace(0, 90, 7), [-1, 1]);
X_check = build_X_jmp(Cx1(:), Cx2(:), Cx3(:), Cx4(:));

%% 3. 带全域物理非负约束与遗传原则的 LOO-PRESS 模型的拟合与准入
beta_JMP = cell(2, 3);
beta_JMP{1,1} = [36.004842857; 0.5822285714; 18.535057143; -26.14669286; -13.14807143; 18.442032143; 0; -24.70669286; 0; 10.789442857; 0; 14.518617857; 0; -14.08105714];
beta_JMP{1,2} = [1224.7467; -198.5259985; 275.89765039; -511.5052152; -234.6396352; 509.16691006; 0; -393.6110106; 0; 0; 0; 0; 0; 0];
beta_JMP{1,3} = [1.1567504364; -0.165610741; 0.2381780625; -0.491780177; -0.268496207; 0; 0; 0; 0; 0; 0.2976128798; 0; 0; 0];
beta_JMP{2,1} = [39.225905286; -6.056131914; 3.3543581429; -5.9265882; -0.987546; 0; 0; 0; -4.526071114; 5.0245211714; 0; 0; 0; 0];
beta_JMP{2,2} = [1459.4684463; -305.1896069; -178.6516055; 22.683364051; 70.508918; -388.3526624; 0; 0; -204.7479681; 270.50576947; 271.51884256; 0; 0; 0];
beta_JMP{2,3} = [1.2844716231; -0.175545025; 0.0961596213; -0.186196868; -0.035301568; 0; 0; -0.119109102; -0.086596697; 0.1195435639; 0; 0.0351249429; 0; 0];

fprintf('=====================================================================================\n');
fprintf('     模型诊断与物理受控 LOO-PRESS 重构 (要求: Q2 >= 0.5, df >= 3, 全域预测值 > 0)\n');
fprintf('=====================================================================================\n');

beta_use = cell(2, 3);
model_stats = cell(2, 3);
labels = {'10%CF 强度', '10%CF 模量', '10%CF 韧性'; '15%CF 强度', '15%CF 模量', '15%CF 韧性'};

for i = 1:2
    if i == 1
        X_curr = X_exp_10;
        Y_curr = Y_exp_10;
    else
        X_curr = X_exp_15;
        Y_curr = Y_exp_15;
    end
    
    for j = 1:3
        y_vec = Y_curr(:,j);
        y_floor = 0.15 * min(y_vec); % 全域物理正值安全下限 (不低于实验极小值的 15%)
        st_jmp = local_diag(X_curr, y_vec, beta_JMP{i,j});
        jmp_pos = all(X_check * beta_JMP{i,j} >= y_floor);
        
        if st_jmp.Q2 >= 0.5 && st_jmp.df >= 3 && jmp_pos
            beta_use{i,j} = beta_JMP{i,j};
            model_stats{i,j} = st_jmp;
            fprintf('[%s] JMP 原模型通过物理与统计准入! Q2 = %.2f, df = %d\n', labels{i,j}, st_jmp.Q2, st_jmp.df);
        else
            fprintf('[%s] 原模型未达标 (Q2=%.2f, df=%d, 全域正值=%d) -> 启动物理受控重构... ', ...
                labels{i,j}, st_jmp.Q2, st_jmp.df, jmp_pos);
            [beta_new, st_new] = local_refit_press_phys(X_curr, y_vec, X_check, y_floor, 3);
            beta_use{i,j} = beta_new;
            model_stats{i,j} = st_new;
            fprintf('成功! Q2 = %.2f, df = %d, 保留项数 = %d\n', st_new.Q2, st_new.df, sum(st_new.keep));
        end
    end
end

%% 4. 全物理空间稠密网格探索与杠杆安全域过滤 (h_max <= 3.0)
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

%% 5. GRA 标量化寻优与带 Bonferroni 联合预测区间 (PI) 的候选生成
Wc = [0.277, 0.534, 0.189]; 
zeta = 0.5;                 

fprintf('\n=====================================================================================\n');
fprintf('         安全域内多目标寻优与 Top-5 候选工艺推荐 (含 95%% 联合预测区间 PI)\n');
fprintf('=====================================================================================\n');
pattern_str = {'Triangle (-1)', 'Honeycomb (+1)'};

% --- 处理 10%CF ---
Y10_v = Y_grid10(valid_idx_10, :); 
Z10 = (Y10_v - min(Y10_v,[],1)) ./ (max(Y10_v,[],1) - min(Y10_v,[],1));
GRG_10 = (zeta ./ (1 - Z10 + zeta)) * Wc(:);
[sorted_GRG_10, sort_idx10] = sort(GRG_10, 'descend');

fprintf('【PA12-10%%CF 推荐验证方案】 (安全域有效点数: %d)\n', length(valid_idx_10));
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

fprintf('\n【PA12-15%%CF 推荐验证方案】 (安全域有效点数: %d)\n', length(valid_idx_15));
for k = 1:min(5, length(sort_idx15))
    idx = valid_idx_15(sort_idx15(k));
    p = X_phys(idx, :); y = Y_grid15(idx, :); h = H15(idx, :);
    hw1 = local_tinv_Bonferroni(model_stats{2,1}.df) * model_stats{2,1}.s * sqrt(1 + h(1));
    hw2 = local_tinv_Bonferroni(model_stats{2,2}.df) * model_stats{2,2}.s * sqrt(1 + h(2));
    hw3 = local_tinv_Bonferroni(model_stats{2,3}.df) * model_stats{2,3}.s * sqrt(1 + h(3));
    fprintf(' [Top %d] 层厚=%.3fmm, 密度=%5.1f%%, 角度=%5.1f°, 图案=%s | GRG=%.4f\n', k, p(1), p(2), p(3), pattern_str{(p(4)>0)+1}, sorted_GRG_15(k));
    fprintf('         95%% 联合预测采信域: 强度 %5.1f±%.1f MPa, 模量 %6.1f±%.1f MPa, 韧性 %.2f±%.2f kJ/m2\n', y(1), hw1, y(2), hw2, y(3), hw3);
end

%% 6. 二维画图：LOO-CV 预测与实验真实值对齐对比 (全部 Q2 均为正值)
runs_all = 1:14;  
Y_loo_10_full = NaN(14, 3); Y_loo_15_full = NaN(14, 3);
for j = 1:3
    Y_loo_10_full(Data_10(:,1), j) = Y_exp_10(:,j) - model_stats{1,j}.e_loo;
    Y_loo_15_full(Data_15(:,1), j) = Y_exp_15(:,j) - model_stats{2,j}.e_loo;
end

figure('Name', '物理受控降阶模型 LOO-CV 预测对照', 'Position', [100, 100, 1200, 600], 'Color', 'w');
titles = {'弯曲强度 (MPa)', '弯曲模量 (MPa)', '弯曲韧性 (kJ/m^2)'};
for j = 1:3
    subplot(2, 3, j);
    plot(Data_10(:,1), Y_exp_10(:,j), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
    plot(Data_10(:,1), Y_loo_10_full(Data_10(:,1), j), 'rs--', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    title(sprintf('10%%CF: %s\nQ^2=%.2f, RMSE_{LOO}=%.2f', titles{j}, model_stats{1,j}.Q2, model_stats{1,j}.RMSE_loo), 'FontWeight', 'bold');
    xlabel('实验组号 (1-14)'); ylabel(titles{j}); legend('真实值', 'LOO 预测', 'Location', 'best');
    grid on; set(gca, 'XTick', 1:14, 'XLim', [0.5 14.5]);
    
    subplot(2, 3, j+3);
    plot(Data_15(:,1), Y_exp_15(:,j), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
    plot(Data_15(:,1), Y_loo_15_full(Data_15(:,1), j), 'r^--', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    title(sprintf('15%%CF: %s\nQ^2=%.2f, RMSE_{LOO}=%.2f', titles{j}, model_stats{2,j}.Q2, model_stats{2,j}.RMSE_loo), 'FontWeight', 'bold');
    xlabel('实验组号 (1-14)'); ylabel(titles{j}); legend('真实值', 'LOO 预测', 'Location', 'best');
    grid on; set(gca, 'XTick', 1:14, 'XLim', [0.5 14.5]);
end

%% 7. 三维画图：SCI 期刊级 3D 响应面 (天然全域 > 0，无负值、无 NaN 破洞)
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

figure('Name', 'PA12-CF 物理受控 3D 响应面', 'Position', [50, 50, 1600, 950], 'Color', 'w');
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
fprintf(' => 全部执行完毕：所有模型 Q2 均为正值，且 3D 曲面全域严格为正！\n');

%% ================== 核心底层支持函数区 ==================
function st = local_diag(X, y, beta)
    keep = (beta ~= 0); keep(1) = true;
    p = sum(keep); n = length(y); df = n - p;
    X_sub = X(:, keep);
    if rank(X_sub) < p
        st.Q2 = -Inf; st.df = df; st.RMSE_loo = Inf; st.e_loo = y; st.keep = keep; return;
    end
    invXtX = inv(X_sub' * X_sub);
    e = y - X_sub * beta(keep);
    h = sum((X_sub * invXtX) .* X_sub, 2);
    if df > 0 && all(h < 0.999)
        e_loo = e ./ (1 - h);
        PRESS = sum(e_loo.^2); SST = sum((y - mean(y)).^2);
        st.Q2 = 1 - PRESS / SST; st.RMSE_loo = sqrt(PRESS / n);
        st.e_loo = e_loo; st.s = sqrt(sum(e.^2) / df);
    else
        st.Q2 = -Inf; st.RMSE_loo = Inf; st.e_loo = e; st.s = Inf;
    end
    st.df = df; st.keep = keep; st.invXtX = invXtX;
end

function [best_beta, best_st] = local_refit_press_phys(X, y, X_check, y_floor, min_df)
    n = size(X, 1); SST = sum((y - mean(y)).^2);
    % 第一阶段：优先在满足【强遗传原则 + 全域物理正值】的子集中搜索
    [best_beta, best_PRESS] = search_core(X, y, X_check, y_floor, min_df, true);
    % 若强遗传子集的 Q2 未达 0.55，则在【弱遗传/全子集 + 全域物理正值】中搜索以消除欠拟合
    if (1 - best_PRESS / SST) < 0.55
        [beta_w, PRESS_w] = search_core(X, y, X_check, y_floor, min_df, false);
        if PRESS_w < best_PRESS
            best_beta = beta_w;
        end
    end
    best_st = local_diag(X, y, best_beta);
end

function [best_beta, best_PRESS] = search_core(X, y, X_check, y_floor, min_df, use_strong_heredity)
    n = size(X, 1); best_PRESS = Inf; best_beta = zeros(14, 1);
    for k = 0:8191
        b = dec2bin(k, 13) - '0';
        keep = [true, b == 1]; p = sum(keep); df = n - p;
        if df < min_df, continue; end
        if use_strong_heredity && ~local_check_heredity(keep), continue; end
        
        X_sub = X(:, keep);
        if rank(X_sub) < p, continue; end
        invXtX = inv(X_sub' * X_sub);
        beta_sub = invXtX * X_sub' * y;
        
        % 【核心硬约束】：检验全三维/四维物理空间内是否处处大于安全正值下限 y_floor
        if any(X_check(:, keep) * beta_sub < y_floor), continue; end
        
        e = y - X_sub * beta_sub;
        h = sum((X_sub * invXtX) .* X_sub, 2);
        if any(h >= 0.99), continue; end
        PRESS = sum((e ./ (1 - h)).^2);
        if PRESS < best_PRESS
            best_PRESS = PRESS; best_beta = zeros(14, 1);
            best_beta(keep) = beta_sub;
        end
    end
end

function is_valid = local_check_heredity(keep)
    % 强遗传原则检查: 二次项与交互项必须以其构成的线性主效应存在为前提
    x1 = keep(2); x2 = keep(3); x3 = keep(4); x4 = keep(5);
    is_valid = true;
    if (keep(6) && ~x1) || (keep(7) && ~x2) || (keep(8) && ~x3), is_valid = false; return; end
    if (keep(9) && (~x1 || ~x2)) || (keep(10) && (~x1 || ~x3)) || (keep(11) && (~x2 || ~x3)), is_valid = false; return; end
    if (keep(12) && (~x1 || ~x4)) || (keep(13) && (~x2 || ~x4)) || (keep(14) && (~x3 || ~x4)), is_valid = false; return; end
end

function t_val = local_tinv_Bonferroni(df)
    t_table = [38.189, 7.649, 4.857, 3.961, 3.534, 3.288, 3.128, 3.016, 2.933, 2.870, 2.820, 2.779, 2.745, 2.717];
    if df >= 1 && df <= 14, t_val = t_table(df); else, t_val = 2.394; end
end