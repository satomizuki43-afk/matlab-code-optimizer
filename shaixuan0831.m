% =========================================================================
% Stage 1: 单轴拉伸本构连续性与夹层同质性数据驱动过滤程序
% =========================================================================
clear; clc; close all;

%% 0. 录入 28 组更新后的原始实验数据
% 列顺序: 拉伸强度, 拉伸模量, 拉伸韧性, 弯曲强度, 弯曲模量, 弯曲韧性
data = [
    23.5496, 1858.1930, 1.70580, 60.9202, 1914.8370, 2.05836;
    11.8469,  710.1259, 1.01367, 24.4026,  820.4030, 0.81116;
    13.6999,  892.7084, 2.69993, 32.1503, 1012.4670, 1.03218;
    27.6140, 2418.5966, 2.45303, 42.6411, 1287.9225, 1.37059;
    21.9482, 1269.8994, 5.66375, 35.7520, 1541.3380, 1.06190;
    29.7193, 2588.9786, 2.06761, 66.4016, 4003.9360, 2.38646;
    18.0449,  984.9952, 2.44300, 22.6772,  726.3775, 0.68853;  % Updated
    23.5591, 1792.2800, 5.07516, 48.6142, 1612.2040, 1.49201;
    21.9963, 1183.9078, 2.23551, 27.4961,  779.5275, 0.79389;  % Updated
    18.1347, 1019.4948, 2.40178, 44.8228, 2037.4010, 1.43979;
    25.5376, 1579.8182, 1.92157, 52.6393, 1568.3960, 1.62597;
    24.9583, 1863.7156, 2.46945, 55.7652, 1964.0900, 1.88028;
    25.3139, 1847.8375, 1.91050, 48.2244, 2001.9680, 1.62866;  % Updated
    17.5295,  947.4689, 4.50832, 20.8346,  726.3775, 0.72190;  % Updated
    32.0983, 2442.6917, 1.55895, 49.4646, 1789.3695, 1.52946;  % Updated
    18.9244, 1114.2925, 1.43821, 25.2638,  903.5430, 0.77829;  % Updated
    23.4652, 1375.2168, 1.84835, 27.3543, 1470.4720, 0.85436;
    35.2360, 2339.5475, 2.25279, 34.1575, 1222.4406, 1.05698;
    32.6064, 1772.2196, 3.09860, 35.9646, 1045.2755, 1.07399;
    44.5228, 2994.8940, 2.76372, 64.1338, 2285.4326, 2.06174;
    29.9705, 1771.2347, 2.13844, 39.1181, 1452.7555, 1.24057;
    33.1223, 2039.7726, 2.50795, 42.8392, 1474.8466, 1.34926;
    28.0917, 1614.6282, 2.39186, 32.0315,  992.1255, 0.93684;  % Updated
    27.9652, 1718.1123, 2.39412, 52.3648, 1507.2655, 1.57701;
    31.6993, 2062.4697, 1.57608, 43.2283, 1328.7400, 1.28116;  % Updated
    30.1741, 2139.5008, 1.43408, 47.3031,  896.7620, 1.47763;
    25.4459, 1697.6526, 1.21201, 28.6202,  911.3840, 0.84979;
    27.2578, 1401.1463, 3.21805, 30.9396, 1023.4895, 0.92327
];

N = size(data, 1);
original_IDs = (1:N)';
sigma_t = data(:, 1);
E_t     = data(:, 2);
U_t     = data(:, 3);

%% 1. 按照三联不等式计算准入阈值
% 抗拉强度下限
mu_sigma    = mean(sigma_t);
s_sigma     = std(sigma_t, 0);
limit_S_low = mu_sigma - 1.28 * s_sigma;

% 拉伸模量四分位距双侧同质窗
E_sorted     = sort(E_t, 'ascend');
Q1_E         = 0.5 * (E_sorted(7)  + E_sorted(8));
Q3_E         = 0.5 * (E_sorted(21) + E_sorted(22));
IQR_E        = Q3_E - Q1_E;
limit_E_low  = Q1_E - 0.5 * IQR_E;
limit_E_high = Q3_E + 0.5 * IQR_E;

% 拉伸韧性下限
mu_U        = mean(U_t);
s_U         = std(U_t, 0);
limit_U_low = mu_U - 1.00 * s_U;

%% 2. 样本状态判别
idx_S     = sigma_t >= limit_S_low;
idx_E_low = E_t >= limit_E_low;
idx_E_up  = E_t <= limit_E_high;
idx_E     = idx_E_low & idx_E_up;
idx_U     = U_t >= limit_U_low;

valid_idx    = idx_S & idx_E & idx_U;
defect_idx   = (~idx_S) | (~idx_E_low) | (~idx_U);
baseline_idx = (~idx_E_up) & idx_S & idx_U;

%% 3. 命令行规范表头输出
fprintf('====================================================================================\n');
fprintf(' %-6s | %-14s | %-14s | %-14s | %-20s\n', ...
        'Run ID', 'Sigma_t (MPa)', 'E_t (MPa)', 'U_t (kJ/m^2)', 'Screening Status');
fprintf('------------------------------------------------------------------------------------\n');
for i = 1:N
    if valid_idx(i)
        status_str = 'Qualified';
    elseif defect_idx(i)
        status_str = 'Defect Eliminated';
    else
        status_str = 'Solidified Baseline';
    end
    fprintf(' %-6d | %-14.4f | %-14.4f | %-14.4f | %-20s\n', ...
            original_IDs(i), sigma_t(i), E_t(i), U_t(i), status_str);
end
fprintf('====================================================================================\n');
fprintf(' Excluded Run IDs: '); fprintf('%d ', original_IDs(~valid_idx)); fprintf('\n\n');

%% 4. 输出对接后续主程序的弯曲参数矩阵
filtered_IDs = original_IDs(valid_idx);
flex_data    = data(valid_idx, 4:6);

fprintf('qualified_data = [\n');
for i = 1:length(filtered_IDs)
    fprintf('    %2d, %8.4f, %9.4f, %8.5f;\n', filtered_IDs(i), flex_data(i,1), flex_data(i,2), flex_data(i,3));
end
fprintf('];\n');

%% 5. 绘制规范化表头与精简图例的准入诊断图
figure('Name', 'Stage 1: Tensile Feasibility Screening', 'Color', 'w', 'Position', [100, 150, 1450, 420]);

% 规范化子图标题与坐标轴表头
titles_str = {'(a) Tensile Strength', ...
              '(b) Tensile Modulus', ...
              '(c) Tensile Toughness'};
y_lbls = {'Tensile Strength, \sigma_t (MPa)', ...
          'Tensile Modulus, E_t (MPa)', ...
          'Tensile Toughness, U_t (kJ/m^2)'};

for k = 1:3
    subplot(1, 3, k); hold on; grid on; box on;
    yk = data(:, k);
    
    if k == 1
        y_min_plot = 8; y_max_plot = 48;
        patch([0.2, 28.8, 28.8, 0.2], [limit_S_low, limit_S_low, y_max_plot, y_max_plot], ...
              [0.88, 0.96, 0.88], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Feasible Domain');
        yline(limit_S_low, 'r--', 'LineWidth', 1.6, 'DisplayName', sprintf('\\mu_\\sigma - 1.28s_\\sigma = %.2f', limit_S_low));
    elseif k == 2
        y_min_plot = 500; y_max_plot = 3200;
        patch([0.2, 28.8, 28.8, 0.2], [limit_E_low, limit_E_low, limit_E_high, limit_E_high], ...
              [0.88, 0.96, 0.88], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Feasible Domain');
        yline(limit_E_low,  'r--', 'LineWidth', 1.6, 'DisplayName', sprintf('Q_1 - 0.5IQR = %.1f', limit_E_low));
        yline(limit_E_high, 'b--', 'LineWidth', 1.6, 'DisplayName', sprintf('Q_3 + 0.5IQR = %.1f', limit_E_high));
    else
        y_min_plot = 0.6; y_max_plot = 6.2;
        patch([0.2, 28.8, 28.8, 0.2], [limit_U_low, limit_U_low, y_max_plot, y_max_plot], ...
              [0.88, 0.96, 0.88], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Feasible Domain');
        yline(limit_U_low, 'r--', 'LineWidth', 1.6, 'DisplayName', sprintf('\\mu_U - 1.00s_U = %.2f', limit_U_low));
    end
    
    xline(14.5, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    
    % 绘制散点，图例去除全部括号说明
    scatter(original_IDs(valid_idx), yk(valid_idx), 50, [0.12, 0.45, 0.70], 'filled', ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'Qualified Samples');
    scatter(original_IDs(defect_idx), yk(defect_idx), 75, [0.85, 0.15, 0.15], 's', 'filled', ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'Defect Eliminated');
    scatter(original_IDs(baseline_idx), yk(baseline_idx), 85, [0.95, 0.60, 0.10], 'd', 'filled', ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'Solidified Baseline');
    
    for id_ex = original_IDs(~valid_idx)'
        text(id_ex + 0.5, yk(id_ex), sprintf('ID %d', id_ex), 'FontName', 'Times New Roman', ...
             'FontSize', 9.5, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
    end
    
    xlim([0.2, 28.8]); ylim([y_min_plot, y_max_plot]);
    xticks([1, 7, 14, 21, 28]);
    xlabel('Experimental Run ID', 'FontWeight', 'bold');
    ylabel(y_lbls{k}, 'FontWeight', 'bold');
    title(titles_str{k}, 'FontWeight', 'bold', 'FontSize', 12);
    legend('Location', 'northwest', 'FontSize', 9, 'Box', 'on');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10.5, 'LineWidth', 1.0);
    hold off;
end