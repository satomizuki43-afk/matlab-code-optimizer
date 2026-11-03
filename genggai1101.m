function results = fdm_paper_analysis(cfg)
% FDM夹层复合材料三点弯曲综合性能分析
% 保存为：fdm_paper_analysis.m
% 运行：results = fdm_paper_analysis;
%
% MATLAB R2020b或更新版本，无需统计或优化工具箱。
%
% 分析顺序：
% 3.1 实验数据
% 3.2 ANOVA与设计可估计性
% 3.3 回归模型与预测验证
% 3.4 残差及影响诊断
% 3.5 多目标综合评价
% 3.6 稳定性与筛选敏感性
% 3.7 相关关系与机理讨论依据
% 3.8 局限、补点与独立确认
%
% 数据：初稿表3、表4中的28组均值。
% 强度、模量：MPa。
% 韧性：按应力-应变积分解释，单位MJ/m^3；
% 原始应变必须为无量纲小数，已有积分数值需要自行核对。
%
% 不自动删除异常点，不虚构逐件SD或纯误差。
% 所有table逐列赋值，不使用Row作为变量名称。

if nargin < 1
    cfg = struct();
end

cfg = settings(cfg);
D = experimentData(cfg.InputCSV);
validateData(D,cfg);

if ~exist(cfg.OutDir,'dir')
    mkdir(cfg.OutDir);
end

writetable(D,fullfile(cfg.OutDir,...
    'S31_all_28_condition_means.csv'));

terms = {'Intercept','Layer','Density','Angle','Pattern',...
    'Layer2','Density2','Angle2','LayerDensity',...
    'LayerAngle','DensityAngle','LayerPattern',...
    'DensityPattern','AnglePattern'};

names = {'FlexStrength','FlexModulus','FlexToughness'};

levels = unique(D.CF)';
totalDiag = table();
summary = table();

results = struct('Config',cfg,'Data',D,...
    'PrimaryScenario','All');

for mi = 1:numel(levels)

    cf = levels(mi);
    B = D(D.CF == cf,:);

    X = features(B,cfg);
    Y = B{:,names};

    tag = sprintf('CF%g',cf);

    fprintf('\n========== %g%% CF ==========\n',cf);

    %% 3.1 实验数据
    % 全部保留；没有逐件结果，不生成平行试件误差条。

    fprintf('3.1 实验数据：%d个工艺条件。\n',height(B));

    dataSummary = table();
    dataSummary.Response = names';
    dataSummary.Minimum = min(Y,[],1)';
    dataSummary.Maximum = max(Y,[],1)';
    dataSummary.MeanAcrossConditions = mean(Y,1)';

    % 注意：这是不同工艺之间的离散程度，不是组内平行试件SD。
    dataSummary.SDAcrossConditions = std(Y,0,1)';

    writetable(dataSummary,...
        fileName(cfg,tag,'S31_data_summary'));

    disp(dataSummary);

    if cfg.MakeFigures
        plotData(B,names,tag,cfg);
    end

    %% 3.2 ANOVA及设计可估计性
    % 完整二阶模型有14列。
    % 先检查秩；不能对不可辨识的完整模型解释全部系数。
    %
    % ANOVA使用预先指定的加性二次模型：
    % 截距 + 4个主效应 + 3个连续因素平方项。
    % 连续因素的线性项和平方项联合进行偏F检验。
    %
    % 误差是工艺均值的模型残差，不是逐件纯误差。

    audit = designAudit(X,terms);

    writetable(audit.Summary,...
        fileName(cfg,tag,'S32_design_rank'));

    writetable(audit.AliasBasis,...
        fileName(cfg,tag,'S32_alias_nullspace'));

    fprintf('3.2 完整二阶设计矩阵：%d列，秩=%d。\n',...
        size(X,2),audit.Rank);

    if audit.Rank < 14
        fprintf('存在混杂，不能独立估计全部二阶系数。\n');
    end

    A = table();

    for j = 1:3
        A = appendTable(A,...
            plannedAnova(X,Y(:,j),names{j}));
    end

    if ~isempty(A)
        A.Q_BH = bhAdjust(A.P_ModelConditional);
    end

    writetable(A,fileName(cfg,tag,'S32_planned_ANOVA'));
    disp(A);

    %% 3.3 回归模型与预测验证
    % 固定主效应模型作为基线。
    % 辅助候选：主效应 + 至多两个二阶项。
    % 所有二阶项保留其主效应，满足强遗传原则。
    %
    % 内层LOO选择模型。
    % 外层每次留出一个完整工艺条件，重新选择模型。
    % 不将同条件平行试件拆分到训练集和验证集中。

    fprintf('3.3 回归与预测验证。\n');

    [models,diagnostics,predictions,frequency] = ...
        fitResponses(B,cfg,terms);

    diagnostics.CF = repmat(cf,3,1);
    diagnostics.Scenario = repmat({'All'},3,1);

    totalDiag = appendTable(totalDiag,diagnostics);

    writetable(diagnostics,...
        fileName(cfg,tag,'S33_model_diagnostics'));

    for j = 1:3

        coef = coefficientTable(models{j},terms);

        writetable(coef,fileName(cfg,tag,...
            ['S33_' names{j} '_coefficients']));

        writetable(predictions{j},fileName(cfg,tag,...
            ['S33_' names{j} '_predictions']));
    end

    disp(diagnostics(:,{'Response','R2','Q2MainEffects',...
        'Q2Nested','RMSENested','Status'}));

    %% 3.4 残差与影响诊断
    % 输出：
    % 实测值-嵌套LOO预测值；
    % 残差-拟合值；
    % Cook距离；
    % 学生化残差正态分位图。
    %
    % 高影响点用于核查，不自动删除。

    fprintf('3.4 残差与影响诊断。\n');

    if cfg.MakeFigures
        drawDiagnostics(B,models,names,tag,cfg);
    end

    %% 3.5 多目标综合评价
    % 先计算实测Pareto集合，再计算GRG。
    % GRG是单一综合评价指标。
    %
    % 全部样本归一化范围作为共同基准，
    % 历史筛选对照使用相同范围与相同权重。
    %
    % 默认权重仅延续原程序，不视为已验证的BWM结果。

    [Z,C,bounds] = greyCoefficients(Y,cfg.Zeta,[]);
    [W,Wcritic] = chooseWeights(Z,cfg);

    wt = table();
    wt.Response = names';
    wt.UsedWeight = W';
    wt.CRITICWeight = Wcritic';

    writetable(wt,fileName(cfg,tag,'S35_weights'));

    ranking = rankMeasured(B,bounds,W,cfg);

    writetable(ranking,...
        fileName(cfg,tag,'S35_measured_ranking'));

    summary = appendTable(summary,...
        scenarioSummary(cf,'All',ranking,X));

    fprintf('3.5 实测Pareto组号：');
    fprintf('%g ',ranking.GroupID(ranking.ParetoEligible));
    fprintf('\n');

    disp(ranking(:,{'GroupID','CandidateEligible',...
        'FlexStrength','FlexModulus','FlexToughness',...
        'GRG_CommonReference','ParetoEligible'}));

    if cfg.MakeFigures
        drawObserved(B,Z,C,W,cf,cfg);
    end

    % 枚举已采用的离散因素水平，生成待确认预测。
    % 图案保持分类编码，不在-1和+1之间插值。
    grid = predictionGrid(B,cfg);

    candidates = predictGrid(...
        grid,B,B,models,bounds,W,cfg);

    % 超出实测归一化范围的预测，不赋予可比较GRG。
    candidates.GRG_ModelReference(...
        candidates.OutsideMeasuredResponseRange) = NaN;

    writetable(candidates,...
        fileName(cfg,tag,'S35_UNCONFIRMED_predictions'));

    %% 3.6 稳定性与筛选敏感性
    % 这里使用离散穷举，不强行生成遗传算法收敛曲线。
    %
    % 分别检查：
    % 权重变化；
    % 外层验证中的模型项选择频率；
    % 历史23组子域；
    % 解除混杂的试验补点。

    fprintf('3.6 稳定性与筛选敏感性。\n');

    eligible = eligibleMask(B,cfg);

    sensitivity = weightSensitivity(...
        C(eligible,:),B.GroupID(eligible));

    writetable(sensitivity,...
        fileName(cfg,tag,'S36_weight_sensitivity'));

    frequencyTable = array2table(frequency,...
        'VariableNames',names);

    frequencyTable.Term = terms';

    writetable(frequencyTable,...
        fileName(cfg,tag,'S36_term_frequency'));

    screened = struct();
    keep = ~ismember(B.GroupID,cfg.LegacyExcludedIDs);

    if any(~keep)

        T = B(keep,:);

        oldRanking = rankMeasured(T,bounds,W,cfg);

        [oldModels,oldDiag] = fitResponses(T,cfg,terms);

        oldDiag.CF = repmat(cf,3,1);
        oldDiag.Scenario = ...
            repmat({'HistoricalSubset'},3,1);

        totalDiag = appendTable(totalDiag,oldDiag);

        writetable(oldRanking,...
            fileName(cfg,tag,'S36_screened_ranking'));

        writetable(oldDiag,...
            fileName(cfg,tag,'S36_screened_diagnostics'));

        summary = appendTable(summary,...
            scenarioSummary(cf,'HistoricalSubset',...
            oldRanking,features(T,cfg)));

        screened = struct('Data',T,...
            'Ranking',oldRanking,...
            'Diagnostics',oldDiag,...
            'Models',{oldModels});
    end

    augmentation = augmentDesign(B,cfg);

    writetable(augmentation,...
        fileName(cfg,tag,'S36_design_augmentation'));

    %% 3.7 相关关系与机理讨论依据
    % 输出响应相关矩阵。
    % 孔隙、熔合颈、纤维取向、界面损伤等机制，
    % 必须结合真实光学/SEM或其他表征人工讨论。
    % 程序不从宏观均值自动生成微观结论。

    correlation = table();
    correlation.Objective = names';

    R = corrcoef(Y);

    for j = 1:3
        correlation.(names{j}) = R(:,j);
    end

    writetable(correlation,...
        fileName(cfg,tag,'S37_response_correlations'));

    fprintf('3.7 已输出响应相关矩阵；相关关系不证明微观机理。\n');

    %% 3.8 局限、补点与独立确认
    % 生成确认模板。
    % 已填写的模板不覆盖。
    % 确认数据独立输入，不参与回归拟合。

    template = confirmationTemplate(...
        B,ranking,augmentation);

    templatePath = fileName(...
        cfg,tag,'S38_confirmation_TEMPLATE');

    if ~isfile(templatePath)
        writetable(template,templatePath);
    end

    confirmation = checkConfirmation(...
        cfg.ConfirmationCSV,cf,B,models,cfg);

    if ~isempty(confirmation)

        writetable(confirmation,...
            fileName(cfg,tag,'S38_confirmation_results'));
    end

    fprintf('3.8 已输出补点建议和确认模板。\n');

    results.Material(mi) = struct(...
        'CF',cf,...
        'Data',B,...
        'Design',audit,...
        'ANOVA',A,...
        'Models',{models},...
        'Diagnostics',diagnostics,...
        'Weights',W,...
        'Ranking',ranking,...
        'WeightSensitivity',sensitivity,...
        'HistoricalSubset',screened,...
        'Augmentation',augmentation,...
        'Confirmation',confirmation);
end

results.Diagnostics = totalDiag;
results.Summary = summary;

writetable(totalDiag,...
    fullfile(cfg.OutDir,'ALL_model_diagnostics.csv'));

writetable(summary,...
    fullfile(cfg.OutDir,'ALL_study_summary.csv'));

save(fullfile(cfg.OutDir,'analysis_results.mat'),'results');

writeLimitations(cfg);

disp(summary);
fprintf('\n输出目录：%s\n',cfg.OutDir);

end


%% ================= 设置与数据 =================

function c = settings(c)

d = struct(...
    'OutDir',fullfile(pwd,'FDM_paper_results'),...
    'InputCSV','',...
    'ConfirmationCSV','',...
    'MakeFigures',true,...
    'ShowFigures',true,...
    'LegacyExcludedIDs',[2 3 6 20 27],...
    'CandidateExcludedIDs',[],...
    'MinTensile',[-Inf -Inf -Inf],...
    'WeightMode','fixed',...
    'Weights',[.277 .534 .189],...
    'SubjectiveWeights',[],...
    'SubjectiveFraction',.5,...
    'Zeta',.5,...
    'MaxParameters',7,...
    'MaxExtraTerms',2,...
    'MinCVResidualDF',3,...
    'MaxAcceptableCVRMSE',[NaN NaN NaN],...
    'AugmentationCount',5,...
    'Center',[.35 60 45],...
    'Scale',[.05 20 45]);

% WeightMode:
% fixed：延续原程序权重。
% equal：等权。
% critic：客观CRITIC权重。
% blend：给定主观权重与CRITIC权重的凸组合。
%
% 尚未提供完整BWM专家偏好数据，因此不伪造BWM计算结果。
%
% MaxAcceptableCVRMSE：
% 按强度、模量、吸能指标的实际工程容许误差填写。
% 保持NaN时，不启用预测候选推荐标志。
%
% MinCVResidualDF、MaxParameters属于建模设置，
% 不是通用SCI判定标准。

assert(isstruct(c) && isscalar(c),...
    'cfg必须是标量结构体。');

assert(all(ismember(fieldnames(c),fieldnames(d))),...
    '存在未知设置字段。');

f = fieldnames(d);

for i = 1:numel(f)
    if ~isfield(c,f{i})
        c.(f{i}) = d.(f{i});
    end
end

assert(c.Zeta > 0 && c.Zeta <= 1,...
    'Zeta须在(0,1]内。');

assert(numel(c.MinTensile) == 3 && ...
    ~any(isnan(c.MinTensile)),...
    '拉伸下限须有3项。');

assert(numel(c.MaxAcceptableCVRMSE) == 3,...
    '预测RMSE限须有3项。');

assert(ismember(c.MaxExtraTerms,0:2),...
    '附加二阶项数须为0、1或2。');

assert(c.MaxParameters >= 5 && ...
    c.MaxParameters <= 14 && ...
    c.MaxParameters == fix(c.MaxParameters),...
    '参数预算无效。');

assert(c.MinCVResidualDF >= 1 && ...
    c.MinCVResidualDF == fix(c.MinCVResidualDF),...
    'CV自由度设置无效。');

assert(numel(c.Center) == 3 && ...
    numel(c.Scale) == 3 && all(c.Scale > 0),...
    '编码设置无效。');

assert(c.AugmentationCount >= 0 && ...
    c.AugmentationCount == fix(c.AugmentationCount),...
    '补点数量无效。');

end


function D = experimentData(file)

% 列顺序：
% ID、CF%、层厚mm、填充%、角度、图案、
% 拉伸强度、拉伸模量、拉伸积分、
% 弯曲强度、弯曲模量、弯曲积分。
%
% 图案：-1=三角形，+1=蜂窝。
% 使用初稿两位小数均值，不与旧高精度子集混用。

M = [
1 10 .30 40 0 -1 23.55 1858.19 1.71 60.92 1914.84 2.06;
2 10 .30 40 90 1 11.85 710.13 1.01 24.40 820.40 .81;
3 10 .35 40 90 1 13.70 892.71 2.70 32.15 1012.47 1.03;
4 10 .40 40 0 1 27.61 2418.60 2.45 42.64 1287.92 1.37;
5 10 .40 40 45 -1 21.95 1269.90 5.66 35.75 1541.34 1.06;
6 10 .30 60 0 1 29.72 2588.98 2.07 66.40 4003.94 2.39;
7 10 .35 60 45 1 18.04 985.00 2.44 22.68 726.38 .69;
8 10 .35 60 45 -1 23.56 1792.28 5.08 48.61 1612.20 1.49;
9 10 .40 60 90 -1 22.00 1183.91 2.24 27.50 779.53 .79;
10 10 .30 80 45 1 18.13 1019.49 2.40 44.82 2037.40 1.44;
11 10 .30 80 90 -1 25.54 1579.82 1.92 52.64 1568.40 1.63;
12 10 .35 80 0 -1 24.96 1863.72 2.47 55.77 1964.09 1.88;
13 10 .40 80 0 -1 25.31 1847.84 1.91 48.22 2001.97 1.63;
14 10 .40 80 90 1 17.53 947.47 4.51 20.83 726.38 .72;
15 15 .30 40 0 -1 32.10 2442.69 1.56 49.46 1789.37 1.53;
16 15 .30 40 90 1 18.92 1114.29 1.44 25.26 903.54 .78;
17 15 .35 40 90 1 23.47 1375.22 1.85 27.35 1470.47 .85;
18 15 .40 40 0 1 35.24 2339.55 2.25 34.16 1222.44 1.06;
19 15 .40 40 45 -1 32.61 1772.22 3.10 35.96 1045.28 1.07;
20 15 .30 60 0 1 44.52 2994.89 2.76 64.13 2285.43 2.06;
21 15 .35 60 45 1 29.97 1771.23 2.14 39.12 1452.76 1.24;
22 15 .35 60 45 -1 33.12 2039.77 2.51 42.84 1474.85 1.35;
23 15 .40 60 90 -1 28.09 1614.63 2.39 32.03 992.13 .94;
24 15 .30 80 45 1 27.97 1718.11 2.39 52.36 1507.27 1.58;
25 15 .30 80 90 -1 31.70 2062.47 1.58 43.23 1328.74 1.28;
26 15 .35 80 0 -1 30.17 2139.50 1.43 47.30 896.76 1.48;
27 15 .40 80 0 -1 25.45 1697.65 1.21 28.62 911.38 .85;
28 15 .40 80 90 1 27.26 1401.15 3.22 30.94 1023.49 .92];

v = {'GroupID','CF','Layer','Density','Angle','Pattern',...
    'TensileStrength','TensileModulus','TensileToughness',...
    'FlexStrength','FlexModulus','FlexToughness'};

if isempty(file)

    D = array2table(M,'VariableNames',v);

else

    D = readtable(file);

    assert(all(ismember(v,D.Properties.VariableNames)),...
        'CSV缺少必要列。');
end

if ~ismember('CandidateEligible',D.Properties.VariableNames)
    D.CandidateEligible = true(height(D),1);
end

end


function validateData(D,c)

v = {'GroupID','CF','Layer','Density','Angle','Pattern',...
    'TensileStrength','TensileModulus','TensileToughness',...
    'FlexStrength','FlexModulus','FlexToughness'};

for i = 1:numel(v)
    assert(isnumeric(D.(v{i})),...
        ['须为数值列：' v{i}]);
end

A = D{:,v};

assert(all(isfinite(A(:))),...
    '输入包含缺失或非有限值。');

assert(numel(unique(D.GroupID)) == height(D),...
    '全局组号不能重复。');

assert(all(ismember(D.Pattern,[-1 1])),...
    '图案须为-1或+1。');

assert(all(ismember(D.CandidateEligible,[0 1])),...
    'CandidateEligible须为0/1。');

assert(all(all(D{:,v(7:12)} > 0)),...
    '力学响应须为正数。');

assert(size(unique(D{:,...
    {'CF','Layer','Density','Angle','Pattern'}},'rows'),1) ...
    == height(D),...
    '同工艺请汇总为一个均值，逐件记录另行保存。');

assert(all(all(abs((D{:,...
    {'Layer','Density','Angle'}}-c.Center)./c.Scale) ...
    <= 1+1e-8)),...
    '数据超出设定因素范围。');

end


%% ================= 回归与综合评价 =================

function [models,diag,P,freq] = fitResponses(T,c,terms)

n = height(T);

% 在外层及内层训练样本中维持统一参数预算。
budget = min(c.MaxParameters,n-2-c.MinCVResidualDF);

assert(budget >= 5,...
    '独立工艺点不足以支持当前参数预算和嵌套验证。');

X = features(T,c);

names = {'FlexStrength','FlexModulus','FlexToughness'};
pool = candidatePool(budget,c.MaxExtraTerms);

models = cell(1,3);
P = cell(1,3);
freq = zeros(14,3);
diag = table();

for j = 1:3

    y = T.(names{j});

    base = fitModel(X,y,1:5);
    assert(base.Valid,'主效应不可估计。');

    m = selectModel(X,y,pool);

    nested = nan(n,1);
    foldTerms = cell(n,1);

    for i = 1:n

        tr = true(n,1);
        tr(i) = false;

        f = selectModel(X(tr,:),y(tr),pool);

        nested(i) = X(i,f.keep)*f.beta;
        foldTerms{i} = f.keep;

        freq(f.keep,j) = freq(f.keep,j)+1/n;
    end

    m.NestedPrediction = nested;
    m.MainEffects = base;

    m.Q2Nested = q2(y,nested);
    m.RMSENested = sqrt(mean((y-nested).^2));

    meanLOO = (sum(y)-y)/(n-1);
    den = sum((y-meanLOO).^2);

    m.Q2VsLOOMean = NaN;

    if den > 0
        m.Q2VsLOOMean = ...
            1-sum((y-nested).^2)/den;
    end

    models{j} = m;

    a = table();
    a.Response = names(j);
    a.N = n;
    a.Parameters = numel(m.keep);
    a.ResidualDF = m.df;
    a.R2 = m.R2;
    a.AdjustedR2 = m.AdjustedR2;

    a.Q2FixedSelected = q2(y,m.loo);
    a.Q2MainEffects = q2(y,base.loo);
    a.Q2Nested = m.Q2Nested;
    a.Q2VsLOOMean = m.Q2VsLOOMean;
    a.RMSENested = m.RMSENested;

    a.MaxLeverage = max(m.h);
    a.ConditionNumber = m.condX;
    a.Terms = {strjoin(terms(m.keep),'+')};

    a.Status = {'Weak_or_unstable'};

    if m.Q2Nested > 0 && m.Q2VsLOOMean > 0
        a.Status = {'Internal_signal_only'};
    end

    limit = c.MaxAcceptableCVRMSE(j);

    if isfinite(limit) && limit > 0 && ...
            m.RMSENested <= limit && ...
            m.Q2VsLOOMean > 0

        a.Status = ...
            {'Internal_error_limit_met_confirmation_needed'};
    end

    diag = appendTable(diag,a);

    p = table();
    p.GroupID = T.GroupID;
    p.Measured = y;
    p.Fitted = m.fit;
    p.MainEffectsLOO = base.loo;
    p.FixedSelectedLOO = m.loo;
    p.NestedLOO = nested;
    p.Residual = m.residual;
    p.Leverage = m.h;
    p.CooksDistance = m.cooks;

    p.OuterFoldTerms = cellfun(...
        @(k)strjoin(terms(k),'+'),...
        foldTerms,'UniformOutput',false);

    P{j} = p;
end

end


function e = eligibleMask(T,c)

e = logical(T.CandidateEligible) & ...
    ~ismember(T.GroupID,c.CandidateExcludedIDs);

v = T{:,{'TensileStrength',...
    'TensileModulus','TensileToughness'}};

e = e & all(v >= c.MinTensile,2);

end


function T = rankMeasured(B,bounds,W,c)

Y = B{:,{'FlexStrength','FlexModulus','FlexToughness'}};

[~,C] = greyCoefficients(Y,c.Zeta,bounds);
[~,Clocal] = greyCoefficients(Y,c.Zeta,[]);

T = B;
T.CandidateEligible = eligibleMask(B,c);

T.GRG_CommonReference = C*W';
T.GRG_LocalReference = Clocal*W';

T.ParetoEligible = false(height(B),1);

e = T.CandidateEligible;

T.ParetoEligible(e) = paretoMax(Y(e,:));

T = sortrows(T,...
    {'CandidateEligible','GRG_CommonReference'},...
    {'descend','descend'});

end


function s = scenarioSummary(cf,name,T,X)

s = table();

s.CF = cf;
s.Scenario = {name};
s.N = height(T);
s.FullQuadraticRank = rank(X);
s.ParetoCount = sum(T.ParetoEligible);

s.BestGroupID = NaN;
s.Strength = NaN;
s.Modulus = NaN;
s.EnergyIntegral = NaN;

k = find(T.CandidateEligible,1);

if ~isempty(k)

    s.BestGroupID = T.GroupID(k);
    s.Strength = T.FlexStrength(k);
    s.Modulus = T.FlexModulus(k);
    s.EnergyIntegral = T.FlexToughness(k);
end

end


%% ================= 独立确认 =================

function V = confirmationTemplate(B,ranking,A)

names = {'Layer','Density','Angle','Pattern'};

k = find(ranking.CandidateEligible,1);

if isempty(k)
    V = B([],names);
else
    V = ranking(k,names);
end

V = [V;A(:,names)];

V.CF = repmat(B.CF(1),height(V),1);
V.FlexStrength = nan(height(V),1);
V.FlexModulus = nan(height(V),1);
V.FlexToughness = nan(height(V),1);

end


function E = checkConfirmation(file,cf,B,models,c)

E = table();

if isempty(file)
    return;
end

assert(isfile(file),'找不到确认数据CSV。');

V = readtable(file);

names = {'FlexStrength','FlexModulus','FlexToughness'};

required = [{'CF','Layer','Density','Angle','Pattern'},names];

assert(all(ismember(required,V.Properties.VariableNames)),...
    '确认CSV列名不完整。');

V = V(V.CF == cf,:);

if isempty(V)
    return;
end

assert(all(ismember(V.Pattern,[-1 1])),...
    '确认数据图案编码无效。');

X = features(V,c);
old = features(B,c);

for i = 1:height(V)

    assert(all(isfinite(X(i,:))),...
        '确认工艺参数存在缺失值。');

    for j = 1:3

        obs = V.(names{j})(i);

        if ~isfinite(obs)
            continue;
        end

        assert(obs > 0,'确认响应须为正数。');

        m = models{j};
        a = X(i,m.keep);

        y = a*m.beta;
        h = sum((a/m.R).^2);

        z = betaincinv(.05,m.df/2,.5);

        half = sqrt(m.df*(1/z-1))*...
            sqrt(m.MSE*(1+h));

        r = table();

        % 不使用Row，避免保留维度名称冲突。
        r.ValidationIndex = i;
        r.Response = names(j);

        r.SeenDesign = any(all(...
            abs(old(:,2:5)-X(i,2:5)) < 1e-10,2));

        r.WithinFactorBounds = ...
            all(abs(X(i,2:4)) <= 1+1e-8);

        r.Measured = obs;
        r.Predicted = y;
        r.RelativeErrorPercent = 100*(y-obs)/obs;

        r.ConditionalPI95Lower = y-half;
        r.ConditionalPI95Upper = y+half;

        r.InConditionalPI95 = ...
            obs >= y-half && obs <= y+half;

        E = appendTable(E,r);
    end
end

end


%% ================= 数值计算 =================

function X = features(T,c)

u = (T{:,{'Layer','Density','Angle'}}-c.Center)./c.Scale;

a = u(:,1);
b = u(:,2);
d = u(:,3);
p = T.Pattern;

X = [ones(height(T),1),u,p,...
    a.^2,b.^2,d.^2,...
    a.*b,a.*d,b.*d,...
    a.*p,b.*p,d.*p];

end


function A = designAudit(X,terms)

A.Rank = rank(X);

[~,s,~] = svd(X,0);
sv = diag(s);

t = table();

t.N = size(X,1);
t.Columns = 14;
t.Rank = A.Rank;

t.FullQuadraticUniqueCoefficients = (A.Rank == 14);
t.MainEffectsRank = rank(X(:,1:5));
t.AdditiveQuadraticRank = rank(X(:,1:8));
t.SmallestSingularValue = sv(end);

A.Summary = t;

N = null(X);

z = table();
z.Term = terms';

for i = 1:size(N,2)
    z.(sprintf('AliasBasis%d',i)) = N(:,i);
end

A.AliasBasis = z;
A.NullSpace = N;

end


function pool = candidatePool(budget,maxExtra)

pool = {1,1:5};
extras = 6:14;

for count = 1:maxExtra

    if 5+count > budget
        continue;
    end

    combinations = nchoosek(extras,count);

    for i = 1:size(combinations,1)

        pool{end+1} = [1:5 combinations(i,:)]; %#ok<AGROW>
    end
end

% 二阶项的主效应全部保留，满足强遗传原则。

end


function m = fitModel(X,y,keep)

A = X(:,keep);

[Q,R] = qr(A,0);

p = numel(keep);
n = numel(y);

m = struct('Valid',false,'keep',keep);

if rank(A) < p || n <= p || cond(R) > 1e8
    return;
end

m.Valid = true;

m.beta = R\(Q'*y);
m.fit = A*m.beta;
m.residual = y-m.fit;

m.h = sum(Q.^2,2);
m.df = n-p;
m.SSE = sum(m.residual.^2);
m.MSE = m.SSE/m.df;

m.condX = cond(R);
m.R = R;

m.loo = nan(n,1);

if all(1-m.h > 1e-9)
    m.loo = y-m.residual./(1-m.h);
end

m.PRESS = sum((y-m.loo).^2);

if any(~isfinite(m.loo))
    m.PRESS = Inf;
end

m.R2 = q2(y,m.fit);
m.AdjustedR2 = 1-(1-m.R2)*(n-1)/m.df;

m.studentized = m.residual./sqrt(...
    max(m.MSE,realmin)*max(1-m.h,realmin));

m.cooks = (m.studentized.^2/p).*...
    m.h./max(1-m.h,realmin);

end


function best = selectModel(X,y,pool)

best = struct('Valid',false,'PRESS',Inf,'keep',1:100);

for i = 1:numel(pool)

    m = fitModel(X,y,pool{i});

    if ~m.Valid || ~isfinite(m.PRESS)
        continue;
    end

    better = ~best.Valid || ...
        m.PRESS < best.PRESS-1e-10*max(1,best.PRESS);

    tied = best.Valid && ...
        abs(m.PRESS-best.PRESS) <= ...
        1e-10*max(1,best.PRESS);

    if better || (tied && numel(m.keep) < numel(best.keep))
        best = m;
    end
end

assert(best.Valid,...
    '没有可估计且具有有效内层LOO预测的候选模型。');

end


function q = q2(y,p)

sst = sum((y-mean(y)).^2);

if sst <= eps(max(1,sum(y.^2)))
    q = NaN;
else
    q = 1-sum((y-p).^2)/sst;
end

end


function T = coefficientTable(m,terms)

V = m.R\eye(numel(m.keep));

se = sqrt(max(m.MSE,0)*sum(V.^2,2));

t = m.beta./max(se,realmin);

p = betainc(...
    m.df./(m.df+t.^2),m.df/2,.5);

T = table();

T.Term = terms(m.keep)';
T.Estimate = m.beta;
T.ConditionalSE = se;
T.ExploratoryP_PostSelection = p;

end


function T = plannedAnova(X,y,response)

keep = 1:8;

full = fitModel(X,y,keep);
T = table();

if ~full.Valid || full.df < 3
    return;
end

groups = {[2 6],[3 7],[4 8],5};

names = {'Layer_linear_and_quadratic',...
    'Density_linear_and_quadratic',...
    'Angle_linear_and_quadratic','Pattern'};

for i = 1:4

    reduced = fitModel(X,y,setdiff(keep,groups{i}));

    if ~reduced.Valid
        continue;
    end

    df = numel(groups{i});

    ss = max(0,reduced.SSE-full.SSE);

    if full.MSE > 0
        F = (ss/df)/full.MSE;
    elseif ss > 0
        F = Inf;
    else
        F = NaN;
    end

    p = betainc(...
        full.df/(full.df+df*F),full.df/2,df/2);

    a = table();

    a.Response = {response};
    a.Factor = names(i);
    a.DF = df;
    a.PartialSS = ss;
    a.F = F;
    a.P_ModelConditional = p;
    a.ResidualDF = full.df;
    a.ResidualSS = full.SSE;

    T = appendTable(T,a);
end

end


function q = bhAdjust(p)

q = nan(size(p));
ok = find(isfinite(p));

if isempty(ok)
    return;
end

[v,ix] = sort(p(ok));

n = numel(v);

a = v.*n./(1:n)';

for i = n-1:-1:1
    a(i) = min(a(i),a(i+1));
end

q(ok(ix)) = min(1,a);

end


function [Z,C,bounds] = greyCoefficients(Y,zeta,bounds)

if isempty(bounds)
    bounds = [min(Y,[],1);max(Y,[],1)];
end

span = bounds(2,:)-bounds(1,:);

spanSafe = span;
spanSafe(span == 0) = 1;

Z = (Y-bounds(1,:))./spanSafe;

Z(:,span == 0) = 1;

C = zeta./(abs(1-Z)+zeta);

end


function [W,wc] = chooseWeights(Z,c)

sd = std(Z,0,1);

R = corrcoef(Z);
R(~isfinite(R)) = 0;
R(1:4:end) = 1;

information = sd.*sum(1-R,2)';
information = max(information,0);

if sum(information) <= 1e-12
    wc = [1 1 1]/3;
else
    wc = information/sum(information);
end

switch lower(c.WeightMode)

    case 'equal'

        W = [1 1 1]/3;

    case 'fixed'

        W = c.Weights;

    case 'critic'

        W = wc;

    case 'blend'

        a = c.SubjectiveFraction;

        assert(a >= 0 && a <= 1,...
            'SubjectiveFraction须在[0,1]内。');

        ws = c.SubjectiveWeights;

        assert(numel(ws) == 3 && all(ws >= 0) && ...
            abs(sum(ws)-1) < 1e-8,...
            '请提供有依据的主观权重。');

        W = a*ws(:)'+(1-a)*wc;

    otherwise

        error('WeightMode须为equal、fixed、critic或blend。');
end

assert(numel(W) == 3 && all(isfinite(W)) && ...
    all(W >= 0) && abs(sum(W)-1) < 1e-8,...
    '权重须非负且和为1。');

W = W(:)';

end


function p = paretoMax(Y)

p = true(size(Y,1),1);

for i = 1:size(Y,1)

    p(i) = ~any(...
        all(Y >= Y(i,:),2) & any(Y > Y(i,:),2));
end

end


function T = weightSensitivity(C,ids)

T = table();

if isempty(ids)
    return;
end

counts = zeros(numel(ids),1);
experiments = 0;

for a = 0:20

    for b = 0:20-a

        w = [a b 20-a-b]/20;

        s = C*w';

        win = abs(s-max(s)) <= 1e-12;

        counts = counts+double(win);
        experiments = experiments+1;
    end
end

T.GroupID = ids;
T.TiedOrSoleWinnerCount = counts;
T.WeightGridPoints = repmat(experiments,numel(ids),1);
T.WinnerFraction = counts/experiments;

% 并列最优时，每个并列组均计入一次。

end


%% ================= 候选预测与补点 =================

function G = predictionGrid(B,c) %#ok<INUSD>

lev = cell(1,3);

names = {'Layer','Density','Angle'};

for k = 1:3
    lev{k} = unique(B.(names{k}));
end

[a,b,d,p] = ndgrid(...
    lev{1},lev{2},lev{3},[-1 1]);

G = table();

G.CF = repmat(B.CF(1),numel(a),1);
G.Layer = a(:);
G.Density = b(:);
G.Angle = d(:);
G.Pattern = p(:);

end


function P = predictGrid(G,T,allData,models,bounds,W,c)

P = G;

Xg = features(G,c);
U = features(T,c);

n = height(G);

P.ExactTrainingCondition = ...
    ismember(Xg(:,2:5),U(:,2:5),'rows');

Xobserved = features(allData,c);

P.PreviouslyMeasured = ...
    ismember(Xg(:,2:5),Xobserved(:,2:5),'rows');

P.PredeclaredEligible = true(n,1);

for i = 1:n

    k = find(all(...
        abs(Xobserved(:,2:5)-Xg(i,2:5)) < 1e-10,2),1);

    if ~isempty(k)

        P.PredeclaredEligible(i) = ...
            logical(allData.CandidateEligible(k)) && ...
            ~ismember(allData.GroupID(k),...
            c.CandidateExcludedIDs);
    end
end

P.WithinPatternHull = false(n,1);

for pat = [-1 1]

    train = U(T.Pattern == pat,2:4);

    ix = G.Pattern == pat;
    query = Xg(ix,2:4);

    inside = ismember(query,train,'rows');

    if size(train,1) >= 4 && ...
            rank(train-mean(train,1)) == 3

        try

            tri = delaunayn(train);

            inside = inside | ...
                ~isnan(tsearchn(train,tri,query));

        catch
            % 凸包数值退化时，仅保留精确训练点的支持标志。
        end
    end

    P.WithinPatternHull(ix) = inside;
end

Yhat = zeros(n,3);
pass = true;

names = {'FlexStrength','FlexModulus','FlexToughness'};

for j = 1:3

    m = models{j};

    A = Xg(:,m.keep);

    Yhat(:,j) = A*m.beta;

    h = sum((A/m.R).^2,2);

    z = betaincinv(.05,m.df/2,.5);
    tc = sqrt(m.df*(1/z-1));

    half = tc*sqrt(m.MSE*(1+h));

    P.(['Predicted_' names{j}]) = Yhat(:,j);

    P.(['ConditionalPI95Lower_' names{j}]) = ...
        Yhat(:,j)-half;

    P.(['ConditionalPI95Upper_' names{j}]) = ...
        Yhat(:,j)+half;

    limit = c.MaxAcceptableCVRMSE(j);

    pass = pass && isfinite(limit) && ...
        limit > 0 && m.RMSENested <= limit && ...
        m.Q2VsLOOMean > 0;
end

P.PositivePredictions = all(Yhat > 0,2);

P.OutsideMeasuredResponseRange = any(...
    Yhat < bounds(1,:) | Yhat > bounds(2,:),2);

[~,C] = greyCoefficients(Yhat,c.Zeta,bounds);

P.GRG_ModelReference = C*W';

P.InternalEngineeringGate = repmat(pass,n,1);

P.TensileFeasibilityVerified = false(n,1);

if ~any(isfinite(c.MinTensile))
    P.TensileFeasibilityVerified = true(n,1);
end

% 有拉伸下限时，只能用实测拉伸结果核实已测工艺。
% 新工艺的弯曲预测不能替代拉伸可行性验证。

if any(isfinite(c.MinTensile))

    required = isfinite(c.MinTensile);

    V = allData{:,{'TensileStrength',...
        'TensileModulus','TensileToughness'}};

    for i = 1:n

        k = find(all(...
            abs(Xobserved(:,2:5)-Xg(i,2:5)) ...
            < 1e-10,2),1);

        if ~isempty(k)

            P.TensileFeasibilityVerified(i) = ...
                all(isfinite(V(k,required))) && ...
                all(V(k,required) >= c.MinTensile(required));
        end
    end
end

P.EligibleForProvisionalRanking = ...
    P.WithinPatternHull & ...
    P.PositivePredictions & ...
    P.PredeclaredEligible & ...
    ~P.OutsideMeasuredResponseRange & ...
    P.InternalEngineeringGate & ...
    P.TensileFeasibilityVerified;

P.ParetoPredictedSupported = false(n,1);

e = P.EligibleForProvisionalRanking;

P.ParetoPredictedSupported(e) = paretoMax(Yhat(e,:));

P.RequiresIndependentConfirmation = true(n,1);

end


function T = augmentDesign(B,c)

G = predictionGrid(B,c);

X = features(B,c);
Q = features(G,c);

available = ~ismember(Q(:,2:5),X(:,2:5),'rows');

T = G([],:);

for step = 1:c.AugmentationCount

    inds = find(available);

    if isempty(inds)
        break;
    end

    ranks = zeros(numel(inds),1);
    score = ranks;

    for k = 1:numel(inds)

        candidate = [X;Q(inds(k),:)];

        s = svd(candidate);
        r = rank(candidate);

        ranks(k) = r;

        score(k) = 2*sum(log(s(1:r)));
    end

    bestRank = max(ranks);

    score(ranks < bestRank) = -Inf;

    [~,k] = max(score);
    ix = inds(k);

    a = G(ix,:);

    a.SequentialOrder = step;
    a.RankAfterAddition = bestRank;

    a.ModelResidualDFIfMeasured = ...
        size(X,1)+1-bestRank;

    T = appendTable(T,a);

    X = [X;Q(ix,:)]; %#ok<AGROW>

    available(ix) = false;
end

end


%% ================= 绘图 =================

function plotData(B,names,tag,c)

visible = 'off';

if c.ShowFigures && usejava('desktop')
    visible = 'on';
end

f = figure('Visible',visible,'Color','w',...
    'Position',[60 60 1150 330]);

for j = 1:3

    subplot(1,3,j);

    plot(B.GroupID,B.(names{j}),'o-',...
        'LineWidth',1);

    xlabel('Group ID');
    ylabel(names{j},'Interpreter','none');
    grid on;
end

print(f,fullfile(c.OutDir,[tag '_S31_data.png']),...
    '-dpng','-r180');

if ~c.ShowFigures || ~usejava('desktop')
    close(f);
end

end


function drawDiagnostics(T,models,names,tag,c)

visible = 'off';

if c.ShowFigures && usejava('desktop')
    visible = 'on';
end

f = figure('Visible',visible,'Color','w',...
    'Name',tag,'Position',[50 50 1250 1100]);

for j = 1:3

    m = models{j};
    y = T.(names{j});

    subplot(4,3,j);

    scatter(y,m.NestedPrediction,30,'filled');
    hold on;

    b = [min([y;m.NestedPrediction]),...
        max([y;m.NestedPrediction])];

    plot(b,b,'k--');

    xlabel('Measured condition mean');
    ylabel('Nested LOO prediction');
    title(names{j},'Interpreter','none');
    grid on;

    subplot(4,3,j+3);

    scatter(m.fit,m.residual,30,'filled');
    yline(0,'k--');

    xlabel('Fitted mean');
    ylabel('Residual');
    grid on;

    subplot(4,3,j+6);

    stem(T.GroupID,m.cooks,'filled');

    xlabel('Group ID');
    ylabel('Cook distance');
    grid on;

    subplot(4,3,j+9);

    r = sort(m.studentized);
    n = numel(r);

    theoretical = -sqrt(2)*...
        erfcinv(2*((1:n)'-.5)/n);

    plot(theoretical,r,'o');

    xlabel('Standard normal quantile');
    ylabel('Studentized residual quantile');
    grid on;
end

print(f,fullfile(c.OutDir,[tag '_diagnostics.png']),...
    '-dpng','-r180');

if ~c.ShowFigures || ~usejava('desktop')
    close(f);
end

end


function drawObserved(B,Z,C,W,cf,c)

visible = 'off';

if c.ShowFigures && usejava('desktop')
    visible = 'on';
end

f = figure('Visible',visible,'Color','w',...
    'Name',sprintf('CF%g measured',cf),...
    'Position',[50 50 1150 420]);

subplot(1,2,1);

imagesc(Z);
colorbar;
caxis([0 1]);

xticks(1:3);
xticklabels({'Strength','Modulus','Energy density'});

yticks(1:height(B));
yticklabels(string(B.GroupID));

ylabel('Group ID');
title('Normalized measured means');

subplot(1,2,2);

y = B{:,{'FlexStrength','FlexModulus','FlexToughness'}};

scatter3(y(:,1),y(:,2),y(:,3),55,C*W','filled');
hold on;

p = false(height(B),1);
e = eligibleMask(B,c);

p(e) = paretoMax(y(e,:));

scatter3(y(p,1),y(p,2),y(p,3),140,'ro');

xlabel('Strength (MPa)');
ylabel('Modulus (MPa)');
zlabel('Energy density (MJ/m^3)');

title('Measured points; circles = eligible Pareto points');

colorbar;
grid on;

print(f,fullfile(c.OutDir,...
    sprintf('CF%g_measured.png',cf)),...
    '-dpng','-r180');

if ~c.ShowFigures || ~usejava('desktop')
    close(f);
end

end


%% ================= 文件与说明 =================

function f = fileName(c,tag,label)

f = fullfile(c.OutDir,[tag '_' label '.csv']);

end


function T = appendTable(T,a)

if isempty(T)
    T = a;
else
    T = [T;a];
end

end


function writeLimitations(c)

fid = fopen(fullfile(c.OutDir,...
    'S37_S38_interpretation.txt'),'w');

assert(fid >= 0,'无法写入说明文件。');

cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>

fprintf(fid,...
    'Condition means only: specimen SD and pure-error lack-of-fit are unavailable.\n');

fprintf(fid,...
    'SDAcrossConditions describes variation across settings, not replicate scatter.\n');

fprintf(fid,...
    'ANOVA uses an additive quadratic model; omitted interactions affect residuals.\n');

fprintf(fid,...
    'Post-selection coefficient p values and intervals exclude selection uncertainty.\n');

fprintf(fid,...
    'Nested CV evaluates model selection; high fitted R2 alone is insufficient.\n');

fprintf(fid,...
    'Historical exclusions define a sensitivity subset, not invalid measurements.\n');

fprintf(fid,...
    'Common-reference normalization and weights are fixed within each material.\n');

fprintf(fid,...
    'Separately normalized GRG scores cannot directly compare materials.\n');

fprintf(fid,...
    'Fixed weights reproduce an earlier setting; BWM provenance is not inferred.\n');

fprintf(fid,...
    'Stress(MPa)*strain(fraction) integral has units MJ/m^3; check percent strain.\n');

fprintf(fid,...
    'Use consistent integration endpoints. The flexural integral is not fracture toughness.\n');

fprintf(fid,...
    'Check actual thickness when changing layer height at a fixed layer count.\n');

fprintf(fid,...
    'Without mass/density, do not claim specific properties or lightweight superiority.\n');

fprintf(fid,...
    'Microscopic mechanisms require optical/SEM evidence, not response correlation alone.\n');

fprintf(fid,...
    'The finite grid is enumerated; continuous global optimality is not established.\n');

fprintf(fid,...
    'Confirmation specimens must be independently printed/tested and excluded from fitting.\n');

fprintf(fid,...
    'Fill the template, save separately, and specify cfg.ConfirmationCSV.\n');

end