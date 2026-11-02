%% PA12-CF 两阶段优化：拉伸筛选后的弯曲响应建模与分类图案验证
% 说明：
% 1) 所有经核查的实验数据均保留，不自动删除异常点；
% 2) 图案采用效应编码：-1=三角形，+1=蜂窝；不在两者之间插值；
% 3) 数据为拉伸筛选后保留的弯曲候选工艺，不将淘汰组视为随机缺失；
% 4) 限制模型项数和残差自由度，避免超小样本过拟合；
% 5) 用PRESS选择模型，同时报告固定模型LOO和嵌套LOO；
% 6) 正值约束只在保留候选工艺上检查，不作为数据修正；
% 7) 已测试保留组的最终GRG采用实测值，模型GRG仅作辅助；
% 8) 输出设计矩阵条件数、杠杆值和95%预测区间；
% 9) 确认试验数据单独输入，不参与模型拟合。

clear; clc; close all;

%% 1. 数据与编码
% 特征顺序：
% [1, layer, density, angle, pattern, layer^2, density^2, angle^2,
%  layer*density, layer*angle, density*angle,
%  layer*pattern, density*pattern, angle*pattern]
% buildX在文件末尾定义为局部函数。脚本中不能同时存在同名变量和局部函数，
% 因此这里不再定义同名匿名函数。

termNames = {'Intercept','Layer','Density','Angle','Pattern', ...
    'Layer^2','Density^2','Angle^2', ...
    'Layer*Density','Layer*Angle','Density*Angle', ...
    'Layer*Pattern','Density*Pattern','Angle*Pattern'};

% 图案编码：-1=三角形，+1=蜂窝
patternNames = ["Triangle", "Honeycomb"];

% [组号, 层厚, 密度, 角度, 图案, 弯曲强度, 弯曲模量, 弯曲韧性]
% 下面数据为拉伸筛选后进入弯曲阶段的候选工艺，全部保留。
Data10 = [
    1,  .30, 40,  0, -1, 60.9202, 1914.8370, 2.058360;
    4,  .40, 40,  0,  1, 42.6411, 1287.9225, 1.370590;
    5,  .40, 40, 45, -1, 35.7520, 1541.3380, 1.061900;
    7,  .35, 60, 45,  1, 22.6772,  726.3775, 0.688530;
    8,  .35, 60, 45, -1, 48.6142, 1612.2040, 1.492010;
    9,  .40, 60, 90, -1, 27.4961,  779.5275, 0.793894;
    10, .30, 80, 45,  1, 44.8228, 2037.4010, 1.439790;
    11, .30, 80, 90, -1, 52.6393, 1568.3960, 1.625970;
    12, .35, 80,  0, -1, 55.7652, 1964.0900, 1.880280;
    13, .40, 80,  0, -1, 48.2244, 2001.9680, 1.628660;
    14, .40, 80, 90,  1, 20.8346,  726.3775, 0.721900];

Data15 = [
    1,  .30, 40,  0, -1, 49.4646, 1789.3695, 1.529460;
    2,  .30, 40, 90,  1, 25.2638,  903.5430, 0.778292;
    3,  .35, 40, 90,  1, 27.3543, 1470.4720, 0.854362;
    4,  .40, 40,  0,  1, 34.1575, 1222.4406, 1.056980;
    5,  .40, 40, 45, -1, 35.9646, 1045.2755, 1.073990;
    7,  .35, 60, 45,  1, 39.1181, 1452.7555, 1.240570;
    8,  .35, 60, 45, -1, 42.8392, 1474.8466, 1.349260;
    9,  .40, 60, 90, -1, 32.0315,  992.1255, 0.936843;
    10, .30, 80, 45,  1, 52.3648, 1507.2655, 1.577010;
    11, .30, 80, 90, -1, 43.2283, 1328.7400, 1.281160;
    12, .35, 80,  0, -1, 47.3031,  896.7620, 1.477630;
    14, .40, 80, 90,  1, 30.9396, 1023.4895, 0.923273];

assert(all(ismember(Data10(:,5),[-1,1])) && all(ismember(Data15(:,5),[-1,1])), ...
    '图案编码必须为-1（三角形）或+1（蜂窝）。');

fprintf('10%% CF：拉伸筛选后进入弯曲建模的候选点数：%d\n',size(Data10,1));
fprintf('15%% CF：拉伸筛选后进入弯曲建模的候选点数：%d\n',size(Data15,1));
fprintf('说明：被拉伸试验筛除的工艺不作为弯曲阶段的缺失数据，也不进入弯曲优化候选域。\n\n');

X10 = buildX(Data10(:,2),Data10(:,3),Data10(:,4),Data10(:,5));
X15 = buildX(Data15(:,2),Data15(:,3),Data15(:,4),Data15(:,5));
Y10 = Data10(:,6:8);
Y15 = Data15(:,6:8);

% 物理审查点：只检查拉伸筛选后保留的候选工艺，避免被淘汰组合影响弯曲模型。
Xcheck = {X10, X15};

%% 2. 模型筛选与诊断
% 超小样本约束：最终模型至少保留5个残差自由度，且最多保留6个参数项。
% 这些限制用于防止少量设计点下保留过多二次项和交互项。
minDf = 5;
maxTerms = 6;
enforcePositive = true;
enforceHeredity = true;
nResp = 3;
respNames = {'Flexural strength (MPa)','Flexural modulus (MPa)', ...
    'Flexural toughness (kJ/m^2)'};

betaUse = cell(2,nResp);
statsUse = cell(2,nResp);
datasets = {X10,Y10; X15,Y15};
materialNames = {'10% CF','15% CF'};

for m = 1:2
    X = datasets{m,1}; Y = datasets{m,2};
    for j = 1:nResp
        y = Y(:,j);
        [b,st] = selectModel(X,y,Xcheck{m},minDf,maxTerms,enforcePositive,enforceHeredity,0);
        if ~st.selected
            warning('%s %s未找到满足正值约束的模型，改用无正值约束模型。', ...
                materialNames{m},respNames{j});
            [b,st] = selectModel(X,y,Xcheck{m},minDf,maxTerms,false,enforceHeredity,0);
        end

        % 固定最终模型的PRESS诊断
        st = diagnoseModel(X,y,b,st.keep,termNames);

        % 严格嵌套LOO：每个外层训练集重新筛选模型
        [q2Nested,rmseNested,yNested] = nestedLOO( ...
            X,y,Xcheck{m},minDf,maxTerms,enforcePositive,enforceHeredity);
        st.Q2_nested = q2Nested;
        st.RMSE_nested = rmseNested;
        st.y_nested = yNested;
        st.Q2_status = q2Nested >= 0.50;

        betaUse{m,j} = b;
        statsUse{m,j} = st;

        fprintf('%s %-24s: terms=%d, df=%d, rank=%d, cond=%.2e, Q2_fixed=%.3f, Q2_nested=%.3f, RMSE_nested=%.3f', ...
            materialNames{m},respNames{j},sum(st.keep),st.df,st.rank,st.condNumber, ...
            st.Q2,st.Q2_nested,st.RMSE_nested);
        if st.Q2_status
            fprintf('  [predictive status: acceptable]\n');
        else
            fprintf('  [predictive status: weak]\n');
        end
    end
end

%% 3. 输出最终保留项及p值（p值仅用于解释，不作为唯一选模标准）
fprintf('\n===== 最终模型保留项与p值 =====\n');
for m = 1:2
    for j = 1:nResp
        st = statsUse{m,j};
        fprintf('\n%s - %s\n',materialNames{m},respNames{j});
        fprintf('  diagnostics: n=%d, terms=%d, df=%d, rank=%d, cond=%.2e, max leverage=%.3f, AICc=%.3f\n', ...
            st.n,st.p,st.df,st.rank,st.condNumber,max(st.leverage),st.AICc);
        fprintf('  note: p-values are exploratory post-selection values, not confirmatory inference.\n');
        for k = find(st.keep(:))'
            fprintf('  %-18s beta=% .6g, p=% .4g\n',termNames{k},st.beta(k),st.pValue(k));
        end
    end
end

%% 4. LOO预测图与残差诊断图
figure('Color','w','Name','LOO prediction and residual diagnostics','Position',[80 80 1350 820]);
for m = 1:2
    X = datasets{m,1}; Y = datasets{m,2};
    for j = 1:nResp
        st = statsUse{m,j};
        subplot(4,3,(m-1)*6+j);
        plot(1:length(Y),Y(:,j),'ko-','LineWidth',1.2,'MarkerFaceColor','k'); hold on;
        plot(1:length(Y),st.y_nested,'rs--','LineWidth',1.2,'MarkerFaceColor','w');
        xlabel('Design point'); ylabel(respNames{j});
        title(sprintf('%s %s: nested Q^2=%.2f, RMSE=%.3g', ...
            materialNames{m},respNames{j},st.Q2_nested,st.RMSE_nested),'Interpreter','none');
        legend('Measured','Nested LOO','Location','best'); grid on;

        subplot(4,3,(m-1)*6+3+j);
        plot(st.yFit,st.residual,'ko','MarkerFaceColor',[.2 .5 .8]); hold on;
        yline(0,'k-'); yline(2*st.sigma,'k--'); yline(-2*st.sigma,'k--');
        xlabel('Fitted value'); ylabel('Residual');
        title(sprintf('%s %s: residual vs fitted',materialNames{m},respNames{j}), ...
            'Interpreter','none'); grid on;
    end
end

%% 5. 拉伸筛选后候选工艺内的多响应搜索
% 当前权重应在优化前预先确定；不要仅为了改变最优点而事后调权重。
% 已有弯曲实测组的最终排序采用实测GRG；模型GRG只用于比较模型趋势。
% 不再对完整连续矩形域进行优化，避免产生未经拉伸筛选和弯曲确认的组合。
W = [0.277,0.534,0.189];
zeta = 0.5;
candidateData = {Data10,Data15};

for m = 1:2
    Dcandidate = candidateData{m};
    [~,ia] = unique(Dcandidate(:,2:5),'rows','stable');
    Dcandidate = Dcandidate(ia,:);
    candidateID = Dcandidate(:,1);
    XphysRaw = Dcandidate(:,2:5);
    Xgrid = buildX(XphysRaw(:,1),XphysRaw(:,2),XphysRaw(:,3),XphysRaw(:,4));
    Ymodel = zeros(size(Xgrid,1),nResp);
    Hgrid = zeros(size(Xgrid,1),nResp);
    for j = 1:nResp
        b = betaUse{m,j}; st = statsUse{m,j};
        Ymodel(:,j) = Xgrid(:,st.keep)*b(st.keep);
        Hgrid(:,j) = sum((Xgrid(:,st.keep)*st.invXtX).*Xgrid(:,st.keep),2);
    end

    % 仅用物理可行且杠杆值不极端的候选点进行模型GRG参考排序。
    valid = all(Ymodel>=0,2) & max(Hgrid,[],2)<=3.0;
    fprintf('\n%s拉伸筛选后候选域有效点数：%d/%d\n', ...
        materialNames{m},sum(valid),size(Xgrid,1));

    % ---------- A. 实测GRG：用于最终已测试候选组排序 ----------
    Ymeas = Dcandidate(:,6:8);
    [GRGmeas,~,~] = computeGRG(Ymeas,W,zeta);
    [GMs,ixm] = sort(GRGmeas,'descend');
    fprintf('%s Top-5实测GRG保留工艺（最终排序依据）：\n',materialNames{m});
    for k = 1:min(5,length(ixm))
        p = Dcandidate(ixm(k),2:5);
        label = patternLabel(p(4));
        fprintf('Measured Top%d (ID%d): layer=%.3f mm, density=%.1f%%, angle=%.1f deg, pattern=%s, GRG=%.4f\n', ...
            k,Dcandidate(ixm(k),1),p(1),p(2),p(3),label,GMs(k));
        fprintf('     measured: strength=%.3f, modulus=%.3f, toughness=%.3f\n', ...
            Ymeas(ixm(k),1),Ymeas(ixm(k),2),Ymeas(ixm(k),3));
    end

    % ---------- B. 模型GRG：仅作为辅助趋势参考 ----------
    Yv = Ymodel(valid,:);
    Xv = XphysRaw(valid,:);
    ymin = min(Yv,[],1); ymax = max(Yv,[],1);
    scale = max(ymax-ymin,eps);
    Z = (Yv-ymin)./scale;
    GRC = zeta./(1-Z+zeta);
    GRG = GRC*W(:);
    [GRGsort,ix] = sort(GRG,'descend');

    fprintf('%s Top-5模型GRG参考排序（非最终结论）：\n',materialNames{m});
    for k = 1:min(5,length(ix))
        p = Xv(ix(k),:);
        id = candidateID(valid);
        id = id(ix(k));
        label = patternLabel(p(4));
        fprintf('Model Top%d (ID%d): layer=%.3f mm, density=%.1f%%, angle=%.1f deg, pattern=%s, GRG=%.4f\n', ...
            k,id,p(1),p(2),p(3),label,GRGsort(k));
        fprintf('     prediction: strength=%.3f, modulus=%.3f, toughness=%.3f\n', ...
            Yv(ix(k),1),Yv(ix(k),2),Yv(ix(k),3));
    end
end

%% 6. 分类变量的直接主效应计算
% 在指定连续参数下比较蜂窝(+1)与三角形(-1)，不对图案进行连续插值。
t0 = .30; rho0 = 80; ang0 = 0;
fprintf('\n===== 指定工艺点的图案效应 =====\n');
for m = 1:2
    fprintf('%s: layer=%.2f, density=%.1f, angle=%.1f\n', ...
        materialNames{m},t0,rho0,ang0);
    for j = 1:nResp
        st = statsUse{m,j}; b = betaUse{m,j};
        Xtri = buildX(t0,rho0,ang0,-1);
        Xhon = buildX(t0,rho0,ang0,+1);
        ytri = Xtri(:,st.keep)*b(st.keep);
        yhon = Xhon(:,st.keep)*b(st.keep);
        fprintf('  %-24s Triangle=%.4f, Honeycomb=%.4f, Delta=%.4f\n', ...
            respNames{j},ytri,yhon,yhon-ytri);
    end
end

%% 7. 外部确认试验接口
% 将NaN替换为确认试验的5根平行试件平均值。
% 当前若只有三角形确认试验，则只能验证三角形工艺点；
% 若要验证图案主效应，应增加同一工艺点下的蜂窝对照。
% 列顺序：[layer density angle pattern strength modulus toughness]
Validation10 = [ ...
    .30,80,0,-1,NaN,NaN,NaN; ... % 三角形，已完成的确认点
    .30,80,0,+1,NaN,NaN,NaN];    % 蜂窝，若要验证图案效应需补做
Validation15 = [ ...
    .30,80,0,-1,NaN,NaN,NaN; ... % 三角形，已完成的确认点
    .30,80,0,+1,NaN,NaN,NaN];    % 蜂窝，若要验证图案效应需补做

fprintf('\n===== 外部确认试验预测 =====\n');
validateTable(Validation10,1,@buildX,betaUse,statsUse,respNames);
validateTable(Validation15,2,@buildX,betaUse,statsUse,respNames);
fprintf('\n===== 外部确认试验实测GRG（仅在填入实测均值后计算） =====\n');
reportValidationGRG(Validation10,Data10,'10% CF',W,zeta);
reportValidationGRG(Validation15,Data15,'15% CF',W,zeta);

%% 8. 条件响应面可视化：仅展示固定角度0°的切片
% 该图用于展示保留工艺子空间内的趋势，不用于重新扩大弯曲优化域。
% 最终候选排序以第5节的实际保留工艺点为准。
[S1,S2] = meshgrid(linspace(.30,.40,40),linspace(40,80,40));
Xtri = buildX(S1(:),S2(:),zeros(numel(S1),1),-ones(numel(S1),1));
Xhon = buildX(S1(:),S2(:),zeros(numel(S1),1), ones(numel(S1),1));

figure('Color','w','Name','Discrete pattern response surfaces','Position',[50 50 1550 900]);
for m = 1:2
    for j = 1:nResp
        st = statsUse{m,j}; b = betaUse{m,j};
        ztri = reshape(Xtri(:,st.keep)*b(st.keep),size(S1));
        zhon = reshape(Xhon(:,st.keep)*b(st.keep),size(S1));
        k = (m-1)*nResp+j;
        subplot(2,3,k);
        surf(S1,S2,ztri,'EdgeColor','none','FaceAlpha',.9); hold on;
        mesh(S1,S2,zhon,'FaceColor','none','EdgeColor',[.1 .1 .1]);
        xlabel('Layer thickness (mm)'); ylabel('Density (%)'); zlabel(respNames{j});
        title(sprintf('%s: conditional %s, angle=0 deg',materialNames{m},respNames{j}), ...
            'Interpreter','none');
        legend('Triangle (-1)','Honeycomb (+1)','Location','best'); grid on; view([-40 25]);
    end
end

fprintf('\n完成：所有拉伸筛选后保留的弯曲数据均已使用；图案仅按-1/+1两个离散水平处理；\n');
fprintf('已测试候选组的最终排序采用实测GRG，模型GRG仅作辅助；确认试验与模型拟合保持分离。\n');

%% ====================== 局部函数 ======================
function X = buildX(x1,x2,x3,x4)
    u1 = (x1-.35)./.05;
    u2 = (x2-60)./20;
    u3 = (x3-45)./45;
    X = [ones(size(x1)),u1,u2,u3,x4,u1.^2,u2.^2,u3.^2, ...
        u1.*u2,u1.*u3,u2.*u3,u1.*x4,u2.*x4,u3.*x4];
end

function label = patternLabel(p)
    if p == -1
        label = 'Triangle';
    elseif p == 1
        label = 'Honeycomb';
    else
        error('Invalid pattern code: %g',p);
    end
end

function [beta,st] = selectModel(X,y,Xcheck,minDf,maxTerms,enforcePositive,enforceHeredity,yFloor)
    n = size(X,1); nTerms = size(X,2);
    bestPRESS = Inf; beta = zeros(nTerms,1); bestKeep = false(nTerms,1);
    for mask = 0:(2^(nTerms-1)-1)
        bits = dec2bin(mask,nTerms-1)-'0';
        keep = [true, logical(bits)];
        p = sum(keep); df = n-p;
        if df < minDf, continue; end
        if p > maxTerms, continue; end
        if enforceHeredity && ~checkHeredity(keep), continue; end
        Xs = X(:,keep);
        if rank(Xs) < p, continue; end
        % 排除数值上严重病态的设计子矩阵。
        if cond(Xs'*Xs) > 1e8, continue; end
        XtXinv = pinv(Xs'*Xs);
        bsub = XtXinv*Xs'*y;
        if enforcePositive && any(Xcheck(:,keep)*bsub < yFloor), continue; end
        e = y-Xs*bsub;
        h = sum((Xs*XtXinv).*Xs,2);
        if any(h>=.99), continue; end
        PRESS = sum((e./(1-h)).^2);
        % PRESS相同或非常接近时优先保留项更少的模型
        if PRESS < bestPRESS-1e-10 || ...
                (abs(PRESS-bestPRESS)<=1e-10 && p<sum(bestKeep))
            bestPRESS = PRESS; beta=zeros(nTerms,1); beta(keep)=bsub; bestKeep=keep;
        end
    end
    st = struct('selected',any(bestKeep),'keep',bestKeep);
end

function ok = checkHeredity(keep)
    % 主效应索引：2,3,4,5；二次项6,7,8；交互项9:14
    ok = true;
    if keep(6) && ~keep(2), ok=false; return; end
    if keep(7) && ~keep(3), ok=false; return; end
    if keep(8) && ~keep(4), ok=false; return; end
    pairs = [2 3 9; 2 4 10; 3 4 11; 2 5 12; 3 5 13; 4 5 14];
    for i=1:size(pairs,1)
        if keep(pairs(i,3)) && (~keep(pairs(i,1)) || ~keep(pairs(i,2)))
            ok=false; return;
        end
    end
end

function st = diagnoseModel(X,y,beta,keep,termNames)
    Xs = X(:,keep); b = beta(keep); n=size(X,1); p=sum(keep); df=n-p;
    invXtX = pinv(Xs'*Xs); yfit=Xs*b; e=y-yfit;
    h=sum((Xs*invXtX).*Xs,2); sigma=sqrt(sum(e.^2)/max(df,1));
    eLoo=e./max(1-h,eps); PRESS=sum(eLoo.^2); SST=sum((y-mean(y)).^2);
    st.keep=keep; st.beta=beta; st.invXtX=invXtX; st.df=df;
    st.rank=rank(Xs); st.condNumber=cond(Xs'*Xs); st.n=n; st.p=p;
    st.yFit=yfit; st.residual=e; st.leverage=h; st.sigma=sigma; st.s=sigma;
    st.e_loo=eLoo; st.Q2=1-PRESS/max(SST,eps); st.RMSE_loo=sqrt(PRESS/n);
    st.AICc = n*log(max(sum(e.^2)/n,eps)) + 2*p + (2*p*(p+1))/max(n-p-1,1);
    st.pValue=nan(size(beta));
    if df>0
        se=sigma*sqrt(max(diag(invXtX),0)); t=b./max(se,eps);
        st.pValue(keep)=2*tcdf(-abs(t),df);
    end
end

function [yhat,lo,hi] = predictWithPI(Xnew,st,beta)
    % 95%%预测区间：包含新观测的随机误差项，因此比均值置信区间更宽。
    Xs = Xnew(:,st.keep);
    yhat = Xs*beta(st.keep);
    hnew = sum((Xs*st.invXtX).*Xs,2);
    if st.df > 0 && isfinite(st.sigma)
        tcrit = tinv(0.975,st.df);
        sePred = st.sigma*sqrt(max(1+hnew,0));
        lo = yhat-tcrit*sePred;
        hi = yhat+tcrit*sePred;
    else
        lo = nan(size(yhat)); hi = nan(size(yhat));
    end
end

function [GRG,ymin,ymax] = computeGRG(Y,W,zeta)
    % Larger-the-better grey relational grade.
    ymin = min(Y,[],1);
    ymax = max(Y,[],1);
    scale = max(ymax-ymin,eps);
    Z = (Y-ymin)./scale;
    GRC = zeta./(1-Z+zeta);
    GRG = GRC*W(:);
end

function [q2,rmse,yPred] = nestedLOO(X,y,Xcheck,minDf,maxTerms,enforcePositive,enforceHeredity)
    n=size(X,1); yPred=nan(n,1);
    for i=1:n
        train=true(n,1); train(i)=false;
        yTrain=y(train); XTrain=X(train,:);
        yFloor=0;
        % 外层留一时，正值约束只使用外层训练集，避免将被留出的试验点
        % 提前用于模型筛选。
        [b,st]=selectModel(XTrain,yTrain,XTrain,minDf,maxTerms,enforcePositive,enforceHeredity,yFloor);
        if ~st.selected
            [b,st]=selectModel(XTrain,yTrain,XTrain,minDf,maxTerms,false,enforceHeredity,yFloor);
        end
        yPred(i)=X(i,st.keep)*b(st.keep);
    end
    PRESS=sum((y-yPred).^2); SST=sum((y-mean(y)).^2);
    q2=1-PRESS/max(SST,eps); rmse=sqrt(PRESS/n);
end

function validateTable(V,matIdx,buildX,betaUse,statsUse,respNames)
    if isempty(V), return; end
    Xv=buildX(V(:,1),V(:,2),V(:,3),V(:,4));
    fprintf('\nMaterial %d validation:\n',matIdx);
    for j=1:3
        st=statsUse{matIdx,j}; b=betaUse{matIdx,j};
        [yp,lo,hi]=predictWithPI(Xv,st,b);
        fprintf('%s predictions: ',respNames{j}); fprintf('%.4g ',yp); fprintf('\n');
        fprintf('  95%% PI lower: '); fprintf('%.4g ',lo); fprintf('\n');
        fprintf('  95%% PI upper: '); fprintf('%.4g ',hi); fprintf('\n');
        for r=1:size(V,1)
            if ~isnan(V(r,4+j))
                err=abs(V(r,4+j)-yp(r))./max(abs(yp(r)),eps)*100;
                fprintf('  row%d measured=%.4g, relative error=%.2f%%\n', ...
                    r,V(r,4+j),err);
            end
        end
    end
end

function reportValidationGRG(V,D,name,W,zeta)
    if isempty(V), return; end
    ok = all(~isnan(V(:,5:7)),2);
    if ~any(ok)
        fprintf('%s：尚未填入确认试验实测均值，跳过实测GRG。\n',name);
        return;
    end
    Vok = V(ok,:);
    % 确认点不参与模型拟合，仅在最终实测候选比较中追加。
    Dconfirm = [900+(1:size(Vok,1))', Vok(:,1:7)];
    Dall = [D;Dconfirm];
    [g,~,~] = computeGRG(Dall(:,6:8),W,zeta);
    [gs,ix] = sort(g,'descend');
    fprintf('%s：加入确认点后的实测GRG排序：\n',name);
    for k=1:min(5,length(ix))
        p=Dall(ix(k),2:5);
        if ix(k)<=size(D,1)
            idText=sprintf('ID%d',Dall(ix(k),1));
        else
            idText=sprintf('Confirm-%d',ix(k)-size(D,1));
        end
        fprintf('  Rank%d (%s): layer=%.3f mm, density=%.1f%%, angle=%.1f deg, pattern=%s, GRG=%.4f\n', ...
            k,idText,p(1),p(2),p(3),patternLabel(p(4)),gs(k));
    end
end
