function results = fdm_sci_revision(dataFile, cfg)
% FDM flexural analysis for small mixed-factor experimental designs.
% MATLAB R2020b+, base MATLAB (DOCX import requires a JVM).
% Example:
%   cfg = struct('LegacyExcludedIDs',[]);
%   R = fdm_sci_revision('manuscript.docx',cfg);
% CSV: GroupID,CF,Layer,Density,Angle,Pattern,FlexStrength,FlexModulus,
% FlexToughness. Pattern: -1=Triangle, +1=Honeycomb. Optional tensile columns:
% TensileStrength,TensileModulus,TensileToughness; CandidateEligible (0/1).
% Strength/modulus: MPa; integral stress(MPa)*strain(fraction): MJ/m^3.
% This program does NOT verify the integration units from condition means.
% No observations are removed automatically. No synthetic replicates.
% Candidate predictions are hypotheses; no continuous optimum is certified.

if nargin < 1, dataFile = 'fdm_input.csv'; end
if nargin < 2, cfg = struct(); end
cfg = defaults(cfg);
D = readInput(dataFile,cfg);
validateInput(D,cfg);
if ~exist(cfg.OutDir,'dir'), mkdir(cfg.OutDir); end
writetable(D,fullfile(cfg.OutDir,'imported_condition_means.csv'));
results = struct('Config',cfg,'Data',D,'Source',char(dataFile));
results.ResponseUnits = {'MPa','MPa','MJ/m^3 (integration must be checked)'};
results.Scope = ['Condition means; internal CV; conditional model inference; ', ...
    'no specimen-level error estimate or external validation is inferred.'];
terms = {'Intercept','Layer','Density','Angle','Pattern','Layer2', ...
    'Density2','Angle2','LayerDensity','LayerAngle','DensityAngle', ...
    'LayerPattern','DensityPattern','AnglePattern'};
responseNames = {'FlexStrength','FlexModulus','FlexToughness'};
cfLevels = unique(D.CF)';
allDiagnostics = table();
for mi = 1:numel(cfLevels)
    cf = cfLevels(mi); B = D(D.CF==cf,:);
    Yref = B{:,responseNames};
    [Zref,Cref,bounds] = greyCoefficients(Yref,cfg.Zeta,[]);
    [W,Wcritic] = chooseWeights(Zref,cfg);
    R = struct('CF',cf,'Weights',W,'CRITICWeights',Wcritic, ...
        'NormalizationBounds',bounds,'Correlation',corrcoef(Yref), ...
        'all',[],'screened',[],'SuggestedAugmentation',[]);
    weightTable = table();
    weightTable.Objective = responseNames';
    weightTable.UsedWeight = W'; weightTable.CRITICWeight = Wcritic';
    writetable(weightTable,fullfile(cfg.OutDir,sprintf('CF%g_weights.csv',cf)));
    scenarios = {'all'};
    if any(ismember(B.GroupID,cfg.LegacyExcludedIDs)), scenarios{end+1}='screened'; end
    for si = 1:numel(scenarios)
        name = scenarios{si}; keep = true(height(B),1);
        if strcmp(name,'screened'), keep = ~ismember(B.GroupID,cfg.LegacyExcludedIDs); end
        T = B(keep,:); n = height(T); assert(n>=8,'Each fit domain requires >=8 unique conditions.');
        tag = sprintf('CF%g_%s',cf,name); X = features(T,cfg);
        Y = T{:,responseNames};
        eligible = logical(T.CandidateEligible) & ~ismember(T.GroupID,cfg.CandidateExcludedIDs);
        if any(isfinite(cfg.MinTensile))
            A = T{:,{'TensileStrength','TensileModulus','TensileToughness'}};
            required = isfinite(cfg.MinTensile);
            assert(all(all(isfinite(A(:,required)))), ...
                'Tensile constraints require actual tensile results.');
            eligible = eligible & all(A(:,required)>=cfg.MinTensile(required),2);
        end
        [~,C,~] = greyCoefficients(Y,cfg.Zeta,bounds);
        [~,Clocal,~] = greyCoefficients(Y,cfg.Zeta,[]);
        ranking = T;
        ranking.CandidateEligible = eligible;
        ranking.GRG_CommonReference = C*W';
        ranking.GRG_LocalReference = Clocal*W';
        ranking.ParetoEligible = false(n,1);
        ranking.ParetoEligible(eligible) = paretoMax(Y(eligible,:));
        ranking = sortrows(ranking,{'CandidateEligible','GRG_CommonReference'},{'descend','descend'});
        writetable(ranking,fullfile(cfg.OutDir,[tag '_measured_ranking.csv']));
        sensitivity = weightSensitivity(C(eligible,:),T.GroupID(eligible));
        writetable(sensitivity,fullfile(cfg.OutDir,[tag '_weight_sensitivity.csv']));
        design = designAudit(X,terms);
        writetable(design.Summary,fullfile(cfg.OutDir,[tag '_design_rank.csv']));
        writetable(design.AliasBasis,fullfile(cfg.OutDir,[tag '_alias_nullspace.csv']));
        % Fixed budget across inner/outer folds within this domain. minDf is
        % an explicit modeling choice, NOT a universal publication standard.
        budget = min(cfg.MaxParameters,n-2-cfg.MinCVResidualDF);
        assert(budget>=5,'Too few conditions for the predeclared main-effects model.');
        pool = candidatePool(budget,cfg.MaxExtraTerms);
        models = cell(1,3); diagnostics = table(); anovaRows = table();
        for j = 1:3
            y = Y(:,j);
            base = fitModel(X,y,1:5);
            assert(base.Valid,'The main-effects design is not independently estimable.');
            selected = selectModel(X,y,pool);
            nested = nan(n,1); foldKeep = cell(n,1);
            for ii = 1:n
                train = true(n,1); train(ii)=false;
                fold = selectModel(X(train,:),y(train),pool);
                nested(ii)=X(ii,fold.keep)*fold.beta;
                foldKeep{ii}=fold.keep;
            end
            selected.NestedPrediction = nested; selected.FoldTerms = foldKeep;
            selected.MainEffects = base;
            selected.Q2Nested = q2(y,nested);
            selected.RMSENested = sqrt(mean((y-nested).^2));
            baselineLOO = (sum(y)-y)/(n-1);
            selected.Q2VsLOOMean = 1-sum((y-nested).^2)/sum((y-baselineLOO).^2);
            models{j} = selected;
            d = table(); d.CF=cf; d.Scenario={name}; d.Response=responseNames(j);
            d.N=n; d.SelectedParameters=numel(selected.keep); d.ResidualDF=selected.df;
            d.DesignFullRank=design.Rank; d.R2=selected.R2; d.AdjustedR2=selected.AdjustedR2;
            d.Q2FixedSelected=q2(y,selected.loo);
            d.Q2Nested=selected.Q2Nested; d.Q2VsLOOMean=selected.Q2VsLOOMean;
            d.RMSENested=selected.RMSENested; d.Q2MainEffects=q2(y,base.loo);
            d.RMSEMainEffects=sqrt(mean((y-base.loo).^2));
            d.MaxLeverage=max(selected.h); d.ConditionNumber=selected.condX;
            d.Terms={strjoin(terms(selected.keep),' + ')};
            d.InnerCVResidualDFMinimum=n-2-numel(selected.keep);
            diagnostics=appendTable(diagnostics,d);
            coefficients = coefficientTable(selected,terms);
            writetable(coefficients,fullfile(cfg.OutDir,[tag '_' responseNames{j} '_coefficients.csv']));
            pred = table(); pred.GroupID=T.GroupID; pred.Observed=y;
            pred.Fitted=selected.fit; pred.FixedModelLOO=selected.loo;
            pred.NestedLOO=nested; pred.MainEffectsLOO=base.loo;
            pred.Residual=y-selected.fit; pred.Leverage=selected.h;
            pred.InternallyStudentized=selected.studentized;
            pred.CooksDistance=selected.cooks;
            pred.OuterFoldTerms=cellfun(@(k)strjoin(terms(k),'+'),foldKeep,'UniformOutput',false);
            writetable(pred,fullfile(cfg.OutDir,[tag '_' responseNames{j} '_predictions.csv']));
            freq = zeros(14,1);
            for ii=1:n, freq(foldKeep{ii})=freq(foldKeep{ii})+1; end
            tf = table(); tf.Term=terms'; tf.OuterFoldFrequency=freq/n;
            writetable(tf,fullfile(cfg.OutDir,[tag '_' responseNames{j} '_term_stability.csv']));
            % Planned additive quadratic ANOVA: four grouped partial F tests.
            % Error is model residual on condition means, NOT pure error.
            ar = plannedAnova(X,y,responseNames{j});
            anovaRows=appendTable(anovaRows,ar);
        end
        if ~isempty(anovaRows), anovaRows.Q_BH=bhAdjust(anovaRows.P_ModelConditional); end
        writetable(anovaRows,fullfile(cfg.OutDir,[tag '_planned_ANOVA.csv']));
        writetable(diagnostics,fullfile(cfg.OutDir,[tag '_model_diagnostics.csv']));
        allDiagnostics=appendTable(allDiagnostics,diagnostics);
        G = predictionGrid(B,cfg);
        proposals = predictGrid(G,T,B,models,bounds,W,cfg);
        writetable(proposals,fullfile(cfg.OutDir,[tag '_UNCONFIRMED_grid_predictions.csv']));
        external = externalCheck(cfg.ConfirmationFile,cf,T,models,cfg);
        if ~isempty(external)
            writetable(external,fullfile(cfg.OutDir,[tag '_external_confirmation.csv']));
        end
        S = struct('Data',T,'Ranking',ranking,'Design',design,'Models',{models}, ...
            'Diagnostics',diagnostics,'ANOVA',anovaRows,'WeightSensitivity',sensitivity, ...
            'UnconfirmedCandidates',proposals,'ExternalConfirmation',external);
        R.(name)=S;
        if cfg.MakeFigures, drawDiagnostics(T,models,responseNames,tag,cfg.OutDir); end
        fprintf('\n%s: full quadratic rank %d/14; observed eligible Pareto IDs: ',tag,design.Rank);
        fprintf('%g ',ranking.GroupID(ranking.ParetoEligible)); fprintf('\n');
        disp(diagnostics(:,{'Response','Q2MainEffects','Q2Nested','RMSENested','Terms'}));
    end
    augmentation = augmentDesign(B,cfg);
    writetable(augmentation,fullfile(cfg.OutDir,sprintf('CF%g_suggested_design_augmentation.csv',cf)));
    R.SuggestedAugmentation=augmentation;
    if cfg.MakeFigures, drawObserved(B,Zref,Cref,W,cf,cfg.OutDir); end
    results.Material(mi)=R;
end
results.Diagnostics=allDiagnostics;
writetable(allDiagnostics,fullfile(cfg.OutDir,'all_model_diagnostics.csv'));
save(fullfile(cfg.OutDir,'analysis_results.mat'),'results');
writeNotes(cfg.OutDir,cfg);
fprintf('\nSaved to %s\n',cfg.OutDir);
fprintf('A high fitted R2 does not certify prediction. Grid results remain unconfirmed.\n');
end

function c = defaults(c)
d = struct('OutDir',fullfile(pwd,'fdm_sci_results'),'MakeFigures',true, ...
    'DocumentCF',[10 15],'LegacyExcludedIDs',[],'CandidateExcludedIDs',[], ...
    'MinTensile',[-Inf -Inf -Inf],'MaxParameters',7,'MaxExtraTerms',2, ...
    'MinCVResidualDF',3,'WeightMode','equal','Weights',[1 1 1]/3, ...
    'SubjectiveWeights',[],'SubjectiveFraction',0.5,'Zeta',0.5, ...
    'MaxAcceptableCVRMSE',[NaN NaN NaN],'ConfirmationFile','', ...
    'AugmentationCount',5,'Center',[.35 60 45],'Scale',[.05 20 45]);
f = fieldnames(c); known = fieldnames(d);
assert(all(ismember(f,known)),'Unknown configuration field; check spelling.');
for i=1:numel(known)
    if ~isfield(c,known{i}), c.(known{i})=d.(known{i}); end
end
assert(c.Zeta>0 && c.Zeta<=1,'Zeta must be in (0,1].');
assert(numel(c.Center)==3 && numel(c.Scale)==3 && all(c.Scale>0),'Invalid coding.');
assert(numel(c.MinTensile)==3 && ~any(isnan(c.MinTensile)),'Invalid tensile constraints.');
assert(numel(c.MaxAcceptableCVRMSE)==3,'Supply three engineering RMSE limits.');
assert(c.MaxParameters>=5 && c.MaxParameters<=14 && c.MaxParameters==fix(c.MaxParameters),'Invalid parameter budget.');
assert(ismember(c.MaxExtraTerms,0:2),'MaxExtraTerms must be 0, 1 or 2.');
assert(c.MinCVResidualDF>=1 && c.MinCVResidualDF==fix(c.MinCVResidualDF),'Invalid minimum CV residual DF.');
assert(c.AugmentationCount>=0 && c.AugmentationCount==fix(c.AugmentationCount),'Invalid augmentation count.');
end

function D = readInput(file,c)
[~,~,ext] = fileparts(file);
if strcmpi(ext,'.docx')
    D = readDocxMeans(file,c.DocumentCF);
else
    D = readtable(file,'TextType','string');
end
required = {'GroupID','CF','Layer','Density','Angle','Pattern','FlexStrength','FlexModulus','FlexToughness'};
assert(all(ismember(required,D.Properties.VariableNames)), ...
    'Missing required columns. See the header of fdm_sci_revision.m.');
if ~isnumeric(D.Pattern)
    p = lower(strtrim(string(D.Pattern))); code=nan(height(D),1);
    code(ismember(p,["triangle","triangular","-1"]))=-1;
    code(ismember(p,["honeycomb","1","+1"]))=1; D.Pattern=code;
end
optional = {'TensileStrength','TensileModulus','TensileToughness'};
for i=1:3
    if ~ismember(optional{i},D.Properties.VariableNames), D.(optional{i})=nan(height(D),1); end
end
if ~ismember('CandidateEligible',D.Properties.VariableNames), D.CandidateEligible=true(height(D),1); end
end

function D = readDocxMeans(file,levels)
assert(usejava('jvm'),'DOCX import requires a JVM; alternatively supply a CSV.');
tmp = tempname; mkdir(tmp); cleanup=onCleanup(@()rmdir(tmp,'s')); %#ok<NASGU>
unzip(file,tmp); doc=xmlread(fullfile(tmp,'word','document.xml'));
tables=doc.getElementsByTagName('w:tbl'); data=[]; tableCount=0;
for it=0:tables.getLength-1
    rows=tables.item(it).getElementsByTagName('w:tr'); block=[];
    for ir=0:rows.getLength-1
        cells=rows.item(ir).getElementsByTagName('w:tc');
        if cells.getLength~=11, continue; end
        v=cell(1,11);
        for ic=0:10
            tt=cells.item(ic).getElementsByTagName('w:t'); text='';
            for k=0:tt.getLength-1, text=[text char(tt.item(k).getTextContent)]; end %#ok<AGROW>
            v{ic+1}=strtrim(text);
        end
        id=str2double(v{1}); pat=lower(v{5});
        if ~isfinite(id) || ~ismember(string(pat),["triangle","honeycomb"]), continue; end
        p=-1; if strcmp(pat,'honeycomb'), p=1; end
        nums=cellfun(@str2double,v([1 2 3 4 6 7 8 9 10 11]));
        assert(all(isfinite(nums)),'Non-numeric value in a recognized experimental table.');
        block=[block;nums(1:4),p,nums(5:10)]; %#ok<AGROW>
    end
    if isempty(block), continue; end
    tableCount=tableCount+1;
    assert(tableCount<=numel(levels),'More matching tables than DocumentCF entries. Use CSV for unambiguous mapping.');
    data=[data;block(:,1),repmat(levels(tableCount),size(block,1),1),block(:,2:end)]; %#ok<AGROW>
end
assert(tableCount==numel(levels),'Matching table count differs from DocumentCF; verify the two table layouts or use CSV.');
D=array2table(data,'VariableNames',{'GroupID','CF','Layer','Density','Angle','Pattern', ...
    'TensileStrength','TensileModulus','TensileToughness','FlexStrength','FlexModulus','FlexToughness'});
fprintf('Imported %d condition means from %d DOCX tables; CF mapping follows DocumentCF.\n',height(D),tableCount);
end

function validateInput(D,c)
vars={'GroupID','CF','Layer','Density','Angle','Pattern','FlexStrength','FlexModulus','FlexToughness'};
for k=1:numel(vars), assert(isnumeric(D.(vars{k})),['Expected numeric column: ' vars{k}]); end
A=D{:,vars}; assert(all(isfinite(A(:))),'Required inputs must be finite.');
assert(numel(unique(D.GroupID))==height(D),'GroupID must be globally unique; use global IDs across materials.');
assert(all(ismember(D.Pattern,[-1 1])),'Pattern must be Triangle(-1) or Honeycomb(+1).');
assert(all(all(D{:,{'FlexStrength','FlexModulus','FlexToughness'}}>0)),'Flexural means must be positive.');
assert(all(ismember(D.CandidateEligible,[0 1])),'CandidateEligible must be 0/1.');
assert(all(all(abs((D{:,{'Layer','Density','Angle'}}-c.Center)./c.Scale)<=1+1e-8)), ...
    'Data exceed configured factor bounds. Change Center/Scale explicitly for a different design.');
assert(size(unique(D{:,{'CF','Layer','Density','Angle','Pattern'}},'rows'),1)==height(D), ...
    'Duplicate design conditions: supply one mean per condition; keep specimen records separately.');
end

function X=features(T,c)
u=(T{:,{'Layer','Density','Angle'}}-c.Center)./c.Scale;
a=u(:,1); b=u(:,2); d=u(:,3); p=T.Pattern;
X=[ones(height(T),1),u,p,a.^2,b.^2,d.^2,a.*b,a.*d,b.*d,a.*p,b.*p,d.*p];
end

function A=designAudit(X,terms)
A.Rank=rank(X); [~,s,~]=svd(X,0); sv=diag(s);
t=table(); t.N=size(X,1); t.Columns=14; t.Rank=A.Rank;
t.FullQuadraticUniqueCoefficients=(A.Rank==14); t.MainEffectsRank=rank(X(:,1:5));
t.AdditiveQuadraticRank=rank(X(:,1:8));
t.SmallestSingularValue=sv(end); A.Summary=t;
N=null(X); z=table(); z.Term=terms';
for i=1:size(N,2), z.(sprintf('AliasBasis%d',i))=N(:,i); end
A.AliasBasis=z; A.NullSpace=N;
end

function pool=candidatePool(budget,maxExtra)
pool={1,1:5}; extras=6:14;
for count=1:maxExtra
    if 5+count>budget, continue; end
    combinations=nchoosek(extras,count);
    for i=1:size(combinations,1), pool{end+1}=[1:5 combinations(i,:)]; end %#ok<AGROW>
end
% All nonlinear terms retain their main-effect parents (strong heredity).
end

function m=fitModel(X,y,keep)
A=X(:,keep); [Q,R]=qr(A,0); p=numel(keep); n=numel(y);
m=struct('Valid',false,'keep',keep);
if rank(A)<p || n<=p || cond(R)>1e8, return; end
m.Valid=true; m.beta=R\(Q'*y); m.fit=A*m.beta; m.residual=y-m.fit;
m.h=sum(Q.^2,2); m.df=n-p; m.SSE=sum(m.residual.^2); m.MSE=m.SSE/m.df;
m.condX=cond(R); m.R=R;
m.loo=nan(n,1);
if all(1-m.h>1e-9), m.loo=y-m.residual./(1-m.h); end
m.PRESS=sum((y-m.loo).^2); if any(~isfinite(m.loo)), m.PRESS=Inf; end
m.R2=q2(y,m.fit); m.AdjustedR2=1-(1-m.R2)*(n-1)/m.df;
m.studentized=m.residual./sqrt(max(m.MSE,realmin)*max(1-m.h,realmin));
m.cooks=(m.studentized.^2/p).*m.h./max(1-m.h,realmin);
end

function best=selectModel(X,y,pool)
best=struct('Valid',false,'PRESS',Inf,'keep',1:100);
for i=1:numel(pool)
    m=fitModel(X,y,pool{i}); if ~m.Valid || ~isfinite(m.PRESS), continue; end
    better=~best.Valid || m.PRESS<best.PRESS-1e-10*max(1,best.PRESS);
    tied=best.Valid && abs(m.PRESS-best.PRESS)<=1e-10*max(1,best.PRESS);
    if better || (tied && numel(m.keep)<numel(best.keep)), best=m; end
end
assert(best.Valid,'No estimable candidate model with valid inner leave-one-out predictions.');
end

function q=q2(y,p)
sst=sum((y-mean(y)).^2);
if sst<=eps(max(1,sum(y.^2))), q=NaN; else, q=1-sum((y-p).^2)/sst; end
end

function T=coefficientTable(m,terms)
V=m.R\eye(numel(m.keep)); se=sqrt(max(m.MSE,0)*sum(V.^2,2));
t=m.beta./max(se,realmin); p=betainc(m.df./(m.df+t.^2),m.df/2,.5);
T=table(); T.Term=terms(m.keep)'; T.Estimate=m.beta; T.ConditionalSE=se;
T.ExploratoryP_PostSelection=p;
end

function T=plannedAnova(X,y,response)
% Additive curvature model is predeclared; interactions are omitted assumptions.
keep=1:8; full=fitModel(X,y,keep); T=table();
if ~full.Valid || full.df<3, return; end
groups={[2 6],[3 7],[4 8],5}; names={'Layer_linear_and_quadratic','Density_linear_and_quadratic', ...
    'Angle_linear_and_quadratic','Pattern'};
for i=1:4
    reduced=fitModel(X,y,setdiff(keep,groups{i}));
    if ~reduced.Valid, continue; end
    df=numel(groups{i}); ss=max(0,reduced.SSE-full.SSE);
    if full.MSE>0, F=(ss/df)/full.MSE; elseif ss>0, F=Inf; else, F=NaN; end
    p=betainc(full.df/(full.df+df*F),full.df/2,df/2);
    a=table(); a.Response={response}; a.Factor=names(i); a.DF=df;
    a.PartialSS=ss; a.F=F; a.P_ModelConditional=p;
    a.ResidualDF=full.df; a.ResidualSS=full.SSE;
    T=appendTable(T,a);
end
end

function q=bhAdjust(p)
[v,ix]=sort(p); n=numel(v); a=v.*n./(1:n)';
for i=n-1:-1:1, a(i)=min(a(i),a(i+1)); end
q=zeros(size(p)); q(ix)=min(1,a);
end

function [Z,C,bounds]=greyCoefficients(Y,zeta,bounds)
if isempty(bounds), bounds=[min(Y,[],1);max(Y,[],1)]; end
span=bounds(2,:)-bounds(1,:); spanSafe=span; spanSafe(span==0)=1;
Z=(Y-bounds(1,:))./spanSafe; Z(:,span==0)=1;
% Fixed ideal and range. Out-of-reference predictions are flagged upstream.
C=zeta./(abs(1-Z)+zeta);
end

function [W,wc]=chooseWeights(Z,c)
sd=std(Z,0,1); R=corrcoef(Z); R(~isfinite(R))=0; R(1:4:end)=1;
information=sd.*sum(1-R,2)'; information=max(information,0);
if sum(information)<=1e-12, wc=[1 1 1]/3; else, wc=information/sum(information); end
switch lower(c.WeightMode)
    case 'equal', W=[1 1 1]/3;
    case 'fixed', W=c.Weights;
    case 'critic', W=wc;
    case 'blend'
        a=c.SubjectiveFraction;
        assert(a>=0 && a<=1,'SubjectiveFraction must be in [0,1].');
        ws=c.SubjectiveWeights;
        assert(numel(ws)==3 && all(ws>=0) && abs(sum(ws)-1)<1e-8, ...
            'Supply documented subjective weights; BWM preferences are not inferred.');
        W=a*ws(:)'+(1-a)*wc;
    otherwise, error('WeightMode: equal, fixed, critic or blend.');
end
assert(numel(W)==3 && all(isfinite(W)) && all(W>=0) && abs(sum(W)-1)<1e-8,'Invalid weights.');
W=W(:)';
end

function p=paretoMax(Y)
p=true(size(Y,1),1);
for i=1:size(Y,1), p(i)=~any(all(Y>=Y(i,:),2)&any(Y>Y(i,:),2)); end
end

function T=weightSensitivity(C,ids)
T=table(); if isempty(ids), return; end
counts=zeros(numel(ids),1); experiments=0;
for a=0:20
    for b=0:20-a
        w=[a b 20-a-b]/20; s=C*w'; win=abs(s-max(s))<=1e-12;
        counts=counts+double(win); experiments=experiments+1;
    end
end
T.GroupID=ids; T.TiedOrSoleWinnerCount=counts;
T.WeightGridPoints=repmat(experiments,numel(ids),1);
T.WinnerFraction=counts/experiments; % Ties counted for every tied condition.
end

function G=predictionGrid(B,c)
lev=cell(1,3); names={'Layer','Density','Angle'};
for k=1:3, lev{k}=unique(B.(names{k})); end
[a,b,d,p]=ndgrid(lev{1},lev{2},lev{3},[-1 1]);
G=table(); G.CF=repmat(B.CF(1),numel(a),1); G.Layer=a(:);
G.Density=b(:); G.Angle=d(:); G.Pattern=p(:);
% Discrete tested factor levels. No rounding of a categorical continuous gene.
end

function P=predictGrid(G,T,allData,models,bounds,W,c)
P=G; Xg=features(G,c); U=features(T,c); n=height(G);
P.ExactTrainingCondition=ismember(Xg(:,2:5),U(:,2:5),'rows');
Xobserved=features(allData,c);
P.PreviouslyMeasured=ismember(Xg(:,2:5),Xobserved(:,2:5),'rows');
P.PredeclaredEligible=true(n,1);
for i=1:n
    k=find(all(abs(Xobserved(:,2:5)-Xg(i,2:5))<1e-10,2),1);
    if ~isempty(k)
        P.PredeclaredEligible(i)=logical(allData.CandidateEligible(k)) && ...
            ~ismember(allData.GroupID(k),c.CandidateExcludedIDs);
    end
end
P.WithinPatternHull=false(n,1);
for pat=[-1 1]
    train=U(T.Pattern==pat,2:4); ix=G.Pattern==pat; query=Xg(ix,2:4);
    inside=ismember(query,train,'rows');
    if size(train,1)>=4 && rank(train-mean(train,1))==3
        try
            tri=delaunayn(train); inside=inside|~isnan(tsearchn(train,tri,query));
        catch
            % Degenerate numerical hull: retain only exact observed support.
        end
    end
    P.WithinPatternHull(ix)=inside;
end
Yhat=zeros(n,3); pass=true; names={'FlexStrength','FlexModulus','FlexToughness'};
for j=1:3
    m=models{j}; A=Xg(:,m.keep); Yhat(:,j)=A*m.beta;
    h=sum((A/m.R).^2,2); z=betaincinv(.05,m.df/2,.5); tc=sqrt(m.df*(1/z-1));
    half=tc*sqrt(m.MSE*(1+h));
    P.(['Predicted_' names{j}])=Yhat(:,j);
    P.(['ConditionalPI95Lower_' names{j}])=Yhat(:,j)-half;
    P.(['ConditionalPI95Upper_' names{j}])=Yhat(:,j)+half;
    limit=c.MaxAcceptableCVRMSE(j);
    pass=pass && isfinite(limit) && limit>0 && m.RMSENested<=limit && m.Q2VsLOOMean>0;
end
P.PositivePredictions=all(Yhat>0,2);
P.OutsideMeasuredResponseRange=any(Yhat<bounds(1,:) | Yhat>bounds(2,:),2);
[~,C,~]=greyCoefficients(Yhat,c.Zeta,bounds);
P.GRG_ModelReference=C*W';
P.InternalEngineeringGate=repmat(pass,n,1);
P.TensileFeasibilityVerified=false(n,1);
if ~any(isfinite(c.MinTensile)), P.TensileFeasibilityVerified=true(n,1); end
% Actual tensile results may certify the constraints only at measured points.
if any(isfinite(c.MinTensile))
    required=isfinite(c.MinTensile);
    V=allData{:,{'TensileStrength','TensileModulus','TensileToughness'}};
    for i=1:n
        k=find(all(abs(Xobserved(:,2:5)-Xg(i,2:5))<1e-10,2),1);
        if ~isempty(k), P.TensileFeasibilityVerified(i)=all(isfinite(V(k,required))) && all(V(k,required)>=c.MinTensile(required)); end
    end
end
P.EligibleForProvisionalRanking=P.WithinPatternHull & P.PositivePredictions & P.PredeclaredEligible & ...
    ~P.OutsideMeasuredResponseRange & P.InternalEngineeringGate & P.TensileFeasibilityVerified;
P.ParetoPredictedSupported=false(n,1); e=P.EligibleForProvisionalRanking;
P.ParetoPredictedSupported(e)=paretoMax(Yhat(e,:));
P.RequiresIndependentConfirmation=true(n,1);
end

function E=externalCheck(file,cf,train,models,c)
E=table(); if isempty(file), return; end
V=readInput(file,c); validateInput(V,c); V=V(V.CF==cf,:);
X=features(V,c); A=features(train,c); names={'FlexStrength','FlexModulus','FlexToughness'};
for i=1:height(V)
    for j=1:3
        m=models{j}; x=X(i,m.keep); h=sum((x/m.R).^2,2);
        z=betaincinv(.05,m.df/2,.5); tc=sqrt(m.df*(1/z-1));
        pred=x*m.beta; half=tc*sqrt(m.MSE*(1+h)); measured=V.(names{j})(i);
        a=table(); a.ValidationIndex=i; a.GroupID=V.GroupID(i); a.Response=names(j);
        a.SeenDesign=any(all(abs(A(:,2:5)-X(i,2:5))<1e-10,2));
        a.Measured=measured; a.Predicted=pred;
        a.RelativeErrorPercent=100*(pred-measured)/measured;
        a.ConditionalPI95Lower=pred-half; a.ConditionalPI95Upper=pred+half;
        a.InConditionalPI95=(measured>=pred-half && measured<=pred+half);
        E=appendTable(E,a);
    end
end
% All columns assigned by name; no variable named Row/Variables/Properties.
end

function T=augmentDesign(B,c)
G=predictionGrid(B,c); X=features(B,c); Q=features(G,c);
available=~ismember(Q(:,2:5),X(:,2:5),'rows'); T=G([],:);
for step=1:c.AugmentationCount
    inds=find(available); if isempty(inds), break; end
    ranks=zeros(numel(inds),1); score=ranks;
    for k=1:numel(inds)
        s=svd([X;Q(inds(k),:)]); r=rank([X;Q(inds(k),:)]);
        ranks(k)=r; score(k)=2*sum(log(s(1:r)));
    end
    bestRank=max(ranks); score(ranks<bestRank)=-Inf; [~,k]=max(score); ix=inds(k);
    a=G(ix,:); a.SequentialOrder=step; a.RankAfterAddition=bestRank;
    a.ModelResidualDFIfMeasured=size(X,1)+1-bestRank;
    T=appendTable(T,a); X=[X;Q(ix,:)]; available(ix)=false; %#ok<AGROW>
end
end

function drawDiagnostics(T,models,names,tag,outDir)
f=figure('Visible','off','Color','w','Position',[50 50 1250 1100]);
for j=1:3
    m=models{j}; y=T.(names{j});
    subplot(4,3,j); scatter(y,m.NestedPrediction,30,'filled'); hold on;
    b=[min([y;m.NestedPrediction]),max([y;m.NestedPrediction])]; plot(b,b,'k--');
    xlabel('Measured condition mean'); ylabel('Nested LOO prediction'); title(names{j},'Interpreter','none'); grid on;
    subplot(4,3,j+3); scatter(m.fit,m.residual,30,'filled'); yline(0,'k--');
    xlabel('Fitted mean'); ylabel('Residual'); grid on;
    subplot(4,3,j+6); stem(T.GroupID,m.cooks,'filled');
    xlabel('Group ID (not randomized run order)'); ylabel('Cook distance'); grid on;
    subplot(4,3,j+9); r=sort(m.studentized); n=numel(r);
    theoretical=-sqrt(2)*erfcinv(2*((1:n)'-.5)/n);
    plot(theoretical,r,'o'); xlabel('Standard normal quantile');
    ylabel('Studentized residual quantile'); grid on;
end
print(f,fullfile(outDir,[tag '_diagnostics.png']),'-dpng','-r180'); close(f);
end

function drawObserved(B,Z,C,W,cf,outDir)
f=figure('Visible','off','Color','w','Position',[50 50 1150 420]);
subplot(1,2,1); imagesc(Z); colorbar; caxis([0 1]);
xticks(1:3); xticklabels({'Strength','Modulus','Energy density'});
yticks(1:height(B)); yticklabels(string(B.GroupID)); ylabel('Group ID'); title('Normalized measured means');
subplot(1,2,2); y=B{:,{'FlexStrength','FlexModulus','FlexToughness'}};
scatter3(y(:,1),y(:,2),y(:,3),55,C*W','filled'); hold on;
p=paretoMax(y); scatter3(y(p,1),y(p,2),y(p,3),140,'ro');
xlabel('Strength (MPa)'); ylabel('Modulus (MPa)'); zlabel('Energy density (MJ/m^3)');
title('Measured points; circles = Pareto points'); colorbar; grid on;
print(f,fullfile(outDir,sprintf('CF%g_measured.png',cf)),'-dpng','-r180'); close(f);
end

function T=appendTable(T,a)
if isempty(T), T=a; else, T=[T;a]; end
end

function writeNotes(outDir,c)
f=fopen(fullfile(outDir,'READ_RESULTS.txt'),'w'); assert(f>=0,'Cannot write result notes.');
cleanup=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'Primary outputs: *_all_measured_ranking.csv and all_model_diagnostics.csv.\n');
fprintf(f,'Screened outputs describe a restricted domain, not proof of invalid excluded observations.\n');
fprintf(f,'CommonReference holds within-material normalization and weights fixed across screening scenarios.\n');
fprintf(f,'LocalReference reproduces scenario-specific normalization; scores are not comparable across scenarios/materials.\n');
fprintf(f,'Nested LOO repeats all response-based model selection. Fixed-model Q2 is conditional and may be optimistic after selection.\n');
fprintf(f,'MainEffects Q2 is for a predeclared linear model. Q2VsLOOMean compares with honest training-mean predictions.\n');
fprintf(f,'The ANOVA uses a predeclared additive quadratic model on condition means. Residual SS includes model lack of fit.\n');
fprintf(f,'Partial SS are adjusted tests; they need not sum to total SS. BH q values cover the 12 tests per material/scenario.\n');
fprintf(f,'Coefficient p values and prediction intervals condition on the selected model; neither includes model-selection uncertainty.\n');
fprintf(f,'Prediction intervals concern comparable condition means under model assumptions, not individual specimens.\n');
fprintf(f,'No pure-error lack-of-fit test, specimen SD, robustness S/N or bootstrap significance is fabricated.\n');
fprintf(f,'No mass/density normalization, fracture toughness or structural EI is inferred from flexural modulus.\n');
fprintf(f,'Prediction grid is discrete and exploratory. Hull membership is geometric support, not a prediction guarantee.\n');
fprintf(f,'New designs with active tensile constraints need tensile measurements; flexural models cannot certify them.\n');
fprintf(f,'Engineering RMSE limits are [%g %g %g]; NaN leaves candidate recommendations disabled.\n',c.MaxAcceptableCVRMSE);
fprintf(f,'Augmentation is sequential rank/D-optimal design advice, NOT an optimum-performance prediction.\n');
fprintf(f,'DOCX import reads numeric table text only; it executes no embedded scripts or instructions. Verify imported CF mapping.\n');
fprintf(f,'Integral stress(MPa) d strain(fraction) has units MJ/m^3; if strain is percent, divide by 100 before integration.\n');
fprintf(f,'If stopping at 5%% strain or first cracking, name the energy endpoint explicitly; it is not fracture toughness.\n');
end
