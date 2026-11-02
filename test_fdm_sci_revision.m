function test_fdm_sci_revision()
% Native MATLAB smoke test using synthetic data only (no research data).
% Run after putting fdm_sci_revision.m on the MATLAB path.
root=tempname; mkdir(root); cleanup=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
[a,b,c,p]=ndgrid([.3 .35 .4],[40 60 80],[0 45 90],[-1 1]);
T=table(); T.GroupID=(1001:1000+numel(a))'; T.CF=repmat(10,numel(a),1);
T.Layer=a(:); T.Density=b(:); T.Angle=c(:); T.Pattern=p(:);
T=responses(T); U=T; U.GroupID=U.GroupID+1000; U.CF(:)=15;
D=[T;U]; csv=fullfile(root,'synthetic_means.csv'); writetable(D,csv);
V=T(1,:); V.GroupID=9001; V.Layer=.325; V.Density=50;
V.Angle=22.5; V.Pattern=1; V=responses(V);
V2=V; V2.GroupID=9002; V2.CF=15; V=[V;V2];
confirm=fullfile(root,'synthetic_confirmation.csv'); writetable(V,confirm);
cfg=struct('OutDir',fullfile(root,'results'),'MakeFigures',false, ...
    'LegacyExcludedIDs',[1001 1002],'CandidateExcludedIDs',1003, ...
    'ConfirmationFile',confirm);
R=fdm_sci_revision(csv,cfg);
assert(numel(R.Material)==2);
assert(height(R.Diagnostics)==9); % CF10 all+screened; CF15 all only.
assert(R.Material(1).all.Design.Rank==14);
assert(isempty(R.Material(2).screened));
S=R.Material(1).all;
assert(~S.Ranking.CandidateEligible(S.Ranking.GroupID==1003));
assert(~any(S.UnconfirmedCandidates.InternalEngineeringGate));
assert(height(S.ExternalConfirmation)==3);
assert(ismember('ValidationIndex',S.ExternalConfirmation.Properties.VariableNames));
assert(~ismember('Row',S.ExternalConfirmation.Properties.VariableNames));
assert(~any(S.ExternalConfirmation.SeenDesign));
assert(max(abs(S.ExternalConfirmation.RelativeErrorPercent))<1e-6);
assert(all(isfinite(R.Diagnostics.Q2Nested)));
fprintf('PASS: synthetic CSV, table outputs, nested CV, differing scenario fields, and confirmation table.\n');
end

function T=responses(T)
u=(T.Layer-.35)/.05; v=(T.Density-60)/20;
w=(T.Angle-45)/45; p=T.Pattern;
T.FlexStrength=45-4*u+3*v-5*w+2*p+1.5*u.*w;
T.FlexModulus=1500-100*u+70*v-120*w+30*p+50*u.*w;
T.FlexToughness=1.5-.10*u+.08*v-.15*w+.04*p+.02*u.*w;
end
