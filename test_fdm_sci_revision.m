function test_fdm_sci_revision()
% Native MATLAB smoke test using synthetic data only (no research data).
% Run after putting revised1031.m on the MATLAB path.
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
R=revised1031(csv,cfg);
assert(numel(R.Material)==2);
assert(height(R.Diagnostics)==9); % CF10 all+screened; CF15 all only.
assert(height(R.Summary)==3);
assert(strcmp(R.Version,'2026-11-02-total'));
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
if usejava('jvm')
    docx=syntheticDocx(root,T,U);
    docCfg=cfg; docCfg.OutDir=fullfile(root,'docx_results'); docCfg.MaxExtraTerms=0;
    docResult=revised1031(docx,docCfg);
    assert(height(docResult.Data)==height(D));
    assert(isequal(docResult.Data.GroupID,D.GroupID));
    assert(max(abs(docResult.Data.FlexStrength-D.FlexStrength))<1e-12);
    assert(height(docResult.Summary)==3);
end
fprintf('PASS: synthetic CSV/DOCX, table outputs, nested CV, differing scenarios, summary and confirmation.\n');
end

function T=responses(T)
u=(T.Layer-.35)/.05; v=(T.Density-60)/20;
w=(T.Angle-45)/45; p=T.Pattern;
T.FlexStrength=45-4*u+3*v-5*w+2*p+1.5*u.*w;
T.FlexModulus=1500-100*u+70*v-120*w+30*p+50*u.*w;
T.FlexToughness=1.5-.10*u+.08*v-.15*w+.04*p+.02*u.*w;
end

function file=syntheticDocx(root,T,U)
folder=fullfile(root,'docx_package'); mkdir(folder); mkdir(fullfile(folder,'word'));
fileXML=fullfile(folder,'word','document.xml'); fid=fopen(fileXML,'w');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>');
blocks={T,U};
for b=1:2
    A=blocks{b}; fprintf(fid,'<w:tbl>');
    for i=1:height(A)
        pattern='Triangle'; if A.Pattern(i)==1, pattern='Honeycomb'; end
        numbers=[A.GroupID(i),A.Layer(i),A.Density(i),A.Angle(i), ...
            30,1800,2,A.FlexStrength(i),A.FlexModulus(i),A.FlexToughness(i)];
        text=arrayfun(@(x)sprintf('%.17g',x),numbers,'UniformOutput',false);
        text=[text(1:4),{pattern},text(5:end)]; fprintf(fid,'<w:tr>');
        for j=1:11, fprintf(fid,'<w:tc><w:p><w:r><w:t>%s</w:t></w:r></w:p></w:tc>',text{j}); end
        fprintf(fid,'</w:tr>');
    end
    fprintf(fid,'</w:tbl>');
end
fprintf(fid,'</w:body></w:document>'); clear cleanup;
archive=fullfile(root,'synthetic.docx.zip'); zip(archive,{'word/document.xml'},folder);
file=fullfile(root,'synthetic.docx'); movefile(archive,file);
end
