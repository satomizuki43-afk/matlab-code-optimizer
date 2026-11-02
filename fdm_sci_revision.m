function results = fdm_sci_revision(dataFile, cfg)
% Compatibility entry. Main analysis and all helpers are in revised1031.m.
% This entry preserves the earlier generic defaults (equal weights, no
% historical exclusions, hidden figures). For the study defaults run:
%   results = revised1031;
if nargin<1 || isempty(dataFile), dataFile='fdm_input.csv'; end
if nargin<2 || isempty(cfg), cfg=struct(); end
if ~isfield(cfg,'LegacyExcludedIDs'), cfg.LegacyExcludedIDs=[]; end
if ~isfield(cfg,'WeightMode'), cfg.WeightMode='equal'; end
if ~isfield(cfg,'Weights'), cfg.Weights=[1 1 1]/3; end
if ~isfield(cfg,'ShowFigures'), cfg.ShowFigures=false; end
results=revised1031(dataFile,cfg);
end
