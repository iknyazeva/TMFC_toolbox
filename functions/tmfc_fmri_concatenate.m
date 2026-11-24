function SPM = tmfc_fmri_concatenate(SPM,scans)

% This is an in-memory adaptation of the original SPM25 function
% spm_fmri_concatenate. It accepts and returns an SPM structure directly,
% without loading, backing up, or saving SPM.mat files.
%
% FORMAT SPM = tmfc_fmri_concatenate(SPM,scans)
% SPM   - SPM structure with concatenated sessions
% scans - [1 x n] vector with the original number of scans in each session
%
% Original documentation:
%
% Adjust an SPM.mat which has concatenated sessions
% FORMAT spm_fmri_concatenate(P, scans)
% Session regressors are added and the high-pass filter and non-sphericity
% estimates adjusted as if sessions are separate.
%
% P     - filename of the SPM.mat file to adjust
% scans - [1 x n] vector with the original number of scans in each session
%
% The expected workflow is:
%
% 1. Manually specify a GLM with timeseries and onsets concatenated
% 2. Run spm_post_concatenate on the saved SPM.mat.
% 3. Estimate the SPM.mat in the normal way.
%
% Tips:
%
% - The BOLD-response may overhang from one session to the next. To reduce
%   this, acquire additional volumes at the end of each session and / or
%   add regressors to model the trials at the session borders.
%__________________________________________________________________________
%
% Guillaume Flandin & Peter Zeidman
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging
%
%
% TMFC modifications:
% Copyright (C) 2026 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if ~isfield(SPM,'Sess') || numel(SPM.nscan) ~= 1
    error('Input is not a single-session SPM structure.');
end
if sum(scans) ~= SPM.nscan
    error('Number of scans does not match.');
end
SPM.nscan = scans;
if numel(SPM.xX.iB) == 1 && SPM.xX.iB == size(SPM.xX.X,2)
    Xb = [];
    Bn = {};
    for s=1:numel(SPM.nscan)
        Xb = blkdiag(Xb,ones(SPM.nscan(s),1));
        Bn{s} = sprintf('Sn(%i) constant',s);
    end
    SPM.xX.X = [SPM.xX.X(:,1:end-1) Xb];
    SPM.xX.iB = SPM.xX.iB:(SPM.xX.iB+size(Xb,2)-1);
    SPM.xX.name = {SPM.xX.name{1:end-1} Bn{:}};
end
s = cumsum([0 SPM.nscan]);
for i=1:numel(SPM.nscan)
    K(i) = struct('HParam',SPM.xX.K(1).HParam,...
                  'row',s(i)+(1:SPM.nscan(i)),...
                  'RT',SPM.xY.RT);
end
SPM.xX.K = spm_filter(K);
switch lower(SPM.xVi.form)
    case {'ar(1)','ar(0.2)'}
        SPM.xVi.Vi = spm_Ce('ar',SPM.nscan,0.2);
        SPM.xVi.form = 'AR(0.2)';
    case 'fast'
        SPM.xVi.Vi = spm_Ce('fast',SPM.nscan,SPM.xY.RT);
    case {'i.i.d','none'}
    otherwise
        warning('Unhandled temporal non-sphericity.');
end
end
