function SPM = tmfc_LSS_design(job,concat_scans)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Build the design, filter and covariance components of one LSS GLM in
% memory, from the same fmri_spec job that tmfc_LSS of main passes to
% spm_jobman (no SPM.mat, no directory, no data access). Used by the
% accelerated tmfc_LSS (branch accelerated_lss) for every trial model.
%
% FORMAT SPM = tmfc_LSS_design(job)
% FORMAT SPM = tmfc_LSS_design(job,concat_scans)
%
%   job           - matlabbatch{1}.spm.stats.fmri_spec of one LSS trial
%   concat_scans  - number of scans in each original session if the parent
%                   SPM.mat has concatenated sessions (SPM.nscan), else []
%
% The job -> SPM conversion mirrors spm_run_fmri_spec.m:54-393 for the job
% fields that tmfc_LSS can produce (HRF without derivatives, no time or
% parametric modulation, no multiple-condition or multiple-regressor files,
% no factors). The design matrix is made by spm_fMRI_design(SPM,0) (save_SPM = 0:
% spm_fmri_spm_ui would write an SPM.mat to the current directory), the filter
% and covariance components mirror spm_fmri_spm_ui.m:198-253. The concatenation
% adjustment mirrors spm_fmri_concatenate.m:60-113. Compared to main, only the
% in-memory specification differs: the resulting design is identical.
%
% =========================================================================
% Derived from tmfc_LSS.m (Copyright (C) 2026 Ruslan Masharipov) and SPM25
% License: GPL-3.0-or-later

if nargin < 2, concat_scans = []; end

if numel(job.sess) ~= 1 || ~isempty(job.fact) || ~isequal(job.bases.hrf.derivs,[0 0]) ...
        || ~isempty(job.sess.multi{1}) || ~isempty(job.sess.multi_reg{1})
    error('tmfc_LSS_design: job contains fields that are not produced by tmfc_LSS.');
end
sess = job.sess;

% Timing parameters and basis functions (spm_run_fmri_spec.m:59-103)
SPM.xY.RT     = job.timing.RT;
SPM.xBF.UNITS = job.timing.units;
SPM.xBF.T     = job.timing.fmri_t;
SPM.xBF.T0    = job.timing.fmri_t0;
SPM.xBF.name  = 'hrf';
SPM.xBF.Volterra = job.volt;

SPM.nscan(1) = numel(sess.scans);
if SPM.nscan(1) == 1
    error('Not enough scans in session 1.');
end

% Conditions (spm_run_fmri_spec.m:225-287)
U = [];
for j = 1:numel(sess.cond)
    cond      = sess.cond(j);
    if cond.tmod > 0 || ~isempty(cond.pmod)
        error('tmfc_LSS_design: modulators are not produced by tmfc_LSS.');
    end
    if isempty(cond.onset) || isempty(cond.duration)
        % spm_jobman rejects empty onsets/durations (Size mismatch, required [Inf 1]) in tmfc_LSS
        error('tmfc_LSS_design:emptyOnsets','Condition ''%s'' has no onsets: LSS model needs at least one other event.',cond.name);
    end
    U(j).name = {cond.name}; %#ok<AGROW>
    U(j).ons  = cond.onset(:); %#ok<AGROW>
    U(j).dur  = cond.duration(:); %#ok<AGROW>
    U(j).orth = cond.orth; %#ok<AGROW>
    if isempty(U(j).orth), U(j).orth = true; end %#ok<AGROW>
    if length(U(j).dur) == 1 %#ok<ISCL>
        U(j).dur = repmat(U(j).dur,size(U(j).ons)); %#ok<AGROW>
    elseif numel(U(j).dur) ~= numel(U(j).ons)
        error('Mismatch between number of onset and number of durations.');
    end
    if any(isnan(U(j).ons)) || any(isnan(U(j).dur))
        error('Onsets or Durations contain NaN values.');
    end
    P.name = 'none';
    P.h    = 0;
    U(j).P = P; %#ok<AGROW>
    clear P
end
SPM.Sess(1).U = U;

% User-specified regressors (spm_run_fmri_spec.m:291-299)
C     = [];
Cname = {};
if isfield(sess,'regress')
    Cname = cell(1,numel(sess.regress));
    for q = 1:numel(sess.regress)
        Cname{q} = sess.regress(q).name;
        if numel(sess.regress(q).val(:)) ~= SPM.nscan(1)
            error('Length of regressor is not commensurate with data points.');
        end
        C = [C, sess.regress(q).val(:)]; %#ok<AGROW>
    end
end
SPM.Sess(1).C.C    = C;
SPM.Sess(1).C.name = Cname;

SPM.factor      = [];
SPM.xGX.iGXcalc = job.global;
SPM.xM.gMT      = job.mthresh;
SPM.xX.K(1).HParam = sess.hpf;
SPM.xVi.form    = job.cvi;

% Design matrix (spm_fmri_spm_ui.m:180); save_SPM = 0: no SPM.mat is written
SPM = spm_fMRI_design(SPM,0);
nscan = SPM.nscan;
nsess = length(nscan);

% High-pass filter structure (spm_fmri_spm_ui.m:198-214)
HParam = [SPM.xX.K(:).HParam];
if length(HParam) == 1 %#ok<ISCL>
    HParam = repmat(HParam,1,nsess);
end
for i = 1:nsess
    K(i) = struct('HParam', HParam(i),...
                  'row',    SPM.Sess(i).row,...
                  'RT',     SPM.xY.RT); %#ok<AGROW>
end
SPM.xX.K = spm_filter(K);

% Intrinsic autocorrelations Vi (spm_fmri_spm_ui.m:219-253)
cVi = SPM.xVi.form;
switch lower(cVi)
    case {'i.i.d', 'none'}
        SPM.xVi.V  = speye(sum(nscan));
        cVi        = 'i.i.d';
    case 'fast'
        SPM.xVi.Vi = spm_Ce('fast',nscan,SPM.xY.RT);
        cVi        = upper(cVi);
    otherwise
        SPM.xVi.Vi = spm_Ce('ar',nscan,0.2);
        cVi        = 'AR(0.2)';
end
SPM.xVi.form = cVi;

% Concatenated sessions (spm_fmri_concatenate.m:60-113)
if ~isempty(concat_scans)
    if sum(concat_scans) ~= SPM.nscan
        error('Number of scans does not match.');
    end
    SPM.nscan = concat_scans;
    if numel(SPM.xX.iB) == 1 && SPM.xX.iB == size(SPM.xX.X,2) %#ok<ISCL>
        Xb = [];
        Bn = {};
        for s = 1:numel(SPM.nscan)
            Xb    = blkdiag(Xb,ones(SPM.nscan(s),1));
            Bn{s} = sprintf('Sn(%i) constant',s); %#ok<AGROW>
        end
        SPM.xX.X    = [SPM.xX.X(:,1:end-1) Xb];
        SPM.xX.iB   = SPM.xX.iB:(SPM.xX.iB+size(Xb,2)-1);
        SPM.xX.name = {SPM.xX.name{1:end-1} Bn{:}}; %#ok<CCAT>
    end
    s = cumsum([0 SPM.nscan]);
    for i = 1:numel(SPM.nscan)
        K(i) = struct('HParam', SPM.xX.K(1).HParam,...
                      'row',    s(i) + (1:SPM.nscan(i)),...
                      'RT',     SPM.xY.RT);
    end
    SPM.xX.K = spm_filter(K);
    switch lower(SPM.xVi.form)
        case {'ar(1)','ar(0.2)'}
            SPM.xVi.Vi   = spm_Ce('ar',SPM.nscan,0.2);
            SPM.xVi.form = 'AR(0.2)';
        case 'fast'
            SPM.xVi.Vi   = spm_Ce('fast',SPM.nscan,SPM.xY.RT);
        case {'i.i.d', 'none'}
        otherwise
            warning('Unhandled temporal non-sphericity.');
    end
end
end
