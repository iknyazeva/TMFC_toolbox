function [sub_check,diag_info] = tmfc_LSS_after_FIR_fwl_batch(tmfc,start_sub,opts)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% LSS-after-FIR batched FWL solver. Same model, estimator and output
% files as the fixed reference tmfc_LSS_after_FIR
% and tmfc_LSS_after_FIR_in_memory; only the linear algebra on the data differs:
%
% Once per run (tmfc_LSS_after_FIR_data + the run-level part):
%   - the residual data Res_*.nii are read once, WITHOUT the global
%     re-scaling of the original main (FIR residuals are not re-scaled: gSF = 1 for both 'None' and
%     'Scaling'), with xM.TH = -Inf, the explicit mask.nii of the parent,
%     any(diff(Y,1)) and maxmem chunking;
%   - the nuisance part N of the after-FIR design is only the session
%     constant(s): ones(nScan,1), or blkdiag(ones(n_1,1),...,ones(n_S,1))
%     for concatenated runs (spm_fmri_concatenate.m:79-89); the filter K
%     (HParam = Inf, per original session for concat runs) leaves the data
%     unchanged, but spm_filter is called as in main;
%   - sN = spm_sp('Set',spm_filter(K,N)) and, for every chunk,
%     KY = spm_filter(K,Y), rN = spm_sp('r',sN,KY).
%
% Per trial:
%   - the design (tmfc_LSS_after_FIR_design) gives X(:,1:2) (Current_trial,
%     Other_trials); KXt = spm_filter(K,X(:,1:2));
%   - Fw = tmfc_LSS_after_FIR_fwl(KXt,sN), and the rank consistency check
%     against spm_sp: xKXs = spm_sp('Set',spm_filter(K,X)). If Fw is
%     deficient or rank(K*X) ~= rank(K*N) + rank(Z), the trial is
%     estimated exactly as the in-memory estimator (pKX = spm_sp('x-',xKXs); beta_1 =
%     pKX(1,:)*KY) - opts.fallback = true by default;
%   - for valid (full-rank) trials, the FWL weights of the first column
%     w_k = (Fw.Vs(1,:)/Fw.sc(1))*Fw.Q' (1 x nScan) and beta_1 = w_k*rN
%     (tmfc_LSS_after_FIR_fwl_apply); opts.multiply = 'batch' (default)
%     applies the weights of many trials with one matrix product
%     W_trials*rN per chunk (rows capped at nScan so the product stays
%     within the maxmem budget of one data chunk), 'trial' uses one
%     w_k*rN product per trial.
%
% Output (same names and format as tmfc_LSS_after_FIR):
%   LSS_regression_after_FIR/<subject>/Betas/Beta_[Run_r]_[Cond_c]_[name]_[Trial_n].nii
%   LSS_regression_after_FIR/<subject>/GLM_batches/GLM_[...].mat
%   (matlabbatch in GLM_batches is formed exactly as in the in-memory estimator and reference)
%
% Differences from tmfc_LSS_after_FIR:
%   - tmfc.defaults.parallel = 1 is ignored with a warning (trials run
%     sequentially);
%   - the in-mask residuals of one run are kept in RAM as in the in-memory estimator.
%
% FORMAT [sub_check,diag_info] = tmfc_LSS_after_FIR_fwl_batch(tmfc)
% FORMAT [sub_check,diag_info] = tmfc_LSS_after_FIR_fwl_batch(tmfc,start_sub)
% FORMAT [sub_check,diag_info] = tmfc_LSS_after_FIR_fwl_batch(tmfc,start_sub,opts)
%
%   opts.multiply - 'batch' (default) or 'trial', see above
%   opts.fallback - true (default): rank-inconsistent trials use the exact
%                   in-memory estimator
%   diag_info(iSub).session(r).trial(k) - nMask, fallback, rk_fwl, rk_spm,
%                   min_sv (smallest scaled singular value of Z)
%
% Gate: maximum beta difference from the reference <= 1e-10 in double.
%
% =========================================================================
% Derived from tmfc_LSS_after_FIR.m (Copyright (C) 2026 Ruslan Masharipov)
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if nargin < 3, opts = struct(); end
if nargin == 1
   start_sub = 1;
end
if ~isfield(opts,'multiply'), opts.multiply = 'batch'; end
if ~isfield(opts,'fallback'), opts.fallback = true; end
if ~strcmp(opts.multiply,'batch') && ~strcmp(opts.multiply,'trial')
    error('tmfc_LSS_after_FIR_fwl_batch: unknown opts.multiply ''%s''.',opts.multiply);
end

% Check subject names
if ~isfield(tmfc.subjects,'name')
    for iSub = 1:length(tmfc.subjects)
        tmfc.subjects(iSub).name = ['Subject_' num2str(iSub,'%04.f')];
    end
end

% Check session mode
if ~isfield(tmfc,'session_mode')
    session_mode = 'strict';
else
    session_mode = tmfc.session_mode;
end

try
    main_GUI = guidata(findobj('Tag','TMFC_GUI'));
    set(main_GUI.TMFC_GUI_S10,'String', 'Updating...','ForegroundColor',[0.772, 0.353, 0.067]);
end

spm('defaults','fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);
spm_get_defaults('stats.resmem',tmfc.defaults.resmem);
spm_get_defaults('stats.maxmem',tmfc.defaults.maxmem);
spm_get_defaults('stats.fmri.ufp',1);

if tmfc.defaults.parallel == 1
    warning('tmfc_LSS_after_FIR_fwl_batch:ParallelNotImplemented', ...
        'parallel computing is not implemented by the batched FWL solver, trials run sequentially.');
end

nSub = length(tmfc.subjects);

cond_list = tmfc.LSS_after_FIR.conditions;
nCond = length(cond_list);

% Initialize waitbar
w = waitbar(0,'Please wait...','Name','LSS regression after FIR','Tag','tmfc_waitbar');
start_time = tic;
count_sub = 1;
cleanupObj = onCleanup(@unfreeze_after_ctrl_c);

diag_info = struct('session',{});

% Loop through subjects
for iSub = start_sub:nSub
    SPM = load(tmfc.subjects(iSub).path).SPM;

    if strcmp(session_mode,'repeated')
        sess_num = 1:length(SPM.Sess);
    else
        sess_num = unique([cond_list.sess]);
    end

    nSess = length(sess_num);

    % Check if SPM.mat has concatenated sessions
    % (if spm_fmri_concatenate.m script was used)
    if size(SPM.nscan,2) == size(SPM.Sess,2)
        SPM_concat = 0;
    else
        SPM_concat = 1;
    end

    sub_dir = fullfile(tmfc.project_path,'LSS_regression_after_FIR',tmfc.subjects(iSub).name);
    if isfolder(sub_dir)
        rmdir(sub_dir,'s');
        pause(0.1);
    end
    mkdir(fullfile(sub_dir,'Betas'));
    mkdir(fullfile(sub_dir,'GLM_batches'));

    % Loop through sessions
    for jSess = 1:nSess

        % Trials of interest
        nTrial = 0;
        ons_of_int = [];
        dur_of_int = [];
        cond_of_int = [];
        trial.cond = [];
        trial.number = [];
        for kCond = 1:nCond
            if strcmp(session_mode,'repeated') || cond_list(kCond).sess == sess_num(jSess)
                nTrial = nTrial + length(SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons);
                ons_of_int = [ons_of_int; SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons];
                dur_of_int = [dur_of_int; SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).dur];
                cond_of_int = [cond_of_int cond_list(kCond).number];
                trial.cond = [trial.cond; repmat(cond_list(kCond).number,length(SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons),1)];
                trial.number = [trial.number; (1:length(SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons))'];
            end
        end

        all_trials_number = (1:nTrial)';

        % Trials of no interest
        cond_of_no_int = setdiff((1:length(SPM.Sess(sess_num(jSess)).U)),cond_of_int);
        ons_of_no_int = [];
        dur_of_no_int = [];
        for kCondNoInt = 1:length(cond_of_no_int)
            ons_of_no_int = [ons_of_no_int; SPM.Sess(sess_num(jSess)).U(cond_of_no_int(kCondNoInt)).ons];
            dur_of_no_int = [dur_of_no_int; SPM.Sess(sess_num(jSess)).U(cond_of_no_int(kCondNoInt)).dur];
        end

        % Residual images of the run (tmfc_LSS_after_FIR.m:230-238)
        scans = afir_scans(tmfc,SPM,sess_num(jSess),jSess,SPM_concat,iSub);

        % Read the run once
        job0 = afir_job(tmfc,SPM,sess_num(jSess),iSub,1,scans,[],[],[],[]);
        D    = tmfc_LSS_after_FIR_data(job0);
        if SPM_concat == 1
            concat_scans = SPM.nscan;
        else
            concat_scans = [];
        end
        nChunk = numel(D.Y);

        % ------- per-trial designs, FWL weights and rank consistency -------
        jobs      = cell(nTrial,1);
        descrip   = cell(nTrial,1);
        fallback  = false(nTrial,1);
        pKX_fb    = cell(nTrial,1);
        rk_spm    = zeros(nTrial,1);
        rk_fwl    = zeros(nTrial,1);
        min_sv    = nan(nTrial,1);
        wAll      = [];
        for kTrial = 1:nTrial

            % Current trial vs all other trials (of interest and no interest)
            other_trials = all_trials_number(all_trials_number~=kTrial);
            jobs{kTrial} = afir_job(tmfc,SPM,sess_num(jSess),iSub,kTrial,scans, ...
                ons_of_int(kTrial),dur_of_int(kTrial), ...
                [ons_of_int(other_trials); ons_of_no_int],[dur_of_int(other_trials); dur_of_no_int]);

            % Design of this trial in memory
            SPM_LSS = tmfc_LSS_after_FIR_design(jobs{kTrial},concat_scans);

            if kTrial == 1
                K     = SPM_LSS.xX.K;
                N     = SPM_LSS.xX.X(:,3:end);
                sN    = spm_sp('Set',spm_filter(K,N));
                sN.X  = full(sN.X);
                rkN   = spm_sp('rk',sN);
                nScan = size(SPM_LSS.xX.X,1);
                wAll  = nan(nTrial,nScan);
            else
                % the nuisance part is the session constant(s), same for all trials
                if ~isequal(SPM_LSS.xX.X(:,3:end),N)
                    error('tmfc_LSS_after_FIR_fwl_batch: nuisance part of the design differs from the run-level part.');
                end
            end

            % rank of K*X by spm_sp (also the fallback estimator)
            xKXs   = spm_sp('Set',spm_filter(K,SPM_LSS.xX.X));
            xKXs.X = full(xKXs.X);
            rk_spm(kTrial) = spm_sp('rk',xKXs);

            % FWL basis of the trial columns
            Fw = tmfc_LSS_after_FIR_fwl(spm_filter(K,SPM_LSS.xX.X(:,1:2)),sN);
            rk_fwl(kTrial) = rkN + Fw.rk;
            min_sv(kTrial) = min(Fw.sv);
            fallback(kTrial) = Fw.deficient || rk_spm(kTrial) ~= rk_fwl(kTrial);

            if fallback(kTrial)
                if opts.fallback
                    % exact in-memory estimator for this trial
                    pKX_fb{kTrial} = spm_sp('x-',xKXs);
                else
                    [~,wAll(kTrial,:)] = tmfc_LSS_after_FIR_fwl_apply(Fw,[]);
                end
            else
                [~,wAll(kTrial,:)] = tmfc_LSS_after_FIR_fwl_apply(Fw,[]);
            end
            descrip{kTrial} = sprintf('spm_spm:beta (%04d) - %s',1,SPM_LSS.xX.name{1});
        end

        % ------- data pass: rN once per chunk, then the trial weights -------
        useW      = ~isnan(wAll(:,1));              % trials with FWL weights
        useExact  = fallback & opts.fallback;       % trials estimated as the in-memory estimator
        beta      = cell(nTrial,1);
        for kTrial = 1:nTrial, beta{kTrial} = cell(nChunk,1); end
        for i = 1:nChunk
            if isempty(D.Y{i})
                % chunk entirely outside the analysis mask: spm_sp('r',sN,[])
                % would return the full projector [nScan x nScan], so the
                % trial weights must not be applied here
                for k = 1:nTrial, beta{k}{i} = zeros(1,0); end
                continue
            end
            KY = spm_filter(K,D.Y{i});
            rN = spm_sp('r',sN,KY);

            % FWL trials
            widx = find(useW);
            if strcmp(opts.multiply,'trial')
                for jj = 1:numel(widx)
                    k = widx(jj);
                    beta{k}{i} = wAll(k,:)*rN;
                end
            else
                % batched product, rows capped at nScan so that the result
                % stays within the maxmem budget of one data chunk
                j0 = 1;
                while j0 <= numel(widx)
                    jb = j0:min(j0+nScan-1,numel(widx));
                    B  = wAll(widx(jb),:)*rN;
                    for jj = 1:numel(jb)
                        beta{widx(jb(jj))}{i} = B(jj,:);
                    end
                    j0 = jb(end) + 1;
                end
            end

            % exact fallback trials (rank-deficient / inconsistent rank)
            eidx = find(useExact);
            for jj = 1:numel(eidx)
                k = eidx(jj);
                beta{k}{i} = pKX_fb{k}(1,:)*KY;
            end
        end

        % ------- write the per-trial outputs, as in the in-memory estimator and reference -------
        for kTrial = 1:nTrial
            cond_name = regexprep(char(SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name(1)),' ','_');
            write_beta(fullfile(sub_dir,'Betas', ...
                ['Beta_[Run_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' cond_name ']_[Trial_' num2str(trial.number(kTrial)) '].nii']), ...
                D,beta{kTrial},descrip{kTrial});

            matlabbatch = {struct('spm',struct('stats',struct('fmri_spec',jobs{kTrial})))};
            parsave_batch(fullfile(sub_dir,'GLM_batches', ...
                ['GLM_[Run_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' cond_name ']_[Trial_' num2str(trial.number(kTrial)) '].mat']),matlabbatch);

            condition(trial.cond(kTrial)).trials(trial.number(kTrial)) = 1;

            diag_info(iSub).session(sess_num(jSess)).trial(kTrial) = struct( ...
                'nMask',D.nIn,'fallback',fallback(kTrial), ...
                'rk_fwl',rk_fwl(kTrial),'rk_spm',rk_spm(kTrial),'min_sv',min_sv(kTrial));
        end

        sub_check(iSub).session(sess_num(jSess)).condition = condition;
        pause(0.01);

        clear nTrial ons* dur* cond_of_int cond_of_no_int trial all_trials_number condition D scans job0 jobs beta fallback pKX_fb wAll
    end

    % Update main TMFC GUI
    try
        main_GUI = guidata(findobj('Tag','TMFC_GUI'));
        set(main_GUI.TMFC_GUI_S10,'String', strcat(num2str(iSub), '/', num2str(nSub), ' done'),'ForegroundColor',[0.219, 0.341, 0.137]);
    end

    % Update waitbar
    elapsed_time = toc(start_time);
    time_per_sub = elapsed_time/count_sub;
    count_sub = count_sub + 1;
    time_remaining = (nSub-iSub)*time_per_sub;
    hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
    try
        waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
    end

    clear SPM
end

% Close waitbar
try
    delete(w);
end

function unfreeze_after_ctrl_c()
    try
        delete(findall(0,'type','figure','Tag', 'tmfc_waitbar'));
        GUI = guidata(findobj('Tag','TMFC_GUI'));
        set([GUI.TMFC_GUI_B1, GUI.TMFC_GUI_B2, GUI.TMFC_GUI_B3, GUI.TMFC_GUI_B4,...
           GUI.TMFC_GUI_B5a, GUI.TMFC_GUI_B5b, GUI.TMFC_GUI_B6, GUI.TMFC_GUI_B7,...
           GUI.TMFC_GUI_B8, GUI.TMFC_GUI_B9, GUI.TMFC_GUI_B10, GUI.TMFC_GUI_B11,...
           GUI.TMFC_GUI_B12a,GUI.TMFC_GUI_B12b,GUI.TMFC_GUI_B13a,GUI.TMFC_GUI_B13b,...
           GUI.TMFC_GUI_B14a, GUI.TMFC_GUI_B14b], 'Enable', 'on');
    end
end
end

%% ========================================================================

% Residual scan list of one run (tmfc_LSS_after_FIR.m:230-238)
function scans = afir_scans(tmfc,SPM,sess_num,jSess,SPM_concat,iSub)
    if SPM_concat == 0
        for image = 1:SPM.nscan(sess_num)
            scans{image,1} = fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                ['Res_' num2str(SPM.Sess(sess_num).row(image),'%.4d') '.nii,1']);
        end
    else
        for image = 1:size(SPM.xY.VY,1)
            scans{image,1} = fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                ['Res_' num2str(SPM.Sess(jSess).row(image),'%.4d') '.nii,1']);
        end
    end
end

% fmri_spec job of one LSS-after-FIR trial (tmfc_LSS_after_FIR.m:222-278)
function job = afir_job(tmfc,SPM,sess_num,iSub,kTrial,scans,cur_ons,cur_dur,oth_ons,oth_dur)
    run_dir = fullfile(tmfc.project_path,'LSS_regression_after_FIR',tmfc.subjects(iSub).name, ...
        ['LSS_Run_' num2str(sess_num) '_Trial_' num2str(kTrial)]);
    job.dir = {run_dir};
    job.timing.units = SPM.xBF.UNITS;
    job.timing.RT = SPM.xY.RT;
    job.timing.fmri_t = SPM.xBF.T;
    job.timing.fmri_t0 = SPM.xBF.T0;

    job.sess.scans = scans;

    % Conditions
    job.sess.cond(1).name = 'Current_trial';
    job.sess.cond(1).onset = cur_ons;
    job.sess.cond(1).duration = cur_dur;
    job.sess.cond(1).tmod = 0;
    job.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    job.sess.cond(1).orth = 1;
    job.sess.cond(2).name = 'Other_trials';
    job.sess.cond(2).onset = oth_ons;
    job.sess.cond(2).duration = oth_dur;
    job.sess.cond(2).tmod = 0;
    job.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    job.sess.cond(2).orth = 1;

    % No confounds (already regressed out during the FIR step)
    job.sess.regress = struct('name', {}, 'val', {});
    job.sess.multi = {''};
    job.sess.multi_reg = {''};

    % HPF, HRF, mask
    job.sess.hpf = Inf;
    job.fact = struct('name', {}, 'levels', {});
    job.bases.hrf.derivs = [0 0];
    job.volt = 1;
    job.global = SPM.xGX.iGXcalc;
    job.mthresh = -Inf;
    job.mask = {fullfile(SPM.swd,'mask.nii')};
    job.cvi = 'None';
end

% Write the first beta image as spm_spm does (spm_spm.m:469-484, 611-616)
function write_beta(fname,D,beta1,descrip)
    Vbeta = struct(...
        'fname',   fname,...
        'dim',     D.DIM,...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'mat',     D.M,...
        'pinfo',   [1 0 0]',...
        'descrip', descrip);
    Vbeta = spm_data_hdr_write(Vbeta);
    for i = 1:numel(D.chunk)
        c = NaN(numel(D.chunk{i}),1);
        c(D.cmask{i}) = beta1{i};
        Vbeta = spm_data_write(Vbeta, c, D.chunk{i});
    end
end

% Save batches
function parsave_batch(fname,matlabbatch)
	save(fname, 'matlabbatch')
end
