function [sub_check] = tmfc_LSS(tmfc,start_sub)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% For each individual trial, the Least-Squares Separate (LSS) approach
% estimates a separate GLM with two regressors. The first regressor models
% the expected BOLD response to the current trial of interest, and the 
% second (nuisance) regressor models the BOLD response to all other trials
% (of interest and no interest). For trials of no interest (e.g., errors),
% individual GLMs are not estimated. Trials of no interest are used only
% for the second (nuisance) regressor.
%
% This function uses SPM.mat file (which contains the specification of the
% 1st-level GLM with canonical HRF) to specify and estimate 1st-level GLMs
% for each individual trial of interest (LSS approach).
%
% Note: If the original GLMs contain parametric or time modulators, they
% will be removed from the LSS GLMs. If the original GLMs contain time and
% dispersion derivatives, they will be removed from the LSS GLMs.
%
% FORMAT [sub_check] = tmfc_LSS(tmfc)
% Run a function starting from the first subject in the list.
%
%   tmfc.subjects.path     - Paths to individual SPM.mat files
%   tmfc.subjects.name     - Subject names within the TMFC project
%   tmfc.project_path      - Path where all results will be saved
%   tmfc.defaults.parallel - 0 or 1 (LSS always runs serially; 1 warns)
%   tmfc.defaults.lss_method - 'exact' (default) or 'pooled' (approximate)
%   tmfc.defaults.lss_reml_warm_start - pooled: default true for FAST,
%                            false otherwise; explicit logical overrides
%   tmfc.defaults.lss_reml_every - pooled ReML interval, default 1
%   tmfc.defaults.maxmem   - e.g. 2^31 = 2GB (how much RAM can be used at
%                            the same time during GLM estimation)
%   tmfc.defaults.resmem   - SPM setting; LSS run data are held in memory
%
%   tmfc.session_mode      - 'strict' or 'repeated'
%                            'strict'   - use the selected run-specific conditions
%                            'repeated' - process the same conditions across all runs
%
%   tmfc.LSS.conditions           - List of conditions of interest
%   tmfc.LSS.conditions.sess      - Session number 
%                                   (as specified in SPM.Sess)
%   tmfc.LSS.conditions.number    - Condition number
%                                   (as specified in SPM.Sess.U)
%   tmfc.LSS.conditions.name      - Condition name
%                                   (as specified in SPM.Sess.U.name)
%   tmfc.LSS.conditions.file_name - Condition-specific file names:
%   (['[Run_' num2str(iSess) ']_[Cond_' num2str(jCond) ']_[' ...
%    regexprep(char(SPM.Sess(iSess).U(jCond).name(1)),' ','_') ']'];)
%
% Run (session) number and condition number must match the original SPM.mat file.
% Consider, for example, a task design with two runs (SPM sessions). Both runs 
% contain three task regressors for "Cond A", "Cond B" and "Errors". If
% you are only interested in comparing "Cond A" and "Cond B", the following
% structure must be specified (see tmfc_conditions_GUI and tmfc_get_conditions):
%
%   tmfc.LSS.conditions(1).sess   = 1;     
%   tmfc.LSS.conditions(1).number = 1; 
%   tmfc.LSS.conditions(1).name = 'Cond_A';
%   tmfc.LSS.conditions(1).file_name = '[Run_1]_[Cond_1]_[Cond_A]';  
%
%   tmfc.LSS.conditions(2).sess   = 1;
%   tmfc.LSS.conditions(2).number = 2;
%   tmfc.LSS.conditions(2).name = 'Cond_B';
%   tmfc.LSS.conditions(2).file_name = '[Run_1]_[Cond_2]_[Cond_B]';  
%
%   tmfc.LSS.conditions(3).sess   = 2;
%   tmfc.LSS.conditions(3).number = 1;
%   tmfc.LSS.conditions(3).name = 'Cond_A';
%   tmfc.LSS.conditions(3).file_name = '[Run_2]_[Cond_1]_[Cond_A]';  
%
%   tmfc.LSS.conditions(4).sess   = 2;
%   tmfc.LSS.conditions(4).number = 2;
%   tmfc.LSS.conditions(4).name = 'Cond_B';
%   tmfc.LSS.conditions(4).file_name = '[Run_2]_[Cond_2]_[Cond_B]'; 
%
% FORMAT [sub_check] = tmfc_LSS(tmfc, start_sub)
% Run the function starting from a specific subject in the path list.
%
%   tmfc           - As above
%   start_sub      - Subject number in the list to start computations from
%
% This in-memory implementation replaces the job-based estimation in main.
% Exact estimates retain the per-trial covariance model.
% Pooled uses run-level sample covariance with per-trial ReML and whitening.
% rWLS pooled shares Cy per run and keeps rWLS ReML for each trial.
%
% =========================================================================
% Copyright (C) 2026 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if nargin == 1
   start_sub = 1;
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

try %#ok<TRYNC>
    main_GUI = guidata(findobj('Tag','TMFC_GUI'));
    set(main_GUI.TMFC_GUI_S6,'String', 'Updating...','ForegroundColor',[0.772, 0.353, 0.067]);
end

spm('defaults','fmri');
spm_get_defaults('cmdline',true);
spm_get_defaults('stats.resmem',tmfc.defaults.resmem);
spm_get_defaults('stats.maxmem',tmfc.defaults.maxmem);
spm_get_defaults('stats.fmri.ufp',1);
UFp = spm_get_defaults('stats.fmri.ufp');

method = 'exact';
if isfield(tmfc.defaults,'lss_method'), method = tmfc.defaults.lss_method; end
if ~any(strcmp(method,{'exact','pooled'}))
    error('tmfc_LSS:InvalidMethod','lss_method must be exact or pooled.');
end
warm = false; every = 1;
has_options = isfield(tmfc.defaults,'lss_reml_warm_start') || isfield(tmfc.defaults,'lss_reml_every');
if isfield(tmfc.defaults,'lss_reml_warm_start'), warm = tmfc.defaults.lss_reml_warm_start; end
if isfield(tmfc.defaults,'lss_reml_every'), every = tmfc.defaults.lss_reml_every; end
if has_options && strcmp(method,'exact')
    warning('tmfc_LSS:PooledOptionsIgnored','Pooled ReML options are ignored for exact LSS.');
end
warned_rwls = false;
if tmfc.defaults.parallel == 1
    warning('tmfc_LSS:SerialOnly','LSS runs serially; parallel is ignored.');
end

nSub = length(tmfc.subjects);

cond_list = tmfc.LSS.conditions;
nCond = length(cond_list);

% Initialize waitbar
w = waitbar(0,'Please wait...','Name','LSS regression','Tag','tmfc_waitbar');
start_time = tic;
count_sub = 1;
cleanupObj = onCleanup(@unfreeze_after_ctrl_c);


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

    is_rwls = strcmp(SPM.xVi.form,'wls');
    if ~isfield(tmfc.defaults,'lss_reml_warm_start')
        warm = strcmp(method,'pooled') && strcmpi(SPM.xVi.form,'fast');
    end
    if strcmp(method,'pooled') && ~is_rwls
        if ~islogical(warm) || ~isscalar(warm)
            error('tmfc_LSS:InvalidWarmStart','lss_reml_warm_start must be a logical scalar.');
        end
        if ~isnumeric(every) || ~isreal(every) || ~isscalar(every) || ~isfinite(every) || every < 1 || every ~= fix(every)
            error('tmfc_LSS:InvalidEvery','lss_reml_every must be a positive integer.');
        end
    end
    if is_rwls && strcmp(method,'pooled') && has_options && ~warned_rwls
        warning('tmfc_LSS:PooledRWLS','Pooled ReML options are unsupported for rWLS and are ignored.');
        warned_rwls = true;
    end

    sub_dir = fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name);
    if isdir(sub_dir) %#ok<ISDIR>
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
                ons_of_int = [ons_of_int; SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons]; %#ok<AGROW>
                dur_of_int = [dur_of_int; SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).dur]; %#ok<AGROW>
                cond_of_int = [cond_of_int cond_list(kCond).number]; %#ok<AGROW>
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
            ons_of_no_int = [ons_of_no_int; SPM.Sess(sess_num(jSess)).U(cond_of_no_int(kCondNoInt)).ons]; %#ok<AGROW>
            dur_of_no_int = [dur_of_no_int; SPM.Sess(sess_num(jSess)).U(cond_of_no_int(kCondNoInt)).dur]; %#ok<AGROW>
        end

        % Functional images of the run
        scans = lss_scans(SPM,sess_num(jSess),jSess,SPM_concat);

        % Read the run once (job fields used by tmfc_LSS_load_data are the same for all trials)
        job = lss_job(tmfc,SPM,sess_num(jSess),iSub,1,scans,[],[],[],[]);
        D   = tmfc_LSS_load_data(job);
        if SPM_concat == 1
            concat_scans = SPM.nscan;
        else
            concat_scans = [];
        end

        % Loop through trials of interest
        for kTrial = 1:nTrial

            % Current trial vs all other trials (of interest and no interest)
            other_trials = all_trials_number(all_trials_number~=kTrial);
            job = lss_job(tmfc,SPM,sess_num(jSess),iSub,kTrial,scans, ...
                ons_of_int(kTrial),dur_of_int(kTrial), ...
                [ons_of_int(other_trials); ons_of_no_int],[dur_of_int(other_trials); dur_of_no_int]);

            % Specify and estimate LSS GLM in memory
            SPM_LSS = tmfc_LSS_design(job,concat_scans);
            if is_rwls
                if kTrial == 1
                    Rrwls.N = SPM_LSS.xX.X(:,3:end);
                    Rrwls.K = SPM_LSS.xX.K;
                    Rrwls.KY = cell(size(D.Y));
                    for iChunk = 1:numel(D.Y)
                        Rrwls.KY{iChunk} = spm_filter(Rrwls.K,D.Y{iChunk});
                    end
                    Rrwls.Cy = [];
                    if strcmp(method,'pooled')
                        Rrwls.Cy = tmfc_rwls_LSS_pooled_cy(SPM_LSS,D,Rrwls.KY);
                    end
                end
                if ~isequal(SPM_LSS.xX.X(:,3:end),Rrwls.N) || ~isequal(SPM_LSS.xX.K,Rrwls.K)
                    error('tmfc_LSS:RWLSRunDesignChanged','Run-level nuisance design or filter changed.');
                end
                beta1 = tmfc_rwls_LSS(SPM_LSS,D,UFp,Rrwls.KY,Rrwls.Cy);
            else
                if kTrial == 1
                    R = tmfc_LSS_exact_prepare(SPM_LSS,D);
                    hParent = []; keepParent = []; Wprev = [];
                    if strcmp(method,'pooled') && warm && ~isfield(SPM_LSS.xVi,'V')
                        [hParent,keepParent] = tmfc_LSS_pooled_parent(SPM,SPM_LSS,sess_num(jSess));
                        if ~strcmpi(SPM.xVi.form,'fast'), keepParent = []; end
                    end
                    if strcmp(method,'pooled') && ~isfield(SPM_LSS.xVi,'V')
                        SPM_run = SPM_LSS;
                        SPM_run.xX.X = R.N;
                        Cy = tmfc_LSS_pooled_cy(SPM_run,D,R,UFp);
                    else
                        Cy = [];
                    end
                end
                if strcmp(method,'pooled')
                    if mod(kTrial-1,every) == 0, Wprev = []; end
                    [beta1,Wprev] = tmfc_LSS_pooled_estimate(SPM_LSS,D,Cy,R,hParent,keepParent,Wprev);
                else
                    beta1 = tmfc_LSS_exact_estimate(SPM_LSS,D,R,UFp);
                end
            end

            % Save individual trial beta image
            cond_name = regexprep(char(SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name(1)),' ','_');
            write_beta(fullfile(sub_dir,'Betas', ...
                ['Beta_[Run_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' cond_name ']_[Trial_' num2str(trial.number(kTrial)) '].nii']), ...
                D,beta1,sprintf('spm_spm:beta (%04d) - %s',1,SPM_LSS.xX.name{1}));

            % Save GLM_batch.mat file (marker that the trial is done)
            matlabbatch = {struct('spm',struct('stats',struct('fmri_spec',job)))};
            tmfc_LSS_save_batch(fullfile(sub_dir,'GLM_batches', ...
                ['GLM_[Run_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' cond_name ']_[Trial_' num2str(trial.number(kTrial)) '].mat']),matlabbatch);

            condition(trial.cond(kTrial)).trials(trial.number(kTrial)) = 1; %#ok<AGROW>


        end

        sub_check(iSub).session(sess_num(jSess)).condition = condition; %#ok<AGROW>
        pause(0.01);

        clear nTrial ons* dur* cond_of_int cond_of_no_int trial all_trials_number condition D R Cy scans
    end

    % Update main TMFC GUI
    try %#ok<TRYNC>
        main_GUI = guidata(findobj('Tag','TMFC_GUI'));
        set(main_GUI.TMFC_GUI_S6,'String', strcat(num2str(iSub), '/', num2str(nSub), ' done'),'ForegroundColor',[0.219, 0.341, 0.137]);
    end

    % Update waitbar
    elapsed_time = toc(start_time);
    time_per_sub = elapsed_time/count_sub;
    count_sub = count_sub + 1;
    time_remaining = (nSub-iSub)*time_per_sub;
    hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
    try %#ok<TRYNC>
        waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
    end

    clear SPM
end

% Close waitbar
try %#ok<TRYNC>
    delete(w);
end

function unfreeze_after_ctrl_c()
    try %#ok<TRYNC>
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

% Scan list of one run (tmfc_LSS.m:208-219)
function scans = lss_scans(SPM,sess_num,jSess,SPM_concat)
    if SPM_concat == 0
        for image = 1:SPM.nscan(sess_num)
            scans{image,1} = [SPM.xY.VY(SPM.Sess(sess_num).row(image)).fname ',' ...
                              num2str(SPM.xY.VY(SPM.Sess(sess_num).row(image)).n(1))]; %#ok<AGROW>
        end
    else
        for image = 1:size(SPM.xY.VY,1)
            scans{image,1} = [SPM.xY.VY(SPM.Sess(jSess).row(image)).fname ',' ...
                              num2str(SPM.xY.VY(SPM.Sess(jSess).row(image)).n(1))]; %#ok<AGROW>
        end
    end
end

% fmri_spec job of one LSS trial (tmfc_LSS.m:202-271)
function job = lss_job(tmfc,SPM,sess_num,iSub,kTrial,scans,cur_ons,cur_dur,oth_ons,oth_dur)
    run_dir = fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Run_' num2str(sess_num) '_Trial_' num2str(kTrial)]);
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

    % Confounds
    for conf = 1:length(SPM.Sess(sess_num).C.name)
        job.sess.regress(conf).name = SPM.Sess(sess_num).C.name{1,conf};
        job.sess.regress(conf).val = SPM.Sess(sess_num).C.C(:,conf);
    end

    job.sess.multi = {''};
    job.sess.multi_reg = {''};

    % HPF, HRF, mask
    job.sess.hpf = SPM.xX.K(sess_num).HParam;
    job.fact = struct('name', {}, 'levels', {});
    job.bases.hrf.derivs = [0 0];
    job.volt = 1;
    job.global = SPM.xGX.iGXcalc;
    job.mthresh = SPM.xM.gMT;

    try
        job.mask = {SPM.xM.VM.fname};
    catch
        job.mask = {''};
    end

    if strcmp(SPM.xVi.form,'i.i.d') || strcmp(SPM.xVi.form,'none')
        job.cvi = 'None';
    elseif strcmp(SPM.xVi.form,'fast') || strcmp(SPM.xVi.form,'FAST')
        job.cvi = 'FAST';
    else
        job.cvi = 'AR(1)';
    end
end

% Write the first beta image as spm_spm does (spm_spm.m:459-484, 611-616)
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

% Save the original trial batch.
function tmfc_LSS_save_batch(fname,matlabbatch)
	save(fname, 'matlabbatch')
end
