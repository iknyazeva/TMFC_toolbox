function sub_check = tmfc_LSS(tmfc,start_sub,end_sub)

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
% Functional data are loaded once per run and reused across trials.
% GLMs are specified and estimated in memory, retaining separate
% ReML/rWLS covariance estimation for each trial.
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
%   tmfc.defaults.parallel - 0 or 1 (sequential or parallel computing)
%   tmfc.defaults.headless - true or false (disable GUI and waitbars)
%   tmfc.defaults.maxmem   - Memory limit (bytes) for data caching and GLM estimation
%   tmfc.defaults.ufp      - F-test threshold for ReML voxel selection
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
% FORMAT [sub_check] = tmfc_LSS(tmfc, start_sub, end_sub)
% Run the function for a specified range of subjects (inclusive).
%
%   tmfc           - As above
%   start_sub      - First subject number
%   end_sub        - Last subject number
%
% =========================================================================
% Copyright (C) 2026 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

% Check subject range and defaults
if nargin < 2 || isempty(start_sub)
    start_sub = 1;
end
if ~isfield(tmfc,'defaults'), tmfc.defaults = struct(); end
if ~isfield(tmfc.defaults,'headless'), tmfc.defaults.headless = false; end
if ~isfield(tmfc.defaults,'parallel'), tmfc.defaults.parallel = 0; end
if ~isfield(tmfc.defaults,'maxmem'), tmfc.defaults.maxmem = 2^31; end
if ~isfield(tmfc.defaults,'ufp'), tmfc.defaults.ufp = 1; end

% Check subject names
nSub = numel(tmfc.subjects);
if ~isfield(tmfc.subjects,'name')
    for iSub = 1:nSub
        tmfc.subjects(iSub).name = sprintf('Subject_%04d',iSub);
    end
end

% Check session mode
if ~isfield(tmfc,'session_mode'), tmfc.session_mode = 'strict'; end
if ~any(strcmp(tmfc.session_mode,{'strict','repeated'}))
    error('tmfc_LSS:SessionMode','session_mode must be strict or repeated.');
end

validateattributes(tmfc.defaults.maxmem,{'numeric'},{'scalar','finite','positive'});
validateattributes(tmfc.defaults.ufp,{'numeric'},{'scalar','>',0,'<=',1});
if nargin < 3 || isempty(end_sub)
    end_sub = nSub;
end
validateattributes(start_sub,{'numeric'},{'scalar','integer','>=',1,'<=',nSub});
validateattributes(end_sub,{'numeric'},{'scalar','integer','>=',start_sub,'<=',nSub});
nSubRun = end_sub-start_sub+1;
sub_check = repmat(struct('session',[]),1,nSub);

% Initialize SPM
spm('defaults','fmri');
spm_get_defaults('cmdline',true);

% Initialize parfor
useParSub = false;
useParTrial = false;
hasPCT = (exist('parfor','builtin')==5) && license('test','Distrib_Computing_Toolbox');
if tmfc.defaults.parallel==1 && hasPCT
    p = gcp('nocreate');
    if isempty(p)
        parpool;
        p = gcp('nocreate');
    end
    nW = p.NumWorkers;

    % Many subjects - parallelize over subjects
    % Few subjects - parallelize over trials
    useParSub = (nSubRun >= nW);
    useParTrial = ~useParSub;
end

% Initialize GUI progress reporting only in interactive mode
showGUI = ~tmfc.defaults.headless && usejava('desktop');
D = []; w = [];
start_time = tic;
if showGUI
    try, figure(findobj('Tag','TMFC_GUI')); end
    try
        main_GUI = guidata(findobj('Tag','TMFC_GUI'));
        set(main_GUI.TMFC_GUI_S6,'String','Updating...', ...
            'ForegroundColor',[0.772 0.353 0.067]);
    end
    cleanupObj = onCleanup(@unfreeze_after_ctrl_c); 

    if ~useParSub
        w = waitbar(0,'Please wait...','Name','LSS regression','Tag','tmfc_waitbar');
    else
        try % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;
            w = waitbar(0,'Please wait...','Name','LSS regression','Tag','tmfc_waitbar');
            afterEach(D,@tmfc_LSS_parfor_waitbar);
            tmfc_LSS_parfor_waitbar(w,nSubRun);
        catch % No waitbar for MATLAB R2016b and earlier
            D = [];
            opts = struct('WindowStyle','non-modal','Interpreter','tex');
            w = warndlg({'\fontsize{12}Sorry, waitbar progress update is not available for parallel computations in MATLAB R2016b and earlier.',[],...
                'Please wait until all computations are completed.',[],...
                'If you want to interrupt computations:',...
                '   1) Do not close this window;',...
                '   2) Select MATLAB main window;',...
                '   3) Press Ctrl+C.'},'Please wait...',opts);
        end
    end
end

% -------------------------------------------------------------------------
% Subject loop
% -------------------------------------------------------------------------
try
    % ===================== Subjects: sequential ==========================
    if ~useParSub
        count_sub = 1;
        for iSub = start_sub:end_sub
            sub_check(iSub) = tmfc_LSS_one_subject(tmfc,iSub,useParTrial);

            % Update GUI progress
            if showGUI
                try
                    main_GUI = guidata(findobj('Tag','TMFC_GUI'));
                    set(main_GUI.TMFC_GUI_S6,'String', ...
                        [num2str(iSub) '/' num2str(nSub) ' done'], ...
                        'ForegroundColor',[0.219 0.341 0.137]);
                end
                elapsed_time = toc(start_time);
                time_per_sub = elapsed_time/count_sub;
                count_sub = count_sub+1;
                time_remaining = (end_sub-iSub)*time_per_sub;
                hms = fix(mod(time_remaining,[0 3600 60])./[3600 60 1]);
                progress = (iSub-start_sub+1)/nSubRun;
                try
                    waitbar(progress,w,[num2str(progress*100,'%.f') '%, ' ...
                        num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' ...
                        num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
                end
            end
            fprintf('LSS: subject %d/%d completed (%.1f min).\n',iSub,end_sub,toc(start_time)/60);
        end
    % ====================== Subjects: parallel ===========================    
    else
        parfor iSub = start_sub:end_sub
            sub_check(iSub) = tmfc_LSS_one_subject(tmfc,iSub,false);
            if showGUI && ~isempty(D)
                try, send(D,[]); end
            end
        end
        if showGUI
            try
                main_GUI = guidata(findobj('Tag','TMFC_GUI'));
                set(main_GUI.TMFC_GUI_S6,'String', ...
                    [num2str(end_sub) '/' num2str(nSub) ' done'], ...
                    'ForegroundColor',[0.219 0.341 0.137]);
            end
        end
        fprintf('LSS: %d subjects completed (%.1f min).\n',nSubRun,toc(start_time)/60);
    end
catch ME
    if showGUI
        try, freeze_GUI(0); catch, end
    end
    fprintf('\n%s\n  TMFC LSS FAILED\n%s\n',repmat('=',1,30),repmat('=',1,30));
    rethrow(ME);
end

% Close waitbar
if showGUI
    try, delete(w); end
end

function unfreeze_after_ctrl_c()
    try
        delete(findall(0,'type','figure','Tag','tmfc_waitbar'));
        GUI = guidata(findobj('Tag','TMFC_GUI'));
        set([GUI.TMFC_GUI_B1, GUI.TMFC_GUI_B2, GUI.TMFC_GUI_B3, GUI.TMFC_GUI_B4,...
            GUI.TMFC_GUI_B5a, GUI.TMFC_GUI_B5b, GUI.TMFC_GUI_B6, GUI.TMFC_GUI_B7,...
            GUI.TMFC_GUI_B8, GUI.TMFC_GUI_B9, GUI.TMFC_GUI_B10, GUI.TMFC_GUI_B11,...
            GUI.TMFC_GUI_B12a,GUI.TMFC_GUI_B12b,GUI.TMFC_GUI_B13a,GUI.TMFC_GUI_B13b,...
            GUI.TMFC_GUI_B14a, GUI.TMFC_GUI_B14b], 'Enable','on');
    end
end
end

%% ========================================================================

%--------------------------------------------------------------------------
% LSS for one subject
%--------------------------------------------------------------------------
function out = tmfc_LSS_one_subject(tmfc,iSub,useParTrial)
    SPM = load(tmfc.subjects(iSub).path).SPM;
    conds = tmfc.LSS.conditions;

    if strcmp(tmfc.session_mode,'repeated')
        sess_list = 1:numel(SPM.Sess);
    else
        sess_list = unique([conds.sess]);
    end

    subdir = fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name);
    % Prepare subject output folders
    if exist(subdir,'dir'), rmdir(subdir,'s'); end
    mkdir(fullfile(subdir,'Betas'));
    mkdir(fullfile(subdir,'GLM_batches'));
    out.session = struct([]);

    % Check if SPM.mat has concatenated sessions 
    % (if spm_fmri_concatenate.m script was used)
    isConcat = numel(SPM.nscan) ~= numel(SPM.Sess);

    % Check rWLS
    isRWLS = strcmpi(SPM.xVi.form,'wls');

    % Loop through sessions
    for sidx = 1:numel(sess_list)
        s = sess_list(sidx);
        [trials,others] = tmfc_LSS_events(SPM,conds,s,tmfc.session_mode);
        nTrials = numel(trials.onset);
        if nTrials == 0, continue, end
        if numel(trials.onset) + numel(others.onset) < 2
            error('tmfc_LSS:Trials','At least two events are required to construct Other_trials.');
        end
        parallel_note = '';
        if useParTrial, parallel_note = ' (parallel trials)'; end
        fprintf('LSS: subject %d, session %d, %d trials%s.\n',iSub,s,nTrials,parallel_note);

        % Shared input file list and nuisance regressors for the full run
        scan_names = tmfc_LSS_scan_names(SPM,s);
        common = tmfc_LSS_spec_common(SPM,s,scan_names);

        % First trial is used to calculate SPM grand-mean scaling, thresholds
        % and image headers exactly once, via the shared helper
        firstjob = tmfc_LSS_trial_job(common,trials,others,1,subdir,s);
        SPM_LSS = tmfc_fmri_spec(firstjob);
        if isConcat, SPM_LSS = tmfc_fmri_concatenate(SPM_LSS,SPM.nscan); end
        fmriData = tmfc_LSS_cache(SPM_LSS,tmfc.defaults.maxmem);

        % Precompute the ESS baseline once per session
        if ~isRWLS
            fmriData = tmfc_LSS_prepare_rssB(fmriData,SPM_LSS);
        end
        
        % Loop through trials of interest
        condition = struct([]);
        ok = zeros(1,nTrials);
        if useParTrial
            parfor k = 1:nTrials
                tmfc_LSS_one_trial(common,trials,others,k,subdir,s, ...
                    SPM.nscan,isConcat,isRWLS,fmriData,tmfc.defaults);
                ok(k) = 1;
            end
        else
            for k = 1:nTrials
                tmfc_LSS_one_trial(common,trials,others,k,subdir,s, ...
                    SPM.nscan,isConcat,isRWLS,fmriData,tmfc.defaults);
                ok(k) = 1;
            end
        end
        for k = 1:nTrials
            c = trials.cond(k); j = trials.number(k);
            condition(c).trials(j) = ok(k);
        end
        out.session(s).condition = condition;
    end
end

%--------------------------------------------------------------------------
% Get trials of interest and trials of no interest
%--------------------------------------------------------------------------
function [trials,other] = tmfc_LSS_events(SPM,conds,s,mode)
    trials = struct('onset',[],'duration',[],'cond',[],'number',[],'names',{{}});
    chosen = [];
    for c = 1:numel(conds)
        if strcmp(mode,'repeated') || conds(c).sess == s
            ix = conds(c).number;
            if ix < 1 || ix > numel(SPM.Sess(s).U)
                error('tmfc_LSS:Condition','Condition %d missing in session %d.',ix,s);
            end
            U = SPM.Sess(s).U(ix);
            ons = U.ons(:);
            dur = U.dur(:);
            if numel(dur) == 1 && numel(ons) > 1, dur = repmat(dur,numel(ons),1); end
            if numel(dur) ~= numel(ons), error('tmfc_LSS:Duration','Incorrect duration count.'); end
            trials.onset = [trials.onset;ons]; 
            trials.duration = [trials.duration;dur];
            trials.cond = [trials.cond;repmat(ix,numel(ons),1)]; 
            trials.number = [trials.number;(1:numel(ons))'];
            trials.names = [trials.names;repmat(U.name(1),numel(ons),1)];
            chosen = [chosen ix]; 
        end
    end
    notChosen = setdiff(1:numel(SPM.Sess(s).U),chosen);
    other.onset = [];
    other.duration = [];
    for ix = notChosen
        ons = SPM.Sess(s).U(ix).ons(:);
        dur = SPM.Sess(s).U(ix).dur(:);
        if numel(dur) == 1 && numel(ons) > 1, dur = repmat(dur,numel(ons),1); end
        if numel(dur) ~= numel(ons), error('tmfc_LSS:Duration','Incorrect duration count.'); end
        other.onset = [other.onset;ons]; 
        other.duration = [other.duration;dur]; 
    end
end

%--------------------------------------------------------------------------
% Functional image filenames for one run
%--------------------------------------------------------------------------
function scans = tmfc_LSS_scan_names(SPM,s)
    row = SPM.Sess(s).row;
    scans = cell(numel(row),1);
    for i=1:numel(row)
        V = SPM.xY.VY(row(i));
        scans{i} = sprintf('%s,%d',V.fname,V.n(1));
    end
end

%--------------------------------------------------------------------------
% Prepare common LSS GLM specification
%--------------------------------------------------------------------------
function common = tmfc_LSS_spec_common(SPM,s,scans)
    common = struct();
    common.timing.units = SPM.xBF.UNITS;
    common.timing.RT = SPM.xY.RT;
    common.timing.fmri_t = SPM.xBF.T;
    common.timing.fmri_t0 = SPM.xBF.T0;
    common.sess.scans = scans;
    common.sess.regress = struct('name',{},'val',{});
    for i=1:numel(SPM.Sess(s).C.name)
        common.sess.regress(i).name = SPM.Sess(s).C.name{i};
        common.sess.regress(i).val = SPM.Sess(s).C.C(:,i);
    end
    common.sess.multi = {''};
    common.sess.multi_reg = {''};
    common.sess.hpf = SPM.xX.K(s).HParam;
    common.fact = struct('name',{},'levels',{});
    common.bases.hrf.derivs = [0 0];
    common.volt = 1;
    common.global = SPM.xGX.iGXcalc;
    common.mthresh = SPM.xM.gMT;
    if isfield(SPM.xM,'VM') && ~isempty(SPM.xM.VM)
        common.mask = {SPM.xM.VM.fname};
    else
        common.mask = {''};
    end
    switch lower(SPM.xVi.form)
        case {'i.i.d','none'}, common.cvi = 'None';
        case 'fast',           common.cvi = 'FAST';
        otherwise,             common.cvi = 'AR(1)';
    end
end

%--------------------------------------------------------------------------
% Specify LSS GLM batch for one trial
%--------------------------------------------------------------------------
function job = tmfc_LSS_trial_job(common,trials,other,k,subdir,s)
    job = common;
    job.dir = {fullfile(subdir,sprintf('LSS_Run_%d_Trial_%d',s,k))};
    job.sess.cond(1) = tmfc_LSS_event_cond('Current_trial',trials.onset(k),trials.duration(k));
    ix = [1:k-1 k+1:numel(trials.onset)];
    ons = [trials.onset(ix);other.onset];
    dur = [trials.duration(ix);other.duration];
    job.sess.cond(2) = tmfc_LSS_event_cond('Other_trials',ons,dur);
end

%--------------------------------------------------------------------------
% Define an LSS condition
%--------------------------------------------------------------------------
function c = tmfc_LSS_event_cond(name,ons,dur)
    c = struct('name',name,'onset',ons,'duration',dur,'tmod',0, ...
        'pmod',struct('name',{},'param',{},'poly',{}),'orth',1);
end

%--------------------------------------------------------------------------
% Specify and estimate LSS GLM for one trial
%--------------------------------------------------------------------------
function tmfc_LSS_one_trial(common,trials,other,k,subdir,s,concatN,isConcat,isRWLS,fmriData,defaults)
    job = tmfc_LSS_trial_job(common,trials,other,k,subdir,s);

    % Use in-memory SPM specification, with no data mapping or global scaling after the first trial of a run
    designJob = rmfield(job,'mask');
    designJob.sess.nscan = numel(job.sess.scans);
    designJob.sess = rmfield(designJob.sess,'scans');
    SPM_LSS = tmfc_fmri_spec(designJob);
    if isConcat, SPM_LSS = tmfc_fmri_concatenate(SPM_LSS,concatN); end
    if isRWLS
        SPM_LSS.xVi.form = 'wls';
        n = size(SPM_LSS.xX.X,1);
        SPM_LSS.xVi.Vi = cell(1,n);
        for j=1:n
            SPM_LSS.xVi.Vi{j} = sparse(j,j,1,n,n);
        end
        if isfield(SPM_LSS.xVi,'V'), SPM_LSS.xVi = rmfield(SPM_LSS.xVi,'V'); end
    end

    name = regexprep(char(trials.names{k}),' ','_');
    stem = sprintf('[Run_%d]_[Cond_%d]_[%s]_[Trial_%d]',s,trials.cond(k),name,trials.number(k));
    betaPath = fullfile(subdir,'Betas',['Beta_' stem '.nii']);
    batchPath = fullfile(subdir,'GLM_batches',['GLM_' stem '.mat']);

    % Estimate in-memory; publish a completed trial only after success
    tmpBeta = fullfile(subdir,'Betas',['_tmp_' stem '.nii']);
    try
        tmfc_LSS_estimate(SPM_LSS,fmriData,defaults,isRWLS,tmpBeta);
        [ok,msg] = movefile(tmpBeta,betaPath,'f');
        if ~ok, error('tmfc_LSS:Publish','Failed to rename beta: %s',msg); end
        matlabbatch = {struct('spm',struct('stats',struct('fmri_spec',job)))};
        save(batchPath,'matlabbatch');
    catch ME
        if exist(tmpBeta,'file'), delete(tmpBeta); end
        rethrow(ME);
    end
end

%--------------------------------------------------------------------------
% Cache functional data for LSS GLM estimation
%--------------------------------------------------------------------------
function fmriData = tmfc_LSS_cache(SPM,maxmem)
    % Read and mask one LSS run, avoiding repeated SPM global estimation.
    % When run does not fit a conservative maxmem budget, cache voxel indices
    % only, then stream the same voxel columns during trial estimation.
    
    V = SPM.xY.VY;
    if spm_mesh_detect(V)
        error('tmfc_LSS:Surface','This LSS implementation supports NIfTI volumes only.');
    end
    for j=1:numel(V)
        if ~spm_existfile(V(j).fname)
            error('tmfc_LSS:MissingFile','File not found: %s',V(j).fname);
        end
        % Same scaling propagation as tmfc_gPPI_estimate / spm_spm
        V(j).private.dat.scl_slope = V(j).pinfo(1);
        V(j).private.dat.scl_inter = V(j).pinfo(2);
    end
    spm_check_orientations(V);
    fmriData.V = V;
    fmriData.dim = V(1).dim;
    fmriData.mat = V(1).mat;
    n = numel(V);
    if numel(SPM.xM.TH) ~= n
        error('tmfc_LSS:Scans','Mask threshold vector does not match scan count.');
    end
    fmriData.nscan = n;
    fmriData.K = SPM.xX.K;
    fmriData.TH = SPM.xM.TH;
    fmriData.implicit = SPM.xM.I;
    fmriData.nanrep = spm_type(V(1).dt(1),'nanrep');
    
    % Explicit mask interpolation, identical to SPM estimation helpers
    mask = true(fmriData.dim);
    for m=1:numel(SPM.xM.VM)
        VM = SPM.xM.VM(m);
        C = spm_bsplinc(VM,[0 0 0 0 0 0]');
        v = true(fmriData.dim);
        [x1,x2] = ndgrid(1:fmriData.dim(1),1:fmriData.dim(2));
        M2 = inv(fmriData.mat\VM.mat);
        for z=1:fmriData.dim(3)
            y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*z+M2(1,4));
            y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*z+M2(2,4));
            y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*z+M2(3,4));
            v(:,:,z) = spm_bsplins(C,y1,y2,y3,[0 0 0 0 0 0]') > 0;
        end
        mask = mask & v;
    end
    
    % Keep workspace for simultaneous design, residuals and ReML matrices
    chunkVox = floor(maxmem/(8*n*8));
    if chunkVox < 1
        error('tmfc_LSS:Memory','maxmem too small even for a single scan-series voxel.');
    end
    chunkVox = min(chunkVox,prod(fmriData.dim));
    N = prod(fmriData.dim);
    nChunks = ceil(N/chunkVox);
    fmriData.blocks = repmat(struct('indices',[],'range',[],'Y',[],'KY',[]),nChunks,1);
    fmriData.inmask = 0;
    cachedBytes = 0;
    cacheLimit = 0.5*maxmem;
    keepCache = true;
    for b=1:nChunks
        ix = (b-1)*chunkVox+1:min(b*chunkVox,N);
        cmask = mask(ix);
        Y = zeros(n,numel(ix));
        for j=1:n
            if ~any(cmask), break; end
            Y(j,cmask) = spm_data_read(V(j),ix(cmask));
            cmask(cmask) = Y(j,cmask) > fmriData.TH(j);
            if fmriData.implicit && ~fmriData.nanrep && fmriData.TH(j) < 0
                cmask(cmask) = abs(Y(j,cmask)) > eps;
            end
        end
        cmask(cmask) = any(diff(Y(:,cmask),1));
        fmriData.blocks(b).indices = ix(cmask);
        fmriData.blocks(b).range = ix;
        fmriData.inmask = fmriData.inmask + nnz(cmask);
        if keepCache
            cachedBytes = cachedBytes + 16*n*nnz(cmask);
            if cachedBytes <= cacheLimit
                fmriData.blocks(b).Y = Y(:,cmask);
                fmriData.blocks(b).KY = spm_filter(fmriData.K,fmriData.blocks(b).Y);
            else
                % Switch the entire run to streaming to keep RAM bounded
                keepCache = false;
                for q=1:b
                    fmriData.blocks(q).Y = [];
                    fmriData.blocks(q).KY = [];
                end
            end
        end
    end
    if fmriData.inmask == 0
        error('tmfc_LSS:EmptyMask','Please check your data: there are no in-mask voxels.');
    end
    fmriData.cached = keepCache;
    if ~keepCache
        fprintf('LSS: memory budget exceeded; using streaming chunks (slower I/O).\n');
    end
end

%--------------------------------------------------------------------------
% Precompute the baseline residual sum of squares for ESS
%--------------------------------------------------------------------------
function D = tmfc_LSS_prepare_rssB(D,SPM)
    xX = SPM.xX;
    
    % Only the default non-rWLS effects-of-interest F-contrast
    if isfield(SPM.xVi,'V') || isfield(SPM.xVi,'Fcontrast')
        return;
    end   
    iX0 = [xX.iB xX.iG];
    if isempty(iX0), return; end
    
    X0 = xX.X(:,iX0);
    sB = spm_sp('Set',spm_filter(xX.K,X0));   
    rssB = cell(numel(D.blocks),1);
    
    for b = 1:numel(D.blocks)   
        ind = D.blocks(b).indices;
        if isempty(ind)
            rssB{b} = zeros(1,0);
            continue;
        end
        if D.cached
            KY = D.blocks(b).KY;
        else
            Y = zeros(D.nscan,numel(ind));
            for t = 1:D.nscan
                Y(t,:) = spm_data_read(D.V(t),ind);
            end
            KY = spm_filter(xX.K,Y);
        end    
        rssB{b} = sum(spm_sp('r',sB,KY).^2,1);
    end    
    D.essBaseline = struct('iX0',iX0,'X0',X0,'rssB',{rssB});
end

%--------------------------------------------------------------------------
% Waitbar for parallel subject mode
%--------------------------------------------------------------------------
function tmfc_LSS_parfor_waitbar(waitbarHandle,iterations)
    persistent w nSub start_time count_sub
    if nargin == 2
        w = waitbarHandle;
        nSub = iterations;
        start_time = tic;
        count_sub = 0;
    else
        count_sub = count_sub+1;
        if ~isempty(w) && ishandle(w)
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            time_remaining = (nSub-count_sub)*time_per_sub;
            hms = fix(mod(time_remaining,[0 3600 60])./[3600 60 1]);
            try
                waitbar(count_sub/nSub,w,[num2str(count_sub/nSub*100,'%.f') '%, ' ...
                    num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' ...
                    num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            end
        end
    end
end
