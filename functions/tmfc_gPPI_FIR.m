function [sub_check,contrasts] = tmfc_gPPI_FIR(tmfc,ROI_set_number,start_sub,end_sub)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Estimates gPPI-FIR GLMs. Saves individual connectivity matrices
% (ROI-to-ROI analysis) and connectivity images (seed-to-voxel analysis)
% for each condition of interest.
%
% The difference between classic gPPI GLM and gPPI-FIR GLM is that the
% latter uses finite impulse response (FIR) functions (instead of canonical
% HRF function) to model activations for conditions of interest and
% conditions of no interest. The FIR model allows to model activations
% with arbitrary hemodynamic response shapes.
%
% Note: gPPI-FIR uses FIR basis functions for psychological regressors. 
% Therefore, temporal and dispersion derivative basis functions from the
% original HRF model are not used in the gPPI-FIR GLM.
% However, PPI regressors can still be HRF-, time-derivative-, 
% or dispersion-derivative-based, depending on the
% explicitly selected condition/basis-function pairs.
% If the original GLMs contain parametric or time modulators, they
% will be included in the gPPI-FIR GLMs. 
%
% FORMAT [sub_check,contrasts] = tmfc_gPPI_FIR(tmfc)
% Run a function starting from the first subject in the list.
%
%   tmfc.subjects.path     - Paths to individual SPM.mat files
%   tmfc.subjects.name     - Subject names within the TMFC project
%   tmfc.project_path      - Path where all results will be saved
%   tmfc.defaults.parallel - 0 or 1 (sequential or parallel computing)
%   tmfc.defaults.headless - true or false (disable GUI and waitbars)
%   tmfc.defaults.maxmem   - e.g. 2^31 = 2GB (how much RAM can be used)
%   tmfc.defaults.resmem   - true or false (store temporary files in RAM)
%   tmfc.defaults.ufp      - Global F threshold for ReML (default: 1)
%   tmfc.defaults.analysis - 1 (Seed-to-voxel and ROI-to-ROI analyses)
%                          - 2 (ROI-to-ROI analysis only)
%                          - 3 (Seed-to-voxel analysis only)
%
%   tmfc.session_mode      - 'strict' or 'repeated'
%                            In repeated mode, run-specific gPPI-FIR
%                            estimates are averaged across all runs
%
%   tmfc.ROI_set(ROI_set_number).gPPI_FIR.window - FIR window length (in seconds)
%   tmfc.ROI_set(ROI_set_number).gPPI_FIR.bins   - Number of FIR time bins
%
%   tmfc.ROI_set                  - List of selected ROIs
%   tmfc.ROI_set.type             - Type of the ROI set
%   tmfc.ROI_set.set_name         - Name of the ROI set
%   tmfc.ROI_set.ROIs.name        - Name of the selected ROI
%   tmfc.ROI_set.ROIs.path_masked - Paths to ROI images masked by
%                                   individual first-level masks
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions        - List of conditions of interest
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.sess   - Run number (SPM session, as specified in SPM.Sess)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.number - Condition number (as specified in SPM.Sess.U)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.pmod   - Parametric/Time modulator number (see SPM.Sess.U.P)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.bf     - Basis function:
%                                                         1 = canonical HRF
%                                                         2 = time derivative
%                                                         3 = dispersion derivative
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.name   - Condition name (as specified in SPM.Sess.U.name(kPmod))
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.file_name - Condition-specific file names:
%
%   Canonical HRF file name:
%
%   ['[Run_' num2str(iSess) ']_[Cond_' num2str(jCond) ']_[' ...
%   regexprep(char(SPM.Sess(iSess).U(jCond).name(kPmod)),' ','_') ']']
%
%   Derivative file name:
%
%   ['[Run_' num2str(iSess) ']_[Cond_' num2str(jCond) ']_[' ...
%   regexprep(char(SPM.Sess(iSess).U(jCond).name(kPmod)),' ','_') ']_[' ...
%   bf_file{kBF} ']']
%
%   where bf_file{kBF} is:
%   'TimeDeriv' - time derivative
%   'DispDeriv' - dispersion derivative
%
% Run number and condition number must match the original SPM.mat file.
% Consider, for example, a task design with two runs (SPM sessions). Both runs 
% contain three task regressors for "Cond A", "Cond B" and "Errors". If
% you are only interested in comparing "Cond A" and "Cond B", the following
% structure must be specified (see tmfc_conditions_GUI and tmfc_get_conditions):
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).sess   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).number = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).bf     = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).name = 'Cond_A';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).file_name = '[Run_1]_[Cond_1]_[Cond_A]';
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).sess   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).bf     = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).name = 'Cond_B';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).file_name = '[Run_1]_[Cond_2]_[Cond_B]';
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).number = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).bf     = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).name = 'Cond_A';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).file_name = '[Run_2]_[Cond_1]_[Cond_A]';
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).bf     = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).name = 'Cond_B';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).file_name = '[Run_2]_[Cond_2]_[Cond_B]';
%
% If you also want to include time and dispersion derivatives for Cond_B
% in run 2, add:
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).bf     = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).name = 'Cond_B';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).file_name = '[Run_2]_[Cond_2]_[Cond_B]_[TimeDeriv]';
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).bf     = 3;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).name = 'Cond_B';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).file_name = '[Run_2]_[Cond_2]_[Cond_B]_[DispDeriv]';
%
% Example of the ROI set (see tmfc_select_ROIs_GUI and tmfc_create_ROI_set):
%
%   tmfc.ROI_set(1).set_name = 'two_ROIs';
%   tmfc.ROI_set(1).type = 'binary_images';
%   tmfc.ROI_set(1).ROIs(1).name = 'ROI_1';
%   tmfc.ROI_set(1).ROIs(2).name = 'ROI_2';
%   tmfc.ROI_set(1).ROIs(1).path_masked(1).subjects = '...\Subject_0001\ROI_1_masked.nii';
%   tmfc.ROI_set(1).ROIs(1).path_masked(2).subjects = '...\Subject_0002\ROI_1_masked.nii';
%
% FORMAT [sub_check,contrasts] = tmfc_gPPI_FIR(tmfc,ROI_set_number,start_sub,end_sub)
% Process an inclusive subject range (e.g. for independent HPC jobs).
%
%   tmfc           - As above
%   ROI_set_number - Number of the ROI set in the tmfc structure
%   start_sub      - First subject to process (default: 1)
%   end_sub        - Last subject to process (default: final subject)
%
% =========================================================================
% Copyright (C) 2026 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if nargin == 1
   ROI_set_number = 1;
   start_sub = 1;
elseif nargin == 2
   start_sub = 1;
end

% Headless mode (default: preserve GUI behavior)
if ~isfield(tmfc.defaults,'headless')
    tmfc.defaults.headless = false;
end
showGUI = ~tmfc.defaults.headless && usejava('desktop');

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

nSub = length(tmfc.subjects);
if nargin < 4 || isempty(end_sub)
    end_sub = nSub;
end
validateattributes(start_sub,{'numeric'},{'scalar','integer','>=',1,'<=',nSub});
validateattributes(end_sub,{'numeric'},{'scalar','integer','>=',start_sub,'<=',nSub});
nSubRun = end_sub - start_sub + 1;
nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);
cond_list = tmfc.ROI_set(ROI_set_number).gPPI.conditions;
nCond = length(cond_list);

sub_check = zeros(1,nSub);
if nargin < 4 && start_sub > 1
    sub_check(1:start_sub) = 1; % Preserve legacy behavior
end

% Prepare gPPI FIR folders
if tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 2
    if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI'))
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI','asymmetrical'));
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI','symmetrical'));
    end
end

if tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 3
    for iROI = 1:nROI
        if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','Seed_to_voxel',tmfc.ROI_set(ROI_set_number).ROIs(iROI).name))
            mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','Seed_to_voxel',tmfc.ROI_set(ROI_set_number).ROIs(iROI).name));
        end
    end
end

for iROI = 1:nROI
    if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','GLM_batches',tmfc.ROI_set(ROI_set_number).ROIs(iROI).name))
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','GLM_batches',tmfc.ROI_set(ROI_set_number).ROIs(iROI).name));
    end
end

% Initialize SPM
spm('defaults','fmri');
spm_get_defaults('cmdline',true);

% Initialize parfor
useParROI = false;
useParSub = false;
hasPCT = (exist('parfor','builtin')==5) && license('test','Distrib_Computing_Toolbox');
if tmfc.defaults.parallel==1 && hasPCT
    p = gcp('nocreate');
    if isempty(p)
        parpool;
        p = gcp('nocreate');
    end
    nW = p.NumWorkers;
    % Single ROI, many subjects - parallelize over subjects
    % Many ROIs - parallelize over ROIs
    if (nSubRun >= nW) && (nROI < nW)
        useParSub = true;
    else
        useParROI  = (nROI >= nW);
    end
end

% Initialize GUI progress reporting only in interactive mode
D = [];
if showGUI
    try, figure(findobj('Tag','TMFC_GUI')); end
    cleanupObj = onCleanup(@unfreeze_after_ctrl_c);

    if ~useParSub
        w = waitbar(0,'Please wait...','Name','gPPI-FIR GLM estimation','Tag', 'tmfc_waitbar');
        start_time = tic;
        count_sub = 1;
    else
        try % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;
            w = waitbar(0,'Please wait...','Name','gPPI-FIR GLM estimation','Tag','tmfc_waitbar');
            afterEach(D, @tmfc_parfor_waitbar);     % Command to update waitbar
            tmfc_parfor_waitbar(w,nSubRun,1);
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
if ~useParSub
    % ===================== Subjects: sequential ==========================  
    for iSub = start_sub:end_sub
        sub_check(iSub) = tmfc_gPPI_FIR_one_sub(tmfc,ROI_set_number,iSub,cond_list,nCond,nROI,useParROI);
        
        % Update waitbar
        if showGUI
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            count_sub = count_sub + 1;
            time_remaining = (end_sub-iSub)*time_per_sub;
            hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
            progress = (iSub-start_sub+1)/nSubRun;
            try
                waitbar(progress, w, [num2str(progress*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            end
        end
    
    end
else
    % ====================== Subjects: parallel ===========================
    parfor iSub = start_sub:end_sub
        sub_check(iSub) = tmfc_gPPI_FIR_one_sub(tmfc,ROI_set_number,iSub,cond_list,nCond,nROI,false);
        if showGUI
            try, send(D,[]); end
        end
    end
end

% Default contrasts info
for cond_PPI = 1:nCond
    if strcmp(session_mode,'repeated')
        contrasts(cond_PPI).title = ...
            ['[All_Runs]_[Cond_' num2str(cond_list(cond_PPI).number) ']_[' ...
            regexprep(cond_list(cond_PPI).name,' ','_') ']'];
        if isfield(cond_list,'bf') && ~isempty(cond_list(cond_PPI).bf)
            if cond_list(cond_PPI).bf == 2
                contrasts(cond_PPI).title = [contrasts(cond_PPI).title '_[TimeDeriv]'];
            elseif cond_list(cond_PPI).bf == 3
                contrasts(cond_PPI).title = [contrasts(cond_PPI).title '_[DispDeriv]'];
            end
        end
    else
        contrasts(cond_PPI).title = cond_list(cond_PPI).file_name;
    end
    contrasts(cond_PPI).weights = zeros(1,nCond);
    contrasts(cond_PPI).weights(1,cond_PPI) = 1;
end

% Close waitbar
if showGUI
    try
        delete(w);
    end
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

%--------------------------------------------------------------------------
% gPPI-FIR for one subject
%--------------------------------------------------------------------------
function ok = tmfc_gPPI_FIR_one_sub(tmfc,ROI_set_number,iSub,cond_list,nCond,nROI,useParROI)

    ok = 0;

    if ~isfield(tmfc,'session_mode')
        session_mode = 'strict';
    else
        session_mode = tmfc.session_mode;
    end

    %=====================[ Specify gPPI-FIR GLM ]=========================
    SPM = load(tmfc.subjects(iSub).path).SPM;

    % Runs and PPI regressors for this subject
    % ---------------------------------------------------------------------
    if strcmp(session_mode,'repeated')

        sess_num = 1:length(SPM.Sess);
        nSess = length(sess_num);

        cond_in_sess = true(nSess,nCond);

    elseif strcmp(session_mode,'strict')

        sess_num = unique([cond_list.sess]);
        nSess = length(sess_num);

        cond_in_sess = false(nSess,nCond);

        for cond_PPI = 1:nCond
            kSess = find(sess_num == cond_list(cond_PPI).sess,1);
            cond_in_sess(kSess,cond_PPI) = true;
        end

    else
        error('session_mode must be ''strict'' or ''repeated''.');
    end

    PPI_num = zeros(nSess,nCond);

    for kSess = 1:nSess
        jConds = find(cond_in_sess(kSess,:));
        PPI_num(kSess,jConds) = 1:length(jConds);
    end

    % Check if SPM.mat has concatenated sessions 
    % (if spm_fmri_concatenate.m script was used)
    if size(SPM.nscan,2) == size(SPM.Sess,2)
        SPM_concat = 0;
    else
        SPM_concat = 1;
    end
    concat_scans = SPM.nscan;

    % Check for rWLS
    if strcmp(SPM.xVi.form,'wls')
        rWLS = 1;
    else
        rWLS = 0;
    end

    % Prepare common GLM specification
    base_spec = struct();
    base_spec.timing.units  = SPM.xBF.UNITS;
    base_spec.timing.RT     = SPM.xY.RT;
    base_spec.timing.fmri_t = SPM.xBF.T;
    base_spec.timing.fmri_t0 = SPM.xBF.T0;
    for kSess = 1:nSess
        % Functional images
        if SPM_concat == 0
            for image = 1:SPM.nscan(sess_num(kSess))
                base_spec.sess(kSess).scans{image,1} = [SPM.xY.VY(SPM.Sess(sess_num(kSess)).row(image)).fname ',' ...
                    num2str(SPM.xY.VY(SPM.Sess(sess_num(kSess)).row(image)).n(1))];
            end
        else
            for image = 1:size(SPM.xY.VY,1)
                base_spec.sess(kSess).scans{image,1} = [SPM.xY.VY(SPM.Sess(kSess).row(image)).fname ',' ...
                    num2str(SPM.xY.VY(SPM.Sess(kSess).row(image)).n(1))];
            end
        end
        % Conditions (including PSY regressors)
        for cond = 1:length(SPM.Sess(sess_num(kSess)).U)
            base_spec.sess(kSess).cond(cond).name = SPM.Sess(sess_num(kSess)).U(cond).name{1};
            base_spec.sess(kSess).cond(cond).onset = SPM.Sess(sess_num(kSess)).U(cond).ons;
            base_spec.sess(kSess).cond(cond).duration = SPM.Sess(sess_num(kSess)).U(cond).dur;
            base_spec.sess(kSess).cond(cond).tmod = 0;
            if length(SPM.Sess(sess_num(kSess)).U(cond).name)>1
                for PM_number = 1:length(SPM.Sess(sess_num(kSess)).U(cond).P)
                    base_spec.sess(kSess).cond(cond).pmod(PM_number).name = SPM.Sess(sess_num(kSess)).U(cond).P(PM_number).name;
                    base_spec.sess(kSess).cond(cond).pmod(PM_number).param = SPM.Sess(sess_num(kSess)).U(cond).P(PM_number).P;
                    base_spec.sess(kSess).cond(cond).pmod(PM_number).poly = SPM.Sess(sess_num(kSess)).U(cond).P(PM_number).h;
                end
            else
                base_spec.sess(kSess).cond(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
            end
            base_spec.sess(kSess).cond(cond).orth = SPM.Sess(sess_num(kSess)).U(cond).orth;
        end
        base_spec.sess(kSess).multi = {''};
        base_spec.sess(kSess).multi_reg = {''};
        base_spec.sess(kSess).hpf = SPM.xX.K(sess_num(kSess)).HParam;
    end
    base_spec.fact = struct('name', {}, 'levels', {});
    base_spec.bases.fir.length = tmfc.ROI_set(ROI_set_number).gPPI_FIR.window;
    base_spec.bases.fir.order  = tmfc.ROI_set(ROI_set_number).gPPI_FIR.bins;
    base_spec.volt = 1;
    base_spec.global = SPM.xGX.iGXcalc;
    base_spec.mthresh = SPM.xM.gMT;
    try
        base_spec.mask = {SPM.xM.VM.fname};
    catch
        base_spec.mask = {''};
    end
    if strcmp(SPM.xVi.form,'i.i.d') || strcmp(SPM.xVi.form,'none')
        base_spec.cvi = 'None';
    elseif strcmp(SPM.xVi.form,'fast') || strcmp(SPM.xVi.form,'FAST')
        base_spec.cvi = 'FAST';
    else
        base_spec.cvi = 'AR(1)';
    end

    % Load VOI time series
    VOI_Yraw = cell(nSess,nROI);
    need_Y = tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 2;
    Y = [];
    for kSess = 1:nSess
        KWY = [];
        for jROI = 1:nROI
            VOI = fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs', ...
                tmfc.subjects(iSub).name,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(jROI).name '_' num2str(sess_num(kSess)) '.mat']);
            tmp = load(VOI);
            if isfield(tmp,'Yraw')
                VOI_Yraw{kSess,jROI} = tmp.Yraw;
            else
                VOI_Yraw{kSess,jROI} = tmp.Y;     % Backward compatibility: old VOI files contain only Y
            end
            if need_Y
                KWY(:,jROI) = tmp.Y;              % Filtered and whitened data
            end
            clear VOI tmp
        end
        if need_Y
            Y = [Y; KWY];
        end
        clear KWY
    end

    % Specify gPPI-FIR GLM batches
    for jROI = 1:nROI
        matlabbatch{1}.spm.stats.fmri_spec = base_spec;
        matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(tmfc.project_path,'ROI_sets',...
            tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR',...
            tmfc.subjects(iSub).name,tmfc.ROI_set(ROI_set_number).ROIs(jROI).name)};

        % Loop through sessions
        for kSess = 1:nSess
   
            % Add PPI regressors
            jConds = find(cond_in_sess(kSess,:));
            for ii = 1:length(jConds)
                cond_PPI = jConds(ii);

                % Repeated runs
                if strcmp(session_mode,'repeated')
                    cond = cond_list(cond_PPI).number;

                    if isfield(cond_list,'pmod') && ~isempty(cond_list(cond_PPI).pmod)
                        pmod = cond_list(cond_PPI).pmod;
                    else
                        pmod = 1;
                    end

                    if isfield(cond_list,'bf') && ~isempty(cond_list(cond_PPI).bf)
                        bf = cond_list(cond_PPI).bf;
                    else
                        bf = 1;
                    end

                    cond_name = char(SPM.Sess(sess_num(kSess)).U(cond).name(pmod));

                    file_name = ['[Run_' num2str(sess_num(kSess)) ']_[Cond_' num2str(cond) ']_[' regexprep(cond_name,' ','_') ']'];

                    if bf == 2
                        file_name = [file_name '_[TimeDeriv]'];
                    elseif bf == 3
                        file_name = [file_name '_[DispDeriv]'];
                    end

                % Strict mode
                else
                    file_name = strrep(cond_list(cond_PPI).file_name,'[Sess_','[Run_');
                end

                file_new = fullfile(tmfc.project_path,'ROI_sets', ...
                    tmfc.ROI_set(ROI_set_number).set_name,'PPIs', ...
                    tmfc.subjects(iSub).name, ...
                    ['PPI_[' regexprep(tmfc.ROI_set(ROI_set_number).ROIs(jROI).name,' ','_') ']_' file_name '.mat']);

                % Backward compatibility
                file_old = strrep(file_new,'[Run_','[Sess_');

                if exist(file_new,'file')
                    PPI = load(file_new);
                elseif exist(file_old,'file')
                    PPI = load(file_old);
                else
                    error('PPI file not found: %s',file_new);
                end

                PPI_id = PPI_num(kSess,cond_PPI);
                matlabbatch{1}.spm.stats.fmri_spec.sess(kSess).regress(PPI_id).name = ['PPI_' PPI.PPI.name];
                matlabbatch{1}.spm.stats.fmri_spec.sess(kSess).regress(PPI_id).val = PPI.PPI.ppi;
                clear PPI
            end

            % Add PHYS regressors
            nPPI_sess = sum(cond_in_sess(kSess,:));
            matlabbatch{1}.spm.stats.fmri_spec.sess(kSess).regress(nPPI_sess+1).name = ['Seed_' tmfc.ROI_set(ROI_set_number).ROIs(jROI).name];
            matlabbatch{1}.spm.stats.fmri_spec.sess(kSess).regress(nPPI_sess+1).val = VOI_Yraw{kSess,jROI};  % Preloaded VOI time series 
            
            % Confounds       
            for conf = 1:length(SPM.Sess(sess_num(kSess)).C.name)
                matlabbatch{1}.spm.stats.fmri_spec.sess(kSess).regress(conf+nPPI_sess+1).name = SPM.Sess(sess_num(kSess)).C.name{1,conf};
                matlabbatch{1}.spm.stats.fmri_spec.sess(kSess).regress(conf+nPPI_sess+1).val = SPM.Sess(sess_num(kSess)).C.C(:,conf);
            end
        end   

        batch{jROI} = matlabbatch;
        save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR', ...
            'GLM_batches',tmfc.ROI_set(ROI_set_number).ROIs(jROI).name,...
            [tmfc.subjects(iSub).name '_gPPI_FIR_GLM.mat']),'matlabbatch');
        clear matlabbatch PPI   
    end
    
    % Specify all ROI GLMs, reuse common design across ROIs
    if SPM_concat == 1
        designs = tmfc_gPPI_FIR_design(batch,concat_scans);
    else
        designs = tmfc_gPPI_FIR_design(batch);
    end
    clear batch

    %=====================[ Estimate gPPI-FIR GLM ]========================
    tmfc_gPPI_FIR_run_estimation(tmfc,ROI_set_number,iSub,designs, ...
        PPI_num,cond_in_sess,sess_num,cond_list,Y,rWLS,useParROI);

    ok = 1;
end

%--------------------------------------------------------------------------
% Specify gPPI-FIR GLM designs for all ROIs
%--------------------------------------------------------------------------
function designs = tmfc_gPPI_FIR_design(batch,concat_scans)
    if nargin < 2, concat_scans = []; end
    designs = cell(size(batch));
    if isempty(batch), return; end
    first_job = batch{1}{1}.spm.stats.fmri_spec;
    template = tmfc_fmri_spec(first_job);
    if ~isempty(concat_scans)
        template = tmfc_fmri_concatenate(template,concat_scans);
    end
    for jROI = 1:numel(batch)
        job = batch{jROI}{1}.spm.stats.fmri_spec;
        designs{jROI} = reuse_design(template,job);
        designs{jROI}.swd = job.dir{1};
    end
    end
    
    function SPM = reuse_design(template,job)
    % Reuse HRF/task design, filtering, covariance, headers and global scaling
    % Only the ROI-dependent nuisance (PPI, PHYS) columns are replaced
    SPM = template;
    if numel(job.sess) ~= numel(SPM.Sess)
        error('tmfc_gPPI_FIR_reuse_design:SessionCount','Session count changed across ROIs.');
    end
    for s = 1:numel(job.sess)
        regress = job.sess(s).regress;
        nR = numel(regress);
        if nR ~= numel(SPM.Sess(s).C.name)
            error('tmfc_gPPI_FIR_reuse_design:RegressorCount', ...
                'Number of user regressors changed for session %d.',s);
        end
        nRows = numel(SPM.Sess(s).row);
        C = [];
        names = cell(1,nR);
        for k = 1:nR
            if numel(regress(k).val) ~= nRows
                error('Length of regressor is not commensurate with data points.');
            end
            C = [C, regress(k).val(:)]; 
            names{k} = regress(k).name;
        end
        cols = SPM.Sess(s).col;
        cCols = cols(end-nR+1:end);
        SPM.Sess(s).C.C = C;
        SPM.Sess(s).C.name = names;
        % Mean-center user regressors, as in spm_fMRI_design
        if nR > 0
            SPM.xX.X(SPM.Sess(s).row,cCols) = spm_detrend(C);
            for k = 1:nR
                SPM.xX.name{cCols(k)} = [sprintf('Sn(%i) ',s) names{k}];
            end
        end
    end
end

%--------------------------------------------------------------------------
% Run gPPI-FIR GLM estimation for all ROIs
%--------------------------------------------------------------------------
function tmfc_gPPI_FIR_run_estimation(tmfc,ROI_set_number,iSub,designs,PPI_num,cond_in_sess,sess_num,cond_list,Y,rWLS,useParROI)

    nROI = numel(designs);
    need_Y = tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 2;
    need_images = tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 3;
    roi_beta = cell(1,nROI);
    
    % Collect final NIfTI paths; estimation writes to temporary _tmp.nii files
    image_paths = cell(1,nROI);
    if need_images
        for jROI = 1:nROI
            [~,beta_paths,~,mean_paths] = tmfc_gPPI_FIR_get_beta_output_info( ...
                tmfc,ROI_set_number,iSub,jROI,PPI_num,cond_in_sess, ...
                sess_num,cond_list,designs{jROI});
            image_paths{jROI} = [beta_paths mean_paths(~cellfun('isempty',mean_paths))];
        end
    end
    
    % Remove previous completion indicators to prevent false GUI progress
    tmfc_gPPI_FIR_invalidate_progress(tmfc,ROI_set_number,iSub,numel(cond_list));
    
    % Run estimation
    if useParROI
        parfor jROI = 1:nROI
            roi_beta{jROI} = tmfc_gPPI_FIR_one_roi(tmfc,ROI_set_number,iSub,jROI, ...
                designs{jROI},PPI_num,cond_in_sess,sess_num,cond_list,Y,rWLS);
        end
    else
        for jROI = 1:nROI
            roi_beta{jROI} = tmfc_gPPI_FIR_one_roi(tmfc,ROI_set_number,iSub,jROI, ...
                designs{jROI},PPI_num,cond_in_sess,sess_num,cond_list,Y,rWLS);
        end
    end
    
    % Save gPPI-FIR matrices
    asym_paths = {};
    sym_paths = {};
    if need_Y
        beta = cat(3,roi_beta{:});
        [asym_paths,sym_paths] = tmfc_save_gPPI_FIR_matrices(tmfc,ROI_set_number,iSub,beta, ...
            PPI_num,cond_in_sess,sess_num,cond_list,designs{1});
    end
    
    % Publish completed results, keeping GUI completion indicators last
    tmfc_gPPI_FIR_publish_files(asym_paths);
    tmfc_gPPI_FIR_publish_files(sym_paths);
    if need_images
        for jROI = 1:nROI-1
            tmfc_gPPI_FIR_publish_files(image_paths{jROI});
        end
        tmfc_gPPI_FIR_publish_files(image_paths{nROI});
    end
end

%--------------------------------------------------------------------------
% Estimate gPPI-FIR GLM for a single ROI
%--------------------------------------------------------------------------
function roi_beta = tmfc_gPPI_FIR_one_roi(tmfc,ROI_set_number,iSub,jROI,SPM,PPI_num,cond_in_sess,sess_num,cond_list,Y,isRWLS)

    roi_beta = [];

    % ROI-to-ROI-only 
    %----------------------------------------------------------------------
    if tmfc.defaults.analysis == 2
        xX = SPM.xX;
        if isfield(SPM.xX,'W')
            SPM.xX = rmfield(SPM.xX,'W');
        end
        if isfield(SPM.xVi,'V')
            SPM.xVi = rmfield(SPM.xVi,'V');
        end
        if ~isRWLS
            xVi = tmfc_spm_est_non_sphericity(SPM,tmfc);
        else
            xVi = tmfc_spm_rwls_est_non_sphericity(SPM,tmfc);
        end
        W = spm_sqrtm(spm_inv(xVi.V));
        W = W.*(abs(W) > 1e-6);
        xKXs = spm_sp('Set',spm_filter(xX.K,W*xX.X));
        xKXs.X = full(xKXs.X);
        pKX = spm_sp('x-',xKXs);

    % Seed-to-voxel  
    %----------------------------------------------------------------------
    else
        [beta_idx,beta_paths,run_paths,mean_paths,mean_names] = ...
            tmfc_gPPI_FIR_get_beta_output_info(tmfc,ROI_set_number,iSub,jROI, ...
                PPI_num,cond_in_sess,sess_num,cond_list,SPM);
        % Use temporary NIfTI filenames until all ROIs have been processed
        for k = 1:numel(beta_paths)
            beta_paths{k} = tmfc_gPPI_FIR_tmp_path(beta_paths{k});
        end
        for c = 1:numel(mean_paths)
            if isempty(mean_paths{c}), continue; end
            for kSess = 1:numel(run_paths{c})
                run_paths{c}{kSess} = tmfc_gPPI_FIR_tmp_path(run_paths{c}{kSess});
            end
            mean_paths{c} = tmfc_gPPI_FIR_tmp_path(mean_paths{c});
        end
        % Estimate
        if isRWLS
            SPM.xVi.form = 'wls';
            nScan = sum(SPM.nscan);
            for jScan = 1:nScan
                SPM.xVi.Vi{jScan} = sparse(nScan,nScan);
                SPM.xVi.Vi{jScan}(jScan,jScan) = 1;
            end
            pKX = tmfc_gPPI_estimate_rwls(SPM,beta_idx,beta_paths,tmfc);
        else
            pKX = tmfc_gPPI_estimate(SPM,beta_idx,beta_paths,tmfc);
        end
        % Mean over runs (session_mode = 'repeated')
        for c = 1:numel(mean_paths)
            if isempty(mean_paths{c}), continue; end
            V = spm_vol(run_paths{c}{1});
            image = spm_read_vols(V);
            for kSess = 2:numel(run_paths{c})
                image = image + spm_read_vols(spm_vol(run_paths{c}{kSess}));
            end
            image = image/numel(run_paths{c});
            V.fname = mean_paths{c};
            V.descrip = ['gPPI-FIR beta: ' mean_names{c}];
            spm_write_vol(V,image);
        end
    end

    % ROI-to-ROI
    %----------------------------------------------------------------------
    if tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 2
        roi_beta = pKX*Y;
    end
end

%--------------------------------------------------------------------------
% Define gPPI beta image indices and output filenames
%--------------------------------------------------------------------------
function [beta_idx,beta_paths,run_paths,mean_paths,mean_names] = ...
    tmfc_gPPI_FIR_get_beta_output_info(tmfc,ROI_set_number,iSub,jROI,PPI_num,cond_in_sess,sess_num,cond_list,SPM)

    base = fullfile(tmfc.project_path,'ROI_sets', ...
        tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','Seed_to_voxel', ...
        tmfc.ROI_set(ROI_set_number).ROIs(jROI).name);
    subject_name = tmfc.subjects(iSub).name;
    nCond = numel(cond_list);
    beta_idx = [];
    beta_paths = {};
    run_paths = cell(1,nCond);
    mean_paths = cell(1,nCond);
    mean_names = cell(1,nCond);
    if isfield(tmfc,'session_mode')
        session_mode = tmfc.session_mode;
    else
        session_mode = 'strict';
    end
    for cond_PPI = 1:nCond
        if strcmp(session_mode,'strict')
            kSess = find(cond_in_sess(:,cond_PPI),1);
            idx = tmfc_gPPI_FIR_ppi_beta_idx(SPM,kSess,PPI_num(kSess,cond_PPI));
            fpath = fullfile(base,[subject_name '_Contrast_' num2str(cond_PPI,'%04.f') '_' cond_list(cond_PPI).file_name '.nii']);
            beta_idx(end+1) = idx;
            beta_paths{end+1} = fpath;
        elseif strcmp(session_mode,'repeated')
            if isfield(cond_list,'bf') && ~isempty(cond_list(cond_PPI).bf)
                bf = cond_list(cond_PPI).bf;
            else
                bf = 1;
            end
            run_paths{cond_PPI} = cell(1,numel(sess_num));
            for kSess = 1:numel(sess_num)
                idx = tmfc_gPPI_FIR_ppi_beta_idx(SPM,kSess,PPI_num(kSess,cond_PPI));
                file_name = ['[Run_' num2str(sess_num(kSess)) ']_[Cond_' ...
                    num2str(cond_list(cond_PPI).number) ']_[' ...
                    regexprep(cond_list(cond_PPI).name,' ','_') ']'];
                if bf == 2
                    file_name = [file_name '_[TimeDeriv]'];
                elseif bf == 3
                    file_name = [file_name '_[DispDeriv]'];
                end
                fpath = fullfile(base,[subject_name '_Contrast_' num2str(cond_PPI,'%04.f') '_' file_name '.nii']);
                beta_idx(end+1) = idx;
                beta_paths{end+1} = fpath;
                run_paths{cond_PPI}{kSess} = fpath;
            end
            mean_name = ['[All_Runs]_[Cond_' ...
                num2str(cond_list(cond_PPI).number) ']_[' ...
                regexprep(cond_list(cond_PPI).name,' ','_') ']'];
            if bf == 2
                mean_name = [mean_name '_[TimeDeriv]'];
            elseif bf == 3
                mean_name = [mean_name '_[DispDeriv]'];
            end
            mean_paths{cond_PPI} = fullfile(base,[subject_name '_Contrast_' num2str(cond_PPI,'%04.f') '_' mean_name '.nii']);
            mean_names{cond_PPI} = mean_name;
        else
            error('session_mode must be ''strict'' or ''repeated''.');
        end
    end
end

%--------------------------------------------------------------------------
% Resolve gPPI-FIR PPI regressor columns from the SPM design 
%--------------------------------------------------------------------------
function idx = tmfc_gPPI_FIR_ppi_beta_idx(SPM,kSess,ppi_id)
    % PPI regressors are the FIRST user-specified regressors in each
    % session. In spm_fMRI_design, user regressors are at the END of
    % SPM.Sess(kSess).col, after all FIR-expanded task regressors.
    % Reuse the same mapping as tmfc_gPPI_FIR_design/reuse_design.
    sess = SPM.Sess(kSess);
    nR = numel(sess.C.name);
    if ppi_id < 1 || ppi_id > nR || ppi_id ~= fix(ppi_id) || nR > numel(sess.col)
        error('tmfc_gPPI_FIR:InvalidPPIIndex', ...
            'Invalid PPI regressor %g in session %d (user regressors: %d).', ...
            ppi_id,kSess,nR);
    end
    local_name = sess.C.name{ppi_id};
    if ~ischar(local_name) || ~strncmp(local_name,'PPI_',4)
        error('tmfc_gPPI_FIR:InvalidPPIRegressor', ...
            'User regressor %d in session %d is not a PPI regressor.',ppi_id,kSess);
    end
    idx = sess.col(end-nR+ppi_id);
    expected_name = [sprintf('Sn(%i) ',kSess) local_name];
    if idx > numel(SPM.xX.name) || idx < 1 || ...
            ~strcmp(SPM.xX.name{idx},expected_name)
        error('tmfc_gPPI_FIR:PPIDesignMismatch', ...
            'PPI regressor mapping mismatch: session %d, regressor %d.', ...
            kSess,ppi_id);
    end
end

%--------------------------------------------------------------------------
% Save gPPI-FIR matrices
%--------------------------------------------------------------------------
function [asym_paths,sym_paths] = tmfc_save_gPPI_FIR_matrices(tmfc,ROI_set_number,iSub,beta,PPI_num,cond_in_sess,sess_num,cond_list,SPM)

    asym_paths = {};
    sym_paths = {};
    gPPI_path = fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR');

    nCond = length(cond_list);

    if ~isfield(tmfc,'session_mode')
        session_mode = 'strict';
    else
        session_mode = tmfc.session_mode;
    end

    for cond_PPI = 1:nCond

        % Strict mode
        % -----------------------------------------------------------------
        if strcmp(session_mode,'strict')

            kSess = find(cond_in_sess(:,cond_PPI),1);

            beta_number = tmfc_gPPI_FIR_ppi_beta_idx(SPM,kSess,PPI_num(kSess,cond_PPI));

            ppi_matrix = squeeze(beta(beta_number,:,:));
            ppi_matrix(1:size(ppi_matrix,1)+1:end) = nan;

            symm_ppi_matrix = (ppi_matrix + ppi_matrix')/2;

            output_path = fullfile(gPPI_path,'ROI_to_ROI','asymmetrical', ...
                [tmfc.subjects(iSub).name '_Contrast_' ...
                num2str(cond_PPI,'%04.f') '_' ...
                cond_list(cond_PPI).file_name '.mat']);
            save(tmfc_gPPI_FIR_tmp_path(output_path),'ppi_matrix');
            asym_paths{end+1} = output_path;

            output_path = fullfile(gPPI_path,'ROI_to_ROI','symmetrical', ...
                [tmfc.subjects(iSub).name '_Contrast_' ...
                num2str(cond_PPI,'%04.f') '_' ...
                cond_list(cond_PPI).file_name '.mat']);
            save(tmfc_gPPI_FIR_tmp_path(output_path),'symm_ppi_matrix');
            sym_paths{end+1} = output_path;

        % Repeated runs
        % -----------------------------------------------------------------
        elseif strcmp(session_mode,'repeated')

            nSess = length(sess_num);
            ppi_matrix_all = [];

            if isfield(cond_list,'bf') && ~isempty(cond_list(cond_PPI).bf)
                bf = cond_list(cond_PPI).bf;
            else
                bf = 1;
            end

            for kSess = 1:nSess

                beta_number = tmfc_gPPI_FIR_ppi_beta_idx(SPM,kSess,PPI_num(kSess,cond_PPI));

                ppi_matrix = squeeze(beta(beta_number,:,:));
                ppi_matrix(1:size(ppi_matrix,1)+1:end) = nan;

                symm_ppi_matrix = (ppi_matrix + ppi_matrix')/2;

                file_name = ...
                    ['[Run_' num2str(sess_num(kSess)) ']_[Cond_' ...
                    num2str(cond_list(cond_PPI).number) ']_[' ...
                    regexprep(cond_list(cond_PPI).name,' ','_') ']'];

                if bf == 2
                    file_name = [file_name '_[TimeDeriv]'];
                elseif bf == 3
                    file_name = [file_name '_[DispDeriv]'];
                end

                output_path = fullfile(gPPI_path,'ROI_to_ROI','asymmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' ...
                    num2str(cond_PPI,'%04.f') '_' file_name '.mat']);
                save(tmfc_gPPI_FIR_tmp_path(output_path),'ppi_matrix');
                asym_paths{end+1} = output_path;

                output_path = fullfile(gPPI_path,'ROI_to_ROI','symmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' ...
                    num2str(cond_PPI,'%04.f') '_' file_name '.mat']);
                save(tmfc_gPPI_FIR_tmp_path(output_path),'symm_ppi_matrix');
                sym_paths{end+1} = output_path;

                if kSess == 1
                    ppi_matrix_all = ppi_matrix;
                else
                    ppi_matrix_all = ppi_matrix_all + ppi_matrix;
                end
            end

            % Equal-weight mean across runs
            ppi_matrix = ppi_matrix_all/nSess;
            symm_ppi_matrix = (ppi_matrix + ppi_matrix')/2;

            file_name = ...
                ['[All_Runs]_[Cond_' num2str(cond_list(cond_PPI).number) ']_[' ...
                regexprep(cond_list(cond_PPI).name,' ','_') ']'];

            if bf == 2
                file_name = [file_name '_[TimeDeriv]'];
            elseif bf == 3
                file_name = [file_name '_[DispDeriv]'];
            end

            output_path = fullfile(gPPI_path,'ROI_to_ROI','asymmetrical', ...
                [tmfc.subjects(iSub).name '_Contrast_' ...
                num2str(cond_PPI,'%04.f') '_' file_name '.mat']);
            save(tmfc_gPPI_FIR_tmp_path(output_path),'ppi_matrix');
            asym_paths{end+1} = output_path;

            output_path = fullfile(gPPI_path,'ROI_to_ROI','symmetrical', ...
                [tmfc.subjects(iSub).name '_Contrast_' ...
                num2str(cond_PPI,'%04.f') '_' file_name '.mat']);
            save(tmfc_gPPI_FIR_tmp_path(output_path),'symm_ppi_matrix');
            sym_paths{end+1} = output_path;

            clear ppi_matrix_all

        else
            error('session_mode must be ''strict'' or ''repeated''.');
        end

        clear ppi_matrix symm_ppi_matrix
    end
end


%--------------------------------------------------------------------------
% Temporary file naming 
%--------------------------------------------------------------------------
function temp_path = tmfc_gPPI_FIR_tmp_path(final_path)
    [folder,name,ext] = fileparts(final_path);
    temp_path = fullfile(folder,[name '_tmp' ext]);
end

%--------------------------------------------------------------------------
% Remove previous gPPI-FIR files used by TMFC.m to check completion
% Preserve custom contrasts (indices greater than nCond)
%-------------------------------------------------------------------------
function tmfc_gPPI_FIR_invalidate_progress(tmfc,ROI_set_number,iSub,nCond)
    gPPI_path = fullfile(tmfc.project_path,'ROI_sets', ...
        tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR');
    if tmfc.defaults.analysis == 2
        folder = fullfile(gPPI_path,'ROI_to_ROI','symmetrical');
        ext = '.mat';
    else
        lastROI = tmfc.ROI_set(ROI_set_number).ROIs(end).name;
        folder = fullfile(gPPI_path,'Seed_to_voxel',lastROI);
        ext = '.nii';
    end
    if ~exist(folder,'dir'), return; end
    prefix = [tmfc.subjects(iSub).name '_Contrast_'];
    listing = dir(folder);
    for j = 1:numel(listing)
        filename = listing(j).name;
        if listing(j).isdir || ~strncmp(filename,prefix,numel(prefix))
            continue;
        end
        rest = filename(numel(prefix)+1:end);
        token = regexp(rest,'^(\d+)_','tokens','once');
        if isempty(token), continue; end
        cond_idx = str2double(token{1});
        if cond_idx < 1 || cond_idx > nCond, continue; end
        if numel(filename) < numel(ext) || ~strcmpi(filename(end-numel(ext)+1:end),ext)
            continue;
        end
        old_path = fullfile(folder,filename);
        delete(old_path);
        if exist(old_path,'file')
            error('tmfc_gPPI_FIR:InvalidateFailed','Cannot invalidate previous gPPI completion file: %s',old_path);
        end
    end
end

%--------------------------------------------------------------------------
% Publish completed gPPI-FIR output files
%--------------------------------------------------------------------------
function tmfc_gPPI_FIR_publish_files(final_paths)
    for j = 1:numel(final_paths)
        final_path = final_paths{j};
        temp_path = tmfc_gPPI_FIR_tmp_path(final_path);
        if exist(temp_path,'file') ~= 2
            error('tmfc_gPPI_FIR:MissingTemporaryFile','Completed gPPI output is missing: %s',temp_path);
        end
        [ok,msg] = movefile(temp_path,final_path,'f');
        if ~ok
            error('tmfc_gPPI_FIR:PublishFailed','Cannot publish gPPI output %s: %s',final_path,msg);
        end
    end
end

%--------------------------------------------------------------------------
% Waitbar for parallel mode
%--------------------------------------------------------------------------
function tmfc_parfor_waitbar(waitbarHandle,iterations,firstsub)
    persistent w nSub start_sub start_time count_sub 
    if nargin == 3
        w = waitbarHandle;
        nSub = iterations;
        start_sub = firstsub - 1;
        start_time = tic;
        count_sub = 1;
    else
        if isvalid(w)         
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            iSub = start_sub + count_sub;
            time_remaining = (nSub-iSub)*time_per_sub;
            hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
            waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            count_sub = count_sub + 1;
        end
    end
end
