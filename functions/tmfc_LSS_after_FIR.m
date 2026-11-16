function [sub_check] = tmfc_LSS_after_FIR(tmfc,start_sub)

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
% It uses residual time-series created after the FIR task regression.
% Residual time-series (Res_*.nii images) are free of (1) co-activations
% and (2) confounds specified in the original SPM.mat file (e.g., motion,
% physiological noise, etc). Using residual time-series, we control for
% spurious inflation of functional connectivity estimates due to
% co-activations.
% 
% This function does not add confounds (e.g., motion, 
% physiological noise, etc) to the LSS models since they have already been
% regressed out during the FIR task regression.
%
% This function does not apply temporal autocorrelation modeling (AR(1) or
% FAST models) since the residuals were already whitened during the
% FIR task regression.  
%
% This function does not apply high-pass filter (HPF) since the residuals
% have already been filtered during the FIR task regression. 
%
% This function uses the binary mask created during estimation of the
% standard GLM (specified in the SPM.mat input files) as an explicit mask.
%
% Note: If the original GLMs contain parametric or time modulators, they
% will be removed from the LSS GLMs. If the original GLMs contain time and
% dispersion derivatives, they will be removed from the LSS GLMs.
%
% FORMAT [sub_check] = tmfc_LSS_after_FIR(tmfc)
% Run the function starting from the first subject in the list.
%
%   tmfc.subjects.path     - Paths to individual SPM.mat files
%   tmfc.subjects.name     - Subject names within the TMFC project
%   tmfc.project_path      - Path where all results will be saved
%   tmfc.defaults.parallel - 0 or 1 (LSS runs serially; 1 warns)
%   tmfc.defaults.maxmem   - e.g. 2^31 = 2GB (how much RAM can be used at
%                            the same time during GLM estimation)
%   tmfc.defaults.resmem   - true or false (store temporary files during
%                            GLM estimation in RAM or on disk)
%
%   tmfc.session_mode      - 'strict' or 'repeated'
%                            'strict'   - use the selected run-specific conditions
%                            'repeated' - process the same conditions across all runs
%
%   tmfc.LSS_after_FIR.conditions        - List of conditions of interest
%   tmfc.LSS_after_FIR.conditions.sess   - Session number
%                                          (as specified in SPM.Sess)
%   tmfc.LSS_after_FIR.conditions.number - Condition number  
%                                          (as specified in SPM.Sess.U)
%   tmfc.LSS_after_FIR.conditions.name   - Condition name
%                                          (as specified in SPM.Sess.U.name)
%   tmfc.LSS_after_FIR.conditions.file_name - Condition-specific file names:
%   (['[Run_' num2str(iSess) ']_[Cond_' num2str(jCond) ']_[' ...
%    regexprep(char(SPM.Sess(iSess).U(jCond).name(1)),' ','_') ']'];)
%
% Run (session) number and condition number must match the original SPM.mat file.
% Consider, for example, a task design with two runs (SPM sessions). Both runs 
% contain three task regressors for "Cond A", "Cond B" and "Errors". If
% you are only interested in comparing "Cond A" and "Cond B", the following
% structure must be specified (see tmfc_conditions_GUI and tmfc_get_conditions):
%
%   tmfc.LSS_after_FIR.conditions(1).sess   = 1;   
%   tmfc.LSS_after_FIR.conditions(1).number = 1; 
%   tmfc.LSS_after_FIR.conditions(1).name = 'Cond_A'; 
%   tmfc.LSS_after_FIR.conditions(1).file_name = '[Run_1]_[Cond_1]_[Cond_A]';  
%
%   tmfc.LSS_after_FIR.conditions(2).sess   = 1;
%   tmfc.LSS_after_FIR.conditions(2).number = 2;
%   tmfc.LSS_after_FIR.conditions(2).name = 'Cond_B';
%   tmfc.LSS_after_FIR.conditions(2).file_name = '[Run_1]_[Cond_2]_[Cond_B]';  
%
%   tmfc.LSS_after_FIR.conditions(3).sess   = 2;
%   tmfc.LSS_after_FIR.conditions(3).number = 1;
%   tmfc.LSS_after_FIR.conditions(3).name = 'Cond_A';
%   tmfc.LSS_after_FIR.conditions(3).file_name = '[Run_2]_[Cond_1]_[Cond_A]';  
%
%   tmfc.LSS_after_FIR.conditions(4).sess   = 2;
%   tmfc.LSS_after_FIR.conditions(4).number = 2;
%   tmfc.LSS_after_FIR.conditions(4).name = 'Cond_B';
%   tmfc.LSS_after_FIR.conditions(4).file_name = '[Run_2]_[Cond_2]_[Cond_B]'; 
%
% FORMAT [sub_check] = tmfc_LSS_after_FIR(tmfc, start_sub)
% Run the function starting from a specific subject in the path list.
%
%   tmfc           - As above
%   start_sub      - Subject number in the list to start computations from
%
% Uses in-memory batch FWL with exact rank-deficient fallback, replacing
% job-based estimation in fixed reference; numerically identical within
% float32 output precision. The FIR residuals are not re-scaled.
%
% =========================================================================
% Copyright (C) 2026 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if nargin < 2
    start_sub = 1;
end
sub_check = tmfc_LSS_after_FIR_fwl_batch(tmfc,start_sub);
end
