function [gp_pooled,gp_warm] = tmfc_gPPI_pooled_options(defaults,cid)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Parse and validate the pooled-Cy / warm-start options shared by
% tmfc_gPPI / tmfc_gPPI_FIR (multi-subject and one_sub flows). Both
% options are explicit opt-in; the default is the exact (byte-identical)
% path. Semantics follow evaluation/gppi/OPTIONS_DESIGN.md ("Revision
% 2026-11-16", R4/R5 and original section 2).
%
%   tmfc.defaults.gppi_method           'exact' (default) | 'pooled'
%   tmfc.defaults.gppi_reml_warm_start  logical scalar, default false
%
% Validation is strict (no silent fallbacks):
%   - invalid values            -> error <cid>:InvalidMethod /
%                                  <cid>:InvalidWarmStart
%   - pooled + analysis 1/3     -> error <cid>:PooledUnsupported
%   - exact + warm set          -> warning <cid>:PooledOptionsIgnored
% (rejected by the caller after the per-subject parent checks, alongside
%  rWLS/concat: pooled + FAST parent -> error <cid>:PooledFastUnsupported
%  (numerically pathological at real scale, supported noise: AR(1));
%  pooled + cvi 'None' a2 fails with the same legacy error as exact/main,
%  R5.)
%
% FORMAT [gp_pooled,gp_warm] = tmfc_gPPI_pooled_options(defaults,cid)
%
%   defaults  - tmfc.defaults struct (fields optional)
%   cid       - error/warning identifier prefix ('tmfc_gPPI' or
%               'tmfc_gPPI_FIR')
%   gp_pooled - true when gppi_method == 'pooled'
%   gp_warm   - true when a warm ReML start is requested AND pooled
%
% License: GPL-3.0-or-later. Derived from the TMFC toolbox (R. Masharipov).

method = 'exact';
if isfield(defaults,'gppi_method')
    method = defaults.gppi_method;
    if ~ischar(method) || ~isrow(method) || ~any(strcmp(method,{'exact','pooled'}))
        error([cid ':InvalidMethod'], ...
            'tmfc.defaults.gppi_method must be ''exact'' or ''pooled''');
    end
end

warm = false;
warm_set = false;
if isfield(defaults,'gppi_reml_warm_start')
    warm = defaults.gppi_reml_warm_start;
    warm_set = true;
    if ~islogical(warm) || numel(warm) ~= 1
        error([cid ':InvalidWarmStart'], ...
            'tmfc.defaults.gppi_reml_warm_start must be a logical scalar');
    end
end

gp_pooled = strcmp(method,'pooled');
gp_warm = warm;

if gp_pooled && defaults.analysis ~= 2
    error([cid ':PooledUnsupported'], ...
        ['pooled applies to ROI-to-ROI (analysis 2); analysis %d noise ' ...
         'operator is estimated inside SPM (spm_jobman) and cannot use ' ...
         'pooled Cy without replacing SPM estimation'],defaults.analysis);
end

if ~gp_pooled && warm_set && warm
    warning([cid ':PooledOptionsIgnored'], ...
        'gppi_reml_warm_start is ignored while gppi_method is ''exact''.');
end
end
