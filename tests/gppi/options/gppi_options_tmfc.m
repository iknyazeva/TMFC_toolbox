function tmfc = gppi_options_tmfc(fx,proj,opts)

% tmfc config for the options tests (gPPI options tests): identical to the config
% run_gppi_exact_tests builds, plus the pooled options passed through
% the public tmfc.defaults interface (gppi_method / gppi_reml_warm_start;
% see evaluation/gppi/OPTIONS_IMPL.md). Shared by run_gppi_options_tests
% and run_gppi_options_negative. tests/gppi/exact is NOT modified.
%
% FORMAT tmfc = gppi_options_tmfc(fx,proj,opts)
%   fx   - fixture bundle from gppi_exact_fixture
%   proj - per-run project root (inputs copied there beforehand)
%   opts - struct; fields gppi_method / gppi_reml_warm_start are copied
%          into tmfc.defaults ONLY when present (absent = exact default)
%
% License: GPL-3.0-or-later.

tmfc = struct();
tmfc.project_path = proj;
tmfc.subjects = struct('path',fx.spmmat,'name',fx.subject);
tmfc.defaults = struct('parallel',0,'maxmem',fx.maxmem, ...
    'resmem',true,'analysis',2);
if isfield(opts,'gppi_method')
    tmfc.defaults.gppi_method = opts.gppi_method;
end
if isfield(opts,'gppi_reml_warm_start')
    tmfc.defaults.gppi_reml_warm_start = opts.gppi_reml_warm_start;
end
if isfield(opts,'analysis')
    tmfc.defaults.analysis = opts.analysis;
end
tmfc.session_mode = fx.session_mode;
tmfc.ROI_set = struct('set_name',fx.set_name,'type','binary_images');
for j = 1:numel(fx.roi_names)
    tmfc.ROI_set(1).ROIs(j).name = fx.roi_names{j};
    tmfc.ROI_set(1).ROIs(j).path_masked(1).subjects = ...
        fullfile(proj,'masks',['ROI_' regexprep(fx.roi_names{j},' ','_') '_masked.nii']);
end
tmfc.ROI_set(1).gPPI.conditions = fx.conditions;
if strcmp(fx.method,'fir')
    tmfc.ROI_set(1).gPPI_FIR.window = fx.window;
    tmfc.ROI_set(1).gPPI_FIR.bins   = fx.bins;
end
end
