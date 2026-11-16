function [h0,keep] = tmfc_LSS_pooled_parent(SPM,SPM_LSS,sess_num)
% Map parent ReML hyperparameters and retained FAST components to one run.
% Replaces default SPM ReML initialisation for optional pooled warm starts.
% License: GPL-3.0-or-later.
if ~isfield(SPM.xVi,'h')
    error('tmfc_LSS_pooled_parent: parent SPM has no xVi.h (unestimated parent?)');
end
q  = 1:size(SPM_LSS.xX.X,1);                   % trial design: one run
qp = SPM.Sess(sess_num).row;                   % rows of the run in the parent
ti = comps_in_block(SPM_LSS.xVi.Vi,q);
pi = comps_in_block(SPM.xVi.Vi,qp);
if numel(ti) ~= numel(pi)
    error('tmfc_LSS_pooled_parent: parent has %d components in run %d, trial design %d', ...
        numel(pi),sess_num,numel(ti));
end
for j = 1:numel(ti)
    if ~isequal(full(SPM_LSS.xVi.Vi{ti(j)}(q,q)),full(SPM.xVi.Vi{pi(j)}(qp,qp)))
        error('tmfc_LSS_pooled_parent: parent component %d differs from trial component %d', ...
            pi(j),ti(j));
    end
end
h0  = zeros(numel(SPM_LSS.xVi.Vi),1);
h0(ti) = SPM.xVi.h(pi);
keep = h0 ~= 0;
end

function c = comps_in_block(Vi,q)
c = [];
for j = 1:numel(Vi)
    if nnz(Vi{j}(q,q)), c(end+1) = j; end %#ok<AGROW>
end
end
