function [sys_info] = equiv_dict_stats(Dbar,W, sys_info)

if nargin < 3
 sys_info = struct('Mu',[],'Mu_av',[],'inf_norm',[],'cond_WDbar',[],'norm_WDbar',[],'norm_W',[],'cond_W',[]);
end


WDbar = W*Dbar;
A = mat_normalize(WDbar);
G = A'*A;
[sys_info.Mu , sys_info.Mu_av] = coherence(A);
sys_info.inf_norm = [sys_info.inf_norm , norm(A(:),'inf')];
sys_info.cond_WDbar = [sys_info.cond_WDbar, cond(WDbar)];
sys_info.norm_WDbar = [sys_info.norm_WDbar, norm(eye - G,'fro')];
sys_info.norm_W = [sys_info.norm_W, norm(W'*W - eye,'fro')];
sys_info.cond_W = [sys_info.cond_W, cond(W)];



end

