%Methods={'MC','MC-S','QMC-S','MLMC','MLMC-A',...    
%   'FFT','FGL','COS',...    
%   'FD','FD-NU','FD-AD',...    
%   'RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT'};
%problems={1,2,3,4}
%for i=1:numel(Methods) 
%	for j=numel(problems)
%		benchop(j,Methods{i})
%	end
%end
res = benchop(4, 'FD');
