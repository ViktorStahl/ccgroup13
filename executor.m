% Copyright (c) 2015, BENCHOP, Slobodan Milovanović
% All rights reserved.
% This MATLAB code has been written for the BENCHOP project.
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%    * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%    * BENCHOP article is properly cited by the user of the BENCHOP codes when publishing/reporting related scientific results.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
function [time,relerr] = executor(rootpath,filepaths,U,par)
%Executing the code and doing performance analysis.
%   2014, Slobodan Milovanovic
M=3;
N=length(filepaths);
codes=cell(N,1);
time=zeros(N,1);
relerr=zeros(N,1);
for ii=1:N
    [pathstr,codes{ii}] = fileparts(filepaths{ii});
    cd(pathstr);
    display(filepaths(ii));
    try
        for jj=1:M+1
            tic
            R=feval(codes{ii},par{:}); %Time is measured here!
            tim=toc;
            t(jj)=tim;
        end
        relerr(ii)=max(abs(R-U)./abs(U));
        if relerr(ii)<=10^-4
            time(ii)=sum(t(2:M+1))/M;
        elseif relerr(ii)<=10^-3
            time(ii)=sum(t(2:M+1))/M+999000000;
        else
            time(ii)=sum(t(2:M+1))/M+888000000;
        end

    catch exception1
        time(ii)=NaN;
        relerr(ii)=NaN;
        cd(rootpath);
        continue
    end

    cd(rootpath);
end

end
