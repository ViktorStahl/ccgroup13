function test = benchop(problem,method)

Methods={'MC','MC-S','QMC-S','MLMC','MLMC-A',...
    'FFT','FGL','COS',...
    'FD','FD-NU','FD-AD',...
    'RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT'};
if (problem == 1)
        display('Problem 1 a) I');
        display(method)
        rootpath=pwd;
        S=[90,100,110]; K=100; T=1.0; r=0.03; sig=0.15;
        U=[2.758443856146076 7.485087593912603 14.702019669720769];

        filepathsBSeuCallUI=getfilenames('./','BSeuCallUI_*.m');
        par={S,K,T,r,sig};
        [timeBSeuCallUI,relerrBSeuCallUI] = executor(rootpath,filepathsBSeuCallUI,U,par);

        tBSeuCallUI=NaN(numel(Methods),1); rBSeuCallUI=tBSeuCallUI;
        for ii=1:numel(Methods)
                for jj=1:numel(filepathsBSeuCallUI)
                        a=filepathsBSeuCallUI{jj}(3:3+numel(Methods{ii}));
            b=[method,'/'];
                        %b=[Methods{ii},'/'];
                        if strcmp(a,b)
                                tBSeuCallUI(ii)=timeBSeuCallUI(jj);
                                rBSeuCallUI(ii)=relerrBSeuCallUI(jj);
                        end
                end
        end

        cd(rootpath);
        for el=1:numel(Methods)
                if strcmp((Methods{el}),method)
                        display('Result for the given problem and method is : ');
                        %display(tBSeuCallUI(el));
            index = el;

                end
    end
test = tBSeuCallUI(index);
end
