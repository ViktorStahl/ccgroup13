%clear
%close all

%format long

% addpath(genpath('./')); %adds all the functions from subfolders to the path
% mfiles=getfilenames('./','BSeuCallU*.m')

%warning off
%pkg load dataframe
%echo off
%Problems = {'1','2','3','4'}
 % Problem 1 a) I
function bo = benchop(problem,method)

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
bo = tBSeuCallUI(index);    


elseif (problem == 2)
	% % Problem 1 b) I
	display('Problem 1 b) I');
	rootpath=pwd;
	S=[90,100,110]; K=100; T=1.0; r=0.03; sig=0.15;
	U=[10.726486710094511 4.820608184813253 1.828207584020458];

	filepathsBSamPutUI=getfilenames('./','BSamPutUI_*.m');
	par={S,K,T,r,sig};
	[timeBSamPutUI,relerrBSamPutUI] = executor(rootpath,filepathsBSamPutUI,U,par)

	tBSamPutUI=NaN(numel(Methods),1); rBSamPutUI=NaN(numel(Methods),1);
	for ii=1:numel(Methods)
	    for jj=1:numel(filepathsBSamPutUI)
	        a=filepathsBSamPutUI{jj}(3:3+numel(Methods{ii}));
            b=[method,'/']
	        %b=[Methods{ii},'/'];
	        if strcmp(a,b)
	            tBSamPutUI(ii)=timeBSamPutUI(jj);
	            rBSamPutUI(ii)=relerrBSamPutUI(jj);
	        end
	    end
	end
	
	cd(rootpath);
	display(m)
	for el=1:numel(Methods)
        if strcmp((Methods{el}),method)
            display('Result for the given problem and method is : ');
            index = el;
        end
    end 
bo = tBSamPutUI(index);
elseif (problem == 3)
	% % % Problem 1 c) I

	display('Problem 1 c) I');
	rootpath=pwd;
	S=[90,100,110]; K=100; T=1.0; r=0.03; sig=0.15; B=1.25*K;
	U=[1.822512255945242 3.294086516281595 3.221591131246868];

	filepathsBSupoutCallI=getfilenames('./','BSupoutCallI_*.m');
	par={S,K,T,r,sig,B};
	[timeBSupoutCallI,relerrBSupoutCallI] = executor(rootpath,filepathsBSupoutCallI,U,par)

	tBSupoutCallI=NaN(numel(Methods),1); rBSupoutCallI=NaN(numel(Methods),1);
	for ii=1:numel(Methods)
	    for jj=1:numel(filepathsBSupoutCallI)
	        a=filepathsBSupoutCallI{jj}(3:3+numel(Methods{ii}));
            b=[method,'/']
	        %b=[Methods{ii},'/'];
	        if strcmp(a,b)
	            tBSupoutCallI(ii)=timeBSupoutCallI(jj);
	            rBSupoutCallI(ii)=relerrBSupoutCallI(jj);
	        end
	    end
	end

	cd(rootpath)
		
	display(m)
	for el=1:numel(Methods)
        if strcmp((Methods{el}),method)
            display('Result for the given problem and method is : ');
            index = el;
        end
    end
bo = tBSupoutCallI(index) 
	
%  skippin Problem 1 a)  II because of octave incompatibility
elseif (problem == 4)
	 % Problem 1 b) II

	display('Problem 1 b) II');
	rootpath=pwd;
	S=[97,98,99]; K=100; T=0.25; r=0.1; sig=0.01;
	U=[3.000000000000682 2.000000000010786   1.000000000010715];

	filepathsBSamPutUII=getfilenames('./','BSamPutUII_*.m');
	par={S,K,T,r,sig};
	[timeBSamPutUII,relerrBSamPutUII] = executor(rootpath,filepathsBSamPutUII,U,par);
	
	tBSamPutUII=NaN(numel(Methods),1); rBSamPutUII=NaN(numel(Methods),1);
	for ii=1:numel(Methods)
	    for jj=1:numel(filepathsBSamPutUII)
	        a=filepathsBSamPutUII{jj}(3:3+numel(Methods{ii}));
	        b=[method,'/']
            %b=[Methods{ii},'/'];
	        if strcmp(a,b)
	            tBSamPutUII(ii)=timeBSamPutUII(jj);
	            rBSamPutUII(ii)=relerrBSamPutUII(jj);
	        end
	    end
	end

	cd(rootpath);
 
	display(method)
	for el=1:numel(Methods)
        if strcmp((Methods{el}),method)
            display('Result for the given problem and method is : ');
            display(tBSamPutUII(el))
            index = el;
        end
    end
bo = tBSamPutUII(index);
 % skipping Problem 1 c) II because of octave incompatibility
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1:numel(Methods)
%	C = {"Methods","Problem1","Problem2","Problem3","Problem4";Methods{i},tBSeuCallUI(i),tBSamPutUI(i),tBSupoutCallI(i),tBSamPutUII(i)};	
%	dataframe(C)
%end

%display('tBSeucallUI')
%Table2 = textable(tBSeuCallUI,tBSamPutUI, tBSupoutCallI,tBSamPutUII)
%'RowNames',Methods}
%Table2=table(tBSeuCallUI,tBSamPutUI,tBSupoutCallI,tBSamPutUII,'RowNames',Methods)
%err=[rBSeuCallUI,rBSamPutUI,rBSupoutCallI,rBSamPutUII];
%err=round(log10(err));

% Now use this table as input in our input struct:
%input.data = Table2;
%input.error = err;

% Set the row format of the data values (in this example we want to use
% integers only):
%input.dataFormat = {'%.1e'};

% Switch transposing/pivoting your table:
%input.transposeTable = 1;

% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
%input.tableColumnAlignment = 'c';

% Switch table borders on/off:
%input.tableBorders = 0;

% Switch to generate a complete LaTex document or just a table:
%input.makeCompleteLatexDocument = 0;

%latex = latexTable(input);
