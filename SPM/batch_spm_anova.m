function [X]=batch_spm_anova(S);

% A general function for N-way mixed (within+between subjects) ANOVAs in SPM5/SPM8
% (though assumes same number of conditions per group)
%                                                                R Henson Oct 2006
%
% The only required argument in S is:
%    imgfiles    - cell array of cell arrays of image filenames for each group and subject
%
% Optional arguments in S are:
%    outdir      - output directory for SPM analysis files
%    mask        - analysis mask, either an image filename or a float for proportion of global (or [] for neither)
%    nsph_flag   - whether nonsphericity correction should be applied
%    sub_effects - whether subject effects included (yes by default)
%    contrasts   - cell array of contrast structures, with fields c
%                  (matrix), type ('F' or 'T') and name (optional)
%    user_regs   - additional user-specified regressors, either:
%                       1) the same number per group, ie a cell-array of regressors for each group, with one value for each condition and subject in that group
%                       2) a single set of matrix that spans all groups, ie a matrix with one row per scan 
%    uUFp        - uncorrected p-value threshold for mask for nonsphericity
%
% The only complicated bit is organising imgfiles correctly, so here's an example:
%
%  imgfiles{1}{1} = ['mydir/grp1_sub1_con1.nii'; 'mydir/grp1_sub1_con2.nii'];
%  imgfiles{1}{2} = ['mydir/grp1_sub2_con1.nii'; 'mydir/grp1_sub2_con2.nii'];
%  imgfiles{2}{1} = ['mydir/grp2_sub1_con1.nii'; 'mydir/grp2_sub1_con2.nii'];
%  imgfiles{2}{2} = ['mydir/grp2_sub2_con1.nii'; 'mydir/grp2_sub2_con2.nii'];
%
% Actually, the organisation of the later-added user_regs is also a bit complex, so here is an example of group-specific regressors:
%
%  user_regs{1} = [[1:4]' rand(4,1)];       % 2 regressors for group 1 
%  user_regs{2} = [[1:4]' rand(4,1)];       % 2 regressors for group 2
%
% where each group has 2 subjects with 2 conditions (ie 4 values), continuing above example
% (currently must be same number of regressors per group), whereas here is an example of 
% user-specified regressors than span all groups:
%
%  user_regs = [[1:8]' rand(8,1)];          % 2 regressors across all groups
%
% Updated 11-3-14 to allow user-specified regressors that span all groups
% Updated 18-6-13 to handle masking properly (note previous S.maskimg is now S.mask)

try 
    imgfiles = S.imgfiles;
catch
    error('Must provide image filenames in S.imgfiles')
end

try 
    outdir = S.outdir;
catch
    outdir = pwd;
end

try 
    nsph_flag = S.nsph_flag;
catch
    nsph_flag = -9;  % Nonsphericity to be determined by design (below)
end

try 
    sub_effects = S.sub_effects;
catch
    sub_effects = 1;
end

try 
    contrasts = S.contrasts;
catch
    contrasts = [];
end

try
    user_regs = S.user_regs;    % user_regs must either be a cell array {grp}{con} of Num_Subs x Num_Regs, or a matrix of Num_Scans X Num_Regs 
catch
    user_regs = cell(1,length(imgfiles));
end

try
    uUFp = S.uUFp;              
catch
    uUFp = 0.001;
end
try
    spm_get_defaults(['stats.' lower(spm_get_defaults('modality')) '.ufp'],uUFp)  % SPM8
catch
    spm_get_defaults(['stats.' lower(spm('CheckModality')) '.ufp'],uUFp) % SPM12
end

ngrp = length(imgfiles);

%% Data files
P={}; cname={};
np=0; nc=0;
for g=1:ngrp
    nsub(g) = length(imgfiles{g});
    ncon    = size(imgfiles{g}{1},1);   % Assumes same ncon for all subjects
    for n=1:ncon
        for s=1:nsub(g)
            np=np+1;
     	    P{np} = imgfiles{g}{s}(n,:);
        end
        nc=nc+1;
        cname{nc} = sprintf('grp %d, con %d, mean',g,n);
        if iscell(user_regs)
            for r=1:length(user_regs{g})
                nc=nc+1;
                cname{nc} = sprintf('grp %d, con %d, reg %d',g,n,r);
            end
        end
    end
end
if ~iscell(user_regs)
    for r=1:size(user_regs,2)
        nc=nc+1;
        cname{nc} = sprintf('glb reg %d',r);
    end
end

try pflag = S.pflag
catch pflag = 0;   %whether you want figures of design matrix/nonsphericity before estimating
end

totsub = sum(nsub);           % (total number of subjects)
nscan  = sum(ncon.*nsub);     % number of rows in X

try eval(sprintf('!mkdir %s',outdir)); end
cd(outdir);

if sub_effects
  for s=1:totsub                  % add subject effects
    cname{end+1} = sprintf('subject %d',s);
  end
end

%% Assemble SPM structure
SPM = [];
SPM.nscan = nscan;
SPM.xY.P  = P;
for i=1:SPM.nscan
    SPM.xY.VY(i) = spm_vol(SPM.xY.P{i});
end

%% Sort out any masking
try 
    mask = S.mask;
catch
    mask     = [];
    SPM.xM   = ones(nscan,1)*-Inf;
    sGXcalc  = 'omit';
end

if ~isempty(mask)
    if isfloat(mask)        % mean voxel value (no checking at moment!)
        glob=[];
        for i=1:SPM.nscan
            glob(i,1) = spm_global(SPM.xY.VY(i));
        end
        SPM.xM  = glob*mask;
        sGXcalc = sprintf('%3.2f of mean voxel value per image',mask);
    else                    % Assume a single, BINARY image
        try
            SPM.xM.VM(1) = spm_vol(mask);   
            SPM.xM.TH    = ones(nscan,1)*-Inf;
            SPM.xM.I     = 1;
            sGXcalc      = sprintf('From file %s',mask);
        catch
            error('Cannot open specified mask image');
        end
    end
end


%% Build design matrix (X), Indices (Ind) and NONSPHERICITY (vi) (inelegant, but gets there...!)
X=[]; Ind=[]; vi={};
nv=0; z=zeros(nscan,nscan); os=0;

for g=1:ngrp
    ns = nsub(g);
    nr = ncon*ns;
    id = [1:ns]';
    
    tmp = kron(eye(ncon),ones(ns,1));
    
    if iscell(user_regs)
        nreg = size(user_regs{g},2);
        
        if ~isempty(user_regs{g});
            tmp = [tmp user_regs{g}];   % user_regs must be (ns x ncon) by nreg
        end
    else
        nreg = 0;
    end
    
    tmpX = [zeros(nr,(ncon+nreg)*(g-1)),...  % assumes same number of user_regs per group
        tmp,...
        zeros(nr,(ncon+nreg)*(ngrp-g))];
    
    % could add constants for group effects if wish
    
    if ncon>1 & sub_effects
      tmpX = [tmpX zeros(size(tmpX,1),sum(nsub(1:(g-1)))) kron(ones(ncon,1),eye(ns))];
    end
    
    if g>1
      if ncon>1 & sub_effects
         X = [X zeros(size(X,1),ns); tmpX];
      else
	     X = [X; tmpX];
      end
    else
      X = tmpX;
    end
    
    % Indices for effects (unnecessary really?)
    Ind = [Ind; ones(nr,1),...                  %kron(ones(ncon,1),id),...
                kron([1:ncon]',ones(ns,1)),...
                kron(ones(ncon,1),id),...      
                ones(nr,1)*g];
    
 % Nonsphericity
    
    if nsph_flag ~= 0      % ie, unless user turns off...

    % unequal covariances (within conditions; independent between groups)
     if nsph_flag==1 | (ncon>2 & ns>1)
      nsph_flag = 1;
      for c1 = 1:ncon
        for c2 = (c1+1):ncon
            nv = nv+1;
            v = z;
            v( os + (c1-1)*ns + id, os + (c2-1)*ns + id )=eye(ns);
            v( os + (c2-1)*ns + id, os + (c1-1)*ns + id )=eye(ns);
            vi{nv} = sparse(v);   
        end
      end
     end

    % unequal variances (need if unequal covariances)
%     if ngrp>1 & all(nsub>1)   % This won't work (need unequal vars too)
     if nsph_flag==1 | (ngrp>1 & ns>1)
      nsph_flag = 1;
      for c1 = 1:ncon
        nv = nv+1;
        v = z;
        v(os + (c1-1)*ns + id, os + (c1-1)*ns + id)=eye(ns);
        vi{nv} = sparse(v);
      end
     end
    
      os = os + ncon*ns;
    end
end

if ~iscell(user_regs)
    if size(user_regs,1) == size(X,1)
        X = [X user_regs];
    else
        error('User-specified regressor does not have enough rows');
    end
end

    
if nsph_flag<0; nsph_flag=0; end

%% If want to peek
if pflag
 figure,imagesc(X),colormap('gray')
 figure,hold on,colormap('gray')
 for pp=1:length(vi)
  subplot(1,length(vi),pp)
  imagesc(vi{pp})
 end
end

nH = (ncon+nreg)*ngrp; % Columns of interest

%Ind = [ones(nscan,1) kron((1:ncon)',ones(totsub,1)) kron(ones(ncon,1),(1:totsub)') ones(nscan,1)];
SPM.xX = struct(...
        'X',X,...
        'iH',[1:nH],'iC',zeros(1,0),'iB',[(nH+1):size(X,2)],'iG',zeros(1,0),...
        'name',{cname},'I',Ind,...
        'sF',{{'repl'  'col'  'dummy'  'grp'}});
SPM.xC  = [];	
% SPM.xGX = struct(...
%         'iGXcalc',1,    'sGXcalc', sGXcalc,                               'rg',[],...
%         'iGMsca',9,     'sGMsca','<no grand Mean scaling>',...
%         'GM',0,         'gSF', SPM.xM,...
%         'iGC',  12,     'sGC',  '(redundant: not doing AnCova)',        'gc',[],...
%         'iGloNorm',9,   'sGloNorm','<no global normalisation>');

if nsph_flag
    SPM.xVi = struct('iid',0,'I',SPM.xX.I,'Vi',{vi} );       
else
    SPM.xVi = struct('iid',1,'V',speye(nscan) );       
end

% Mdes    = struct(...	
%         'Analysis_threshold',   {'None (-Inf)'},...
%         'Implicit_masking',     {'Yes: NaNs treated as missing'},...
%         'Explicit_masking',     {'No'});
% SPM.xM  = struct(...
%         'T',-Inf,'TH',ones(nscan,1)*-Inf,'I',1,'VM',mask,'xs',Mdes);
Pdes    = {{sprintf('%d condition, +%d covariate, +0 block, +0 nuisance',ncon*ngrp,nreg*ngrp); sprintf(', having %d degrees of freedom',rank(X)); sprintf('leaving %d degrees of freedom from %d images',nscan-rank(X),nscan)}};
SPM.xsDes = struct(...
        'Design',               {'Generic ANOVA with pooled error'},...
        'Global_calculation',   {sGXcalc},...
        'Grand_mean_scaling',   {'<no grand Mean scaling>'},...
        'Global_normalisation', {'<no global normalisation>'},...
        'Parameters',           Pdes);
    
save SPM SPM

%return

% Estimate parameters
%===========================================================================
SPM = spm_spm(SPM);


% Always Effects of interest contrast
%===========================================================================
try SPM = rmfield(SPM,'xCon'); end;

cn = 1;
c              = eye(nH);
if size(c,1)>1,  c=detrend(c,0); end
c              = [c zeros(size(c,1),size(SPM.xX.X,2)-nH)];
cname          = 'Unwhitened effects of interest';
SPM.xCon(cn)   = spm_FcUtil('Set',cname,'F','c',c',SPM.xX.xKXs);

% Additional user contrasts
%===========================================================================

if ~isempty(contrasts)
    for n=1:length(contrasts)
        cn = cn+1;
        c  = contrasts{n}.c;
        
        if size(c,2) ~= size(SPM.xX.X,2)
            if size(c,2) ~= nH
                error(sprintf('Contrast %d supplied does not have %d columns',n,nH))
            else
                c = [c zeros(size(c,1),size(SPM.xX.X,2)-nH)];
            end
        end
        
        if ~isfield(contrasts{n},'name')
            cname        = sprintf('User Con %d',n);
        else
            cname        = contrasts{n}.name;
        end
        
        SPM.xCon(cn) = spm_FcUtil('Set',cname,contrasts{n}.type,'c',c',SPM.xX.xKXs);
    end
end

spm_contrasts(SPM);

