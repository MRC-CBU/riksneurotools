
function [Zmat, Bmat, pZmat, pBmat, aY, X0r] = rsfMRI_GLM(S);

% [Zmat, Bmat, pZmat, pBmat, aY, X0r] = rsfMRI_GLM(S);
%
% Function (using SPM8 functions) for estimating linear regressions 
% between fMRI timeseries in each pair of Nr ROIs, adjusting for bandpass 
% filter, confounding timeseries (eg CSF) and (SVD of) various expansions of 
% movement parameters, and properly modelling dfs based on comprehensive 
% model of error autocorrelation.
%
% rik.henson@mrc-cbu.cam.ac.uk, Jan 2013
%
% Many thanks to Linda Geerligs for pointing out improvements!
%
% S.Y = [Ns x Nr] data matrix, where Ns = number of scans (ie resting-state fMRI timeseries) and Nr = number of ROIs
% S.M = [Ns x 6] matrix of 6 movement parameters from realignment (x,y,z,pitch,roll,yaw)
% S.C = [Ns x Nc] matrix of confounding timeseries, Nc = number of confounds (eg extracted from WM, CSF or Global masks)
% S.G = [0/1] - whether to include global over ROIs (ie given current data) (default = 0)
% S.TR = time between volumes (TR), in seconds
%
% (S.HPC (default = 100) = highpass cut-off (in seconds)
% (S.LPC (default = 10) = lowpass cut-off (in seconds)
% (S.CY (default = {}) = precomputed covariance over pooled voxels (optional)
% (S.pflag (default=0) is whether to calculate partial regressions too (takes longer))
% (S.svd_thr (default=.99) is threshold for SVD of confounds)
% (S.SpikeMovAbsThr (default='', ie none) = absolute threshold for outliers based on RMS of Difference of Translations 
% (S.SpikeMovRelThr (default='', ie none) = relative threshold (in SDs) for outliers based on RMS of Difference of Translations or Rotations
% (S.SpikeDatRelThr (default='', ie none) = relative threshold (in SDs) for mean of Y over voxels (outliers, though arbitrary?)
% (S.SpikeLag (default=1) = how many TRs after a spike are modelled out as separate regressors
% (S.StandardiseY (default = 0) = whether to Z-score each ROI's timeseries to get standardised Betas)
% (S.PreWhiten (default = 1) = whether to estimate autocorrelation of error and prewhiten (slower, but better Z-values (less important for Betas?))
% (S.GlobMove = user-specified calculation of global movement)
%
% Zmat  = [Nr x Nr] matrix of Z-statistics for linear regression from seed (row) to target (column) ROI
% Bmat  = [Nr x Nr] matrix of betas for linear regression from seed (row) to target (column) ROI
% pZmat = [Nr x Nr] matrix of Z-statistics for partial linear regression from seed (row) to target (column) ROI
% pZmat = [Nr x Nr] matrix of betas for partial linear regression from seed (row) to target (column) ROI
% aY    = data adjusted for confounds
% X0r   = confounds (filter, confound ROIs, movement expansion)
%
% Note:
%    some Inf Z-values can be returned in Zmat (when p-value so low than inv_Ncdf is infinite)
%       (these could be replaced by the maximum Z-value in matrix?)
%
% Potential improvements:  
%   regularised regression (particularly for pZmat), eg L1 with LASSO - or L2 implementable with spm_reml_sc?

try Y = S.Y; catch error('Need Nscans x Nrois data matrix'); end
try C = S.C; catch error('Need Nscans x Nc matrix of Nc confound timeseries (Nc can be zero)'); end
try TR = S.TR; catch error('Need TR (in secondss)'); end

try HPC = S.HPC; catch 
    HPC = 1/0.01;
    warning('Assuming highpass cut-off of %d',HPC); 
end

try LPC = S.LPC; catch 
    LPC = 1/0.2;
    warning('Assuming lowpass cut-off of %d',LPC); 
end

try CY = S.CY; catch
    CY = [];
end

try GlobalFlag = S.G; catch
    GlobalFlag = 0;
end

try SpikeMovAbsThr = S.SpikeMovAbsThr; catch
    SpikeMovAbsThr = '';
%    SpikeMovAbsThr = 0.5;   % mm from Power et al
%    SpikeMovAbsThr = 0.25;  % mm from Satterthwaite et al 2013
end

try SpikeMovRelThr = S.SpikeMovRelThr; catch
%    SpikeMovRelThr = '';
    SpikeMovRelThr = 5;     % 5 SDs of mean?
end

try SpikeDatRelThr = S.SpikeDatRelThr; catch
    SpikeDatRelThr = '';
%    SpikeDatRelThr = 5;     % 5 SDs of mean?
end

try SpikeLag = S.SpikeLag; catch
    SpikeLag = 1;
%    SpikeLag = 5;     % 5 TRs after spike?
end

try StandardiseY = S.StandardiseY; catch
    StandardiseY = 0;
end

try PreWhiten = S.PreWhiten; catch
    PreWhiten = 1;
end

try VolterraLag = S.VolterraLag; catch
    VolterraLag = 5;  % artifacts can last up to 5 TRs = 10s, according to Power et al (2013)
end

try pflag   = S.pflag; catch     pflag = 0;   end
try svd_thr = S.svd_thr; catch svd_thr = .99; end

Ns = size(Y,1);
Nr = size(Y,2);


%% If want to try Matlab's LASSO (takes ages though)    
% lassoflag = 0;
% if lassoflag
%     matlabpool open
%     opts = statset('UseParallel','always');
% end


%% Create a DCT bandpass filter (so filtering part of model, countering Hallquist et al 2013 Neuroimage)
K   = spm_dctmtx(Ns,Ns);
nHP = fix(2*(Ns*TR)/HPC + 1);
nLP = fix(2*(Ns*TR)/LPC + 1);
K   = K(:,[2:nHP nLP:Ns]);      % Remove initial constant
Nk  = size(K,2);
fprintf('Bandpass filter using %d dfs (%d left)\n',Nk,Ns-Nk)


%% Create comprehensive model of residual autocorrelation (to counter Eklund et al, 2012, Neuroimage)
if PreWhiten
    T     = (0:(Ns - 1))*TR;                    % time
    d     = 2.^(floor(log2(TR/4)):log2(64));    % time constants (seconds)
    Q    = {};                                  % dictionary of components
    for i = 1:length(d)
        for j = 0:1
            Q{end + 1} = toeplitz((T.^j).*exp(-T/d(i)));
        end
    end
end


%% Detect outliers in movement and/or data
%M   = detrend(M,0);
%M(:,4:6)= M(:,4:6)*180/pi;

if ~isfield(S,'GlobMove')
    try M = S.M; catch error('Need Nscans x 6 movement parameter matrix'); end
    if size(M,2) == 6
        dM  = [zeros(1,6); diff(M,1,1)];    % First-order derivatives
        % Combine translations and rotations based on ArtRepair approximation for voxels 65mm from origin?
        %cdM = sqrt(sum(dM(:,1:3).^2,2) + 1.28*sum(dM(:,4:6).^2,2));
        dM(:,4:6) = dM(:,4:6)*50;  % Approximate way of converting rotations to translations, assuming sphere radius 50mm (and rotations in radians) from Linda Geerligs
        cdM       = sum(abs(dM),2);
    else
        error('If not user-specified global movement passed, then must pass 3 translations and 3 rotations')
    end
else
    cdM = S.GlobMove;
    try M = S.M; catch M=[]; end
 end

if ~isempty(SpikeMovAbsThr)  % Absolute movement threshold 
%    rms  = sqrt(mean(dM(:,1:3).^2,2));   % if want translations only    
%    aspk = find(rms > SpikeMovAbsThr);
    aspk = find(cdM > SpikeMovAbsThr);
    fprintf('%d spikes in absolute movement differences (based on threshold of %4.2f)\n',length(aspk),SpikeMovAbsThr)
else
    aspk = []; 
end

if ~isempty(SpikeMovRelThr)  % Relative (SD) threshold for translation and rotation
%     rms  = sqrt(mean(dM(:,1:3).^2,2));    
%     rspk = find(rms > (mean(rms) + SpikeMovRelThr*std(rms)));
%     rms  = sqrt(mean(dM(:,4:6).^2,2));    
%     rspk = [rspk; find(rms > (mean(rms) + SpikeMovRelThr*std(rms)))];
    rspk = find(cdM > (mean(cdM) + SpikeMovRelThr*std(cdM)));
    fprintf('%d spikes in relative movement differences (based on threshold of %4.2f SDs)\n',length(rspk),SpikeMovRelThr)
else
    rspk = [];
end

if ~isempty(SpikeDatRelThr)  % Relative (SD) threshold across all ROIs (dangerous? Arbitrary?)
    dY   = [zeros(1,Nr); diff(Y,1,1)]; 
    rms  = sqrt(mean(dY.^2,2));    
    dspk = find(rms > (mean(rms) + SpikeDatRelThr*std(rms)));
    fprintf('%d spikes in mean data across ROIs (based on threshold of %4.2f SDs)\n',length(dspk),SpikeDatRelThr)
else
    dspk = [];
end


%% Create delta-function regressors for each spike
spk = unique([aspk; rspk; dspk]); lspk = spk;
for q = 2:SpikeLag
    lspk = [lspk; spk + (q-1)];
end
spk = unique(lspk);

if ~isempty(spk)
    RSP = zeros(Ns,length(spk));
    n = 0;
    for p = 1:length(spk)
        if spk(p) <= Ns
            n=n+1;
            RSP(spk(p),n) = 1;
        end
    end 
    fprintf('%d unique spikes in total\n',length(spk))
    RSP = spm_en(RSP,0);
else
    RSP = [];
end


%% Create expansions of movement parameters
% Standard differential + second-order expansion (a la Satterthwaite et al, 2012)
% sM  = []; for m=1:6; for n=m:6; sM  = [sM M(:,m).*M(:,n)];    end; end   % Second-order expansion
% sdM = []; for m=1:6; for n=m:6; sdM = [sdM dM(:,m).*dM(:,n)]; end; end   % Second-order expansion of derivatives
% aM  = [M dM sM sdM];

% Above commented bits are subspace of more general Volterra expansion
%U=[]; for c=1:6; U(c).u = M(:,c); U(c).name{1}=sprintf('m%d',c); end; [aM,aMname] = spm_Volterra(U,[1 0 0; 1 -1 0; 0 1 -1]',2);
%U=[]; for c=1:6; U(c).u = M(:,c); U(c).name{1}=sprintf('m%d',c); end; [aM,aMname] = spm_Volterra(U,[1 0; 1 -1; 0 1]',2); %Only N and N+1 needed according to Satterthwaite et al 2013
bf = eye(VolterraLag);  % artifacts can last up to 5 TRs = 10s, according to Power et al (2013)
%    bf = [1 0 0 0 0; 1 1 0 0 0; 1 1 1 0 0; 1 1 1 1 0; 1 1 1 1 1];  % only leads to a few less modes below, and small compared to filter anyway!
bf = [bf; diff(bf)];
U=[]; for c=1:size(M,2); U(c).u = M(:,c); U(c).name{1}='c'; end; aM = spm_Volterra(U,bf',2);    
aM = spm_en(aM,0);

%% Add Global? (Note: may be passed by User in S.C anyway)
% (recommended by Rik and Power et al, 2013, though will entail negative
% correlations, which can be problem for some graph-theoretic measures)
if GlobalFlag
    C = [C mean(Y,2)];
end
C = spm_en(C,0);

%% Combine all confounds (assumes more data than confounds, ie Ns > Nc) and perform dimension reduction (cf PCA of correlation matrix)
%X0  = [C RSP K(:,2:end) aM];  % exclude constant term from K

XM  = spm_en(ones(Ns,1));  % constant term
X1  = [XM K RSP];          % Regressors that don't want to SVD
X0  = [aM C];              % Regressors that will SVD below
%X1  = [XM K C RSP];        % Regressors that don't want to SVD
%X0  = [aM];                % Regressors that will SVD below
X0r   = X0;
if ~isempty(X0)
    R1  = eye(Ns) - X1*pinv(X1);
    X0  = R1*X0;               % Project out of SVD-regressors what explained by non-SVD regressors!
    X0  = spm_en(X0,0);
    
    % Possible of course that some dimensions tiny part of SVD of X (so excluded) by happen to correlate highly with y...
    % ...could explore some L1 (eg LASSO) or L2 regularisation of over-parameterised model instead,
    % but LASSO takes ages (still working on possible L2 approach with spm_reml)
    if svd_thr < 1
        [U,S] = spm_svd(X0,0);
        S     = diag(S).^2; S = full(cumsum(S)/sum(S));
        Np    = find(S > svd_thr); Np = Np(1);
        X0r   = full(U(:,1:Np));
        fprintf('%d SVD modes left (from %d original terms) - %4.2f%% variance of correlation explained\n',Np,size(X0,2),100*S(Np))
    end
end

X0r = [K RSP X0r XM]; % Reinsert mean
Nc  = size(X0r,2);

if Nc >= Ns; error('Not enough dfs (scans) to estimate'); 
else fprintf('%d confounds for %d scans (%d left)\n',Nc,Ns,Ns-Nc); end


%% Create adjusted data, in case user wants for other metrics, eg, MI, and standardised data, if requested
R  = eye(Ns) - X0r*pinv(X0r);
aY = R*Y;

if StandardiseY
    Y = zscore(Y);
else
    Y = Y/std(Y(:));  % Some normalisation necessary to avoid numerical underflow, eg in cov(Y') below
end

%% Pool data covariance over ROIs (assuming enough of them!), unless specified
if isempty(CY) & PreWhiten
    CY = cov(Y');  % Correct to only do this after any standardisation?
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Main Loop
%% Do GLMs for all target ROIs for a given source ROI
Ar    = [1:Nr];
Zmat  = zeros(Nr);  % Matrix of Z statistics for each pairwise regression
Bmat  = zeros(Nr);

if pflag
    pZmat = zeros(Nr);  % Matrix of Z statistics for each pairwise partial regression
    pBmat = zeros(Nr);  % Matrix of Z statistics for each pairwise partial regression
else
    pZmat = [];
    pBmat = [];
end

%[V h]   = rik_reml(CY,X0r,Q,1,0,4);   % if could assume that error did not depend on seed timeseries
%W       = spm_inv(spm_sqrtm(V));
            
lNpc = 0;  % Just for printing pflag output below
for n = 1:Nr
       Ir = setdiff(Ar,n);
       
       Yr = Y(:,Ir);        % Yr contains all other timeseries, so ANCOVA below run on Nr-2 timeseries in one go
       Xr = Y(:,n);
%       Xr = spm_en(Xr);     % Normalise regressors so Betas can be compared directly across (seed) ROIs (not necessary if Standardised Y already)? 
       X  = [Xr X0r];
              
%        if lassoflag == 1  % takes too long (particularly for cross-validation to determine lambda)
%            for pn = 1:length(Ir)
%                Ir = setdiff(Ar,pn);
%                
%                Yr = Y(:,Ir);        % Yr contains all other timeseries, so ANCOVA below run on Nr-2 timeseries in one go
%                Xr = Y(:,pn);
%                Xr = spm_en(Xr,0);     % Normalise regressors so Betas can be compared directly across (seed) ROIs (not necessary is Standardised Y already)? 
%                
%                fprintf('l');
%                [lB,lfit] = lasso(X0r,Yr(:,pn),'CV',10,'Options',opts);
%                keepX     = find(lB(:,lfit.Index1SE));
%                X         = [Xr X0r(:,keepX)];
%                Nc        = length(keepX);
%                
%                [V h]   = rik_reml(CY,X,Q,1,0,4); % rik_reml is just version of spm_reml with fprintf commented out to speed up
%                W       = spm_inv(spm_sqrtm(V));
%                
%                %% Estimate T-value for regression of first column (seed timeseries)
%                [T(pn),dfall,Ball] = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Nc)]');
%                df(pn) = dfall(2);
%                B(pn)  = Ball(1);
%            end
%            fprintf('\n');
%        else

       %% Estimate autocorrelation of error (pooling across ROIs) and prewhitening matrix
       if PreWhiten
            [V h]   = rik_reml(CY,X,Q,1,0,4);   % rik_reml is just version of spm_reml with fprintf commented out to speed up
            W       = spm_inv(spm_sqrtm(V));
       else
            W       = speye(Ns);
       end
       
       %% Estimate T-value for regression of first column (seed timeseries)       
       [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr,[1 zeros(1,Nc)]');
       
       try 
           Zmat(n,Ir) = norminv(spm_Tcdf(T,df(2)));  % Needs Matlab Stats toolbox, but has larger range of Z (so not so many "Inf"s)
       catch
           Zmat(n,Ir) = spm_invNcdf(spm_Tcdf(T,df(2)));   
       end
       
       Bmat(n,Ir) = B(1,:);
       
       %% Estimate T-value for PARTIAL regression of first column (seed timeseries) - takes ages!
       if pflag               
           for pn = 1:length(Ir)

               pIr = setdiff(Ar,[n Ir(pn)]);
               pY  = spm_en(Y(:,pIr),0);       % Is necessary for SVD below 
               
               % SVD
               if svd_thr < 1
                   [U,S] = spm_svd([pY X0],0);
                   S     = diag(S).^2; S = full(cumsum(S)/sum(S));
                   Np    = find(S > svd_thr); Np = Np(1);
                   XY0   = full(U(:,1:Np));
                   Npc   = Np;
               else
                   XY0   = [pY X0];
               end
               
               XY0 = [XY0 ones(Ns,1)];  % Reinsert mean
               Npc = size(XY0,2);

%                if lassoflag == 1
%                    fprintf('l')
%                    [lB,lfit] = lasso(XY0,Yr(:,pn),'CV',10,'Options',opts);
%                    keepX     = find(lB(:,lfit.Index1SE));
%                    X         = [Xr XY0(:,keepX)];
%                    Nc        = length(keepX);
%                    
%                    [V h]   = rik_reml(CY,X,Q,1,0,4); % rik_reml is just version of spm_reml with fprintf commented out to speed up
%                    W       = spm_inv(spm_sqrtm(V));
%                    
%                    %% Estimate T-value for regression of first column (seed timeseries)
%                    [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Nc)]');
%                else
                    
               if Npc >= (Ns-1)  % -1 because going to add Xr below 
                   warning('Not enough dfs (scans) to estimate - just adjusting data and ignoring loss of dfs');
                   R = eye(Ns) - pY*pinv(pY);    % Residual-forming matrix
                   [T,df,B] = spm_ancova(W*X,speye(Ns,Ns),R*W*Yr(:,pn),[1 zeros(1,Npc)]');  % Ok to assume error and hence W unaffected by addition of pY in X? 
               else
                  
                   if lNpc ~= Npc  % Just to reduce time taken to print to screen
                        fprintf('   partial for seed region %d and target region %d: %d confounds for %d scans (%d left)\n',n,pn,Npc,Ns,Ns-Npc);  
                   end
                   lNpc = Npc;
                   
                   X  = [Xr XY0];    
                   
                   %% Commented out below because ok to assume error and hence W unaffected by addition of pY in X? 
%                   [V h]   = rik_reml(CY,X,Q,1,0,4);   % rik_reml is just version of spm_reml with fprintf commented out to speed up
%                   W       = spm_inv(spm_sqrtm(V));
                   
                   %% Estimate T-value for regression of first column (seed timeseries)
                   [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Npc)]');
               end

               try 
                   pZmat(n,Ir(pn)) = norminv(spm_Tcdf(T,df(2)));
               catch
                   pZmat(n,Ir(pn)) = spm_invNcdf(spm_Tcdf(T,df(2)));
               end
               
               pBmat(n,Ir(pn)) = B(1);
           end
       end
       fprintf('.')
end
fprintf('\n')

return

%figure,imagesc(Zmat); colorbar
%figure,imagesc(pZmat); colorbar



function [V,h,Ph,F,Fa,Fc] = rik_reml(YY,X,Q,N,D,t,hE,hP)
% ReML estimation of [improper] covariance components from y*y'
% FORMAT [C,h,Ph,F,Fa,Fc] = rik_reml(YY,X,Q,N,D,t,hE,hP);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
%
% N   - number of samples                 (default 1)
% D   - Flag for positive-definite scheme (default 0)
% t   - regularisation                    (default 4)
% hE  - hyperprior                        (default 0)
% hP  - hyperprecision                    (default 1e-16)
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of h
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates.
%
% see also: spm_reml_sc for the equivalent scheme using log-normal
% hyperpriors
%__________________________________________________________________________
%
% SPM ReML routines:
%
%      spm_reml:    no positivity constraints on covariance parameters
%      spm_reml_sc: positivity constraints on covariance parameters
%      spm_sp_reml: for sparse patterns (c.f., ARD)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Karl Friston
% $Id: spm_reml.m 5223 2013-02-01 11:56:05Z ged $
% Modified by Rik to remove screen output and turn off warnings
 
% check defaults
%--------------------------------------------------------------------------
try, N;  catch, N  = 1;     end       % assume a single sample if not specified
try, K;  catch, K  = 32;    end       % default number of iterations
try, D;  catch, D  = 0;     end       % default checking
try, t;  catch, t  = 4;     end       % default regularisation
try, hE; catch, hE = 0;     end       % default hyperprior
try, hP; catch, hP = 1e-16; end       % default hyperprecision
 
% catch NaNs
%--------------------------------------------------------------------------
W     = Q;
q     = find(all(isfinite(YY)));
YY    = YY(q,q);
for i = 1:length(Q)
    Q{i} = Q{i}(q,q);
end
 
% dimensions
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
 
% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
else
    X = spm_svd(X(q,:),0);
end
 
% initialise h and specify hyperpriors
%==========================================================================
h   = zeros(m,1);
for i = 1:m
    h(i,1) = any(diag(Q{i}));
end
hE  = sparse(m,1) + hE;
hP  = speye(m,m)*hP;
dF  = Inf;
D   = 8*(D > 0);
 
warning off

% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:K

    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        C = C + Q{i}*h(i);
    end
 
    % positive [semi]-definite check
    %----------------------------------------------------------------------
    for i = 1:D
        if min(real(eig(full(C)))) < 0

            % increase regularisation and re-evaluate C
            %--------------------------------------------------------------
            t     = t - 1;
            h     = h - dh;
            dh    = spm_dx(dFdhh,dFdh,{t});
            h     = h + dh;
            C     = sparse(n,n);
            for i = 1:m
                C = C + Q{i}*h(i);
            end
        else
            break
        end
    end


    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iC     = spm_inv(C);
    iCX    = iC*X;
    if ~isempty(X)
        Cq = spm_inv(X'*iCX);
    else
        Cq = sparse(0);
    end

    % M-step: ReML estimate of hyperparameters
    %======================================================================

    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    P     = iC - iCX*Cq*iCX';
    U     = speye(n) - P*YY/N;
    for i = 1:m

        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        PQ{i}     = P*Q{i};
        dFdh(i,1) = -spm_trace(PQ{i},U)*N/2;

    end

    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m

            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j) = -spm_trace(PQ{i},PQ{j})*N/2;
            dFdhh(j,i) =  dFdhh(i,j);

        end
    end
 
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;

    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = spm_dx(dFdhh,dFdh,{t});
    h     = h + dh;

    % predicted change in F - increase regularisation if increasing
    %----------------------------------------------------------------------
    pF    = dFdh'*dh;
    if pF > dF
        t = t - 1;
    else
        t = t + 1/4;
    end
    
    % revert to SPD checking, if near phase-transition
    %----------------------------------------------------------------------
    if ~isfinite(pF) || abs(pF) > 1e6
        [V,h,Ph,F,Fa,Fc] = rik_reml(YY,X,Q,N,1,t - 2);
        return
    else
        dF = pF;
    end
    
    % Convergence (1% change in log-evidence)
    %======================================================================
%    fprintf('%s %-23d: %10s%e [%+3.2f]\n','  ReML Iteration',k,'...',full(pF),t);
 
    % final estimate of covariance (with missing data points)
    %----------------------------------------------------------------------
    if dF < 1e-1, break, end
end

 
% re-build predicted covariance
%==========================================================================
V     = 0;
for i = 1:m
    V = V + W{i}*h(i);
end
 
% check V is positive semi-definite (if not already checked)
%==========================================================================
if ~D
    if min(eig(V)) < 0
        [V,h,Ph,F,Fa,Fc] = rik_reml(YY,X,Q,N,1,2,hE(1),hP(1));
        return
    end
end
 
% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
Ph    = -dFdhh;
if nargout > 3
 
    % tr(hP*inv(Ph)) - nh + tr(pP*inv(Pp)) - np (pP = 0)
    %----------------------------------------------------------------------
    Ft = trace(hP*inv(Ph)) - length(Ph) - length(Cq);
 
    % complexity - KL(Ph,hP)
    %----------------------------------------------------------------------
    Fc = Ft/2 + e'*hP*e/2 + spm_logdet(Ph*inv(hP))/2 - N*spm_logdet(Cq)/2;
 
    % Accuracy - ln p(Y|h)
    %----------------------------------------------------------------------
    Fa = Ft/2 - trace(C*P*YY*P)/2 - N*n*log(2*pi)/2 - N*spm_logdet(C)/2;
 
    % Free-energy
    %----------------------------------------------------------------------
    F  = Fa - Fc;
 
end

warning on


function [C] = spm_trace(A,B)
% fast trace for large matrices: C = spm_trace(A,B) = trace(A*B)
% FORMAT [C] = spm_trace(A,B)
%
% C = spm_trace(A,B) = trace(A*B) = sum(sum(A'.*B));
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_trace.m 4805 2012-07-26 13:16:18Z karl $

% fast trace for large matrices: C = spm_trace(A,B) = trace(A*B)
%--------------------------------------------------------------------------
C = sum(sum(A'.*B));


