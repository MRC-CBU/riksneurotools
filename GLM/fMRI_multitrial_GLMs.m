function Beta = fMRI_multitrial_GLMs(S);

% Matlab function to return single-trial Betas from GLM using LSA, LSS or L2-regularised LSA
% fit to fMRI timeseries from a scan x ROI data matrix, given onsets and durations of events
%
% rik.henson@mrc-cbu.cam.ac.uk, 2017 (see Abdulrahman & Henson, 2016, Neuroimage)
%
% Inputs (fields of structure S):
%
%    Required:
%           S.d = Ns scan x Nr ROI matrix of fMRI data
%           S.events = event onsets and durations for each condition in SPM format
%
%    Optional:
%           S.method = GLM method, eg LSU, LSA or LSS (see Abdulrahman & Henson, 2016, Neuroimage)
%                      default is LSA
%           S.lambda = regularisation parameter (default = 0 = no regularisation)
%           S.XC = any confounding regressors (eg movement parameters), as
%                  matrix of Ns by number of confounds (default = none)
%           S.coi = "conditions of interest" - indices of conditions for
%                   which single-trial estimates are desired (default=all)
%           S.TR = TR (default=2s)
%           S.units = TR/sec (default=sec)
%           S.T = microtime resolution (default=16)
%           S.T0 = reference timebin (default=8)
%           S.bf = HRF basis functions (default = SPM's canonical HRF)
%           S.HC = highpass cutoff for filtering (default=128s)
%           S.zflag = whether to Z-score regressors (default=0=no) 
%
% Outputs:
%
%    Betas - note - after Z-scoring all columns
%
% Future additions
%   Could add AR(1) modelling


try
    d = S.d;
catch
    error('Please provide ROI fMRI data in S.d')
end

Ns = size(d,1);
Nr = size(d,2);
if Ns <= 1, error('Need at least 2 scans (rows of S.d)'); end

try
    TR = S.TR;
catch
    TR = 2;  % Assume 2s TR unless specified otherwise
end

try
    events = S.events;
catch
    error('Please provide SPMs event-structure in S.events');
end   

try
    units = S.units;
catch
    warning('Assuming units are seconds');
    units = 'secs';
end   

try
    T = S.T;
catch
    T = 16;
end   

dt = TR/T;   % Minimum time for convolution space (seconds)

if strcmp(units,'secs')
    st = T/TR; % If onsets in seconds
elseif strcmp(units,'scans')
    st = T;  % If onsets in scans
else
    error('Units must be secs or scans');
end
Nt = Ns*T;

try
    T0 = S.T0;
catch
    warning('Synchronising events with middle of scan');
    T0 = round(T/2);
end   

if T0<1 | T0>T
    error('T0 %d must lie between 0 and T %d',T0,T);
end


sots = events.ons;
Nj = length(sots);
if Nj < 1, error('Need at least 1 trial-type'); end

try
    coi = S.coi;  % Subset of condition of interest
catch
    coi = 1:Nj;
end

%nams = S.events.nam;

try
    durs = S.events.dur;
catch
    durs = cell(1,Nj)
    for j=1:Nj
        durs{j} = zeros(1,length(sots{j})); % assume all events
    end
end

try
    bf = S.bf;
catch
    bf = 'hrf';  % SPM's canonical HRF
%    bf = 'Finite Impulse Response'  % FIR
end

if isstr(bf)
    xBF.dt     = dt;
    xBF.name   = bf;
    xBF.length = 30;
    xBF.order  = 30/TR;
    bf = spm_get_bf(xBF);
    bf = bf.bf;
    bf = bf/max(bf(:));
else
    % Assume bf is a TxNk matrix of basis functions in units of dt
end
Nk = size(bf,2);

try
    HC = S.HC;
catch
    HC = 128;
end

if HC > 0
    HO = fix(2*(Ns*TR)/HC+1);
    K  = spm_dctmtx(Ns,HO);
%    df = df - size(K,2);
%
else
     K = [];
%    K  = eye(Ns);
end

try 
    XC = S.XC;
catch
    XC = [];
end

try 
    meth = S.method;
catch
    meth = 'LSA';
end

try 
    lambda = S.lambda;
catch
    lambda = 0; % Default to no regularisation of LSU/LSA
end

try 
    zflag = S.zflag;
catch
    zflag = 0;  % Default to Z-scoring design matrix columns
end


%% Create GLMs

s = [T0:T:Nt];

X = []; X0 = []; Beta = {};
 
switch meth
    case 'LSU'      
        for j = 1:Nj
            u = zeros(Nt,1);
            Ni = length(sots{j});
            for i = 1:Ni    
                t1 = round(sots{j}(i)*st)+1;
                t2 = t1+round(durs{j}(i)*st);
                u(t1:t2) = 1;
            end
            for k = 1:Nk
                b = conv(u,bf(:,k));
                if ismember(j,coi)
                    X(:,end+1) = b(s); 
                else
                    X0(:,end+1) = b(s);
                end               
            end
        end            

        Z = [X X0 XC K];
        if zflag         
            Z = zscore(Z);   % zscore works on constant column
        end
%        figure,imagesc(Z);
        
%        Adjusting data and model by X0 before regularisation gives same answer as below
%        
%         X0 = Z(:,(Nji+1):end); X = Z(:,1:Nji);
%         Z = orthog(X,X0);
%         ad = orthog(d,X0);   
%         pX = inv(Z'*Z + lambda*eye(Nji))*Z';
%         B = pX*ad;
        
        Nji = length(coi);
        if lambda == 0
            pX = pinv(Z);  % Needed for LSA below, so use here too
        else    
            Np = size(Z,2);
            R = [[eye(Nji) zeros(Nji,Np-Nji)]; zeros(Np-Nji,Np)];       
            pX = inv(Z'*Z + lambda*R)*Z';   
        end
        
        B = pX*d;
        
        for j = 1:Nji
            Beta{j} = B(j,:);
        end
    
    case 'LSA'
        ri = {}; lastri = 0;
        for j = 1:Nj
            if ismember(j,coi)
                Ni = length(sots{j});
                u  = zeros(Nt,Ni);
                for i = 1:Ni
                    t1 = round(sots{j}(i)*st)+1;
                    t2 = t1+round(durs{j}(i)*st);
                    u(t1:t2,i) = 1;
                end
                ri{find(coi==j)} = [1:(Ni*Nk)] + lastri;
                lastri = ri{find(coi==j)}(end);
                for i = 1:Ni
                    for k = 1:Nk
                        b = conv(u(:,i),bf(:,k));
                        X(:,end+1) = b(s);
                    end
                end              
            else
                Ni = 1;
                u = zeros(Nt,1);
                u(round(sots{j}*st)+1) = 1;
                for i = 1:Ni
                    for k = 1:Nk
                        b = conv(u(:,i),bf(:,k));
                        X0(:,end+1) = b(s);
                    end
                end
            end
        end
      
        Z = [X X0 XC K];
        if zflag         
            Z = zscore(Z);   % zscore works on constant column
        end
%        figure,imagesc(Z);
        
        if lambda == 0
            pX = pinv(Z);  % Needed for LSA below, so use here too
        else    
            Nji = length(coi); Np = size(Z,2);
            R = [[eye(Nji) zeros(Nji,Np-Nji)]; zeros(Np-Nji,Np)];       
            pX = inv(Z'*Z + lambda*R)*Z';   
        end
        
        B = pX*d;
        
        for j = 1:length(coi)
            Beta{j} = B(ri{j},:);
        end
        
    case 'LSS'
        ri = {}; lastri = 0;
        for j = 1:Nj
            if ismember(j,coi)
                Ni = length(sots{j});
                u  = zeros(Nt,Ni);
                for i = 1:Ni
                    t1 = round(sots{j}(i)*st)+1;
                    t2 = t1+round(durs{j}(i)*st);
                    u(t1:t2,i) = 1;
                end
                ri{find(coi==j)} = [1:(Ni*Nk)] + lastri;
                lastri = ri{find(coi==j)}(end);
                for i = 1:Ni
                    for k = 1:Nk
                        b = conv(u(:,i),bf(:,k));
                        X(:,end+1) = b(s);
                    end
                end              
            else
                Ni = 1;
                u = zeros(Nt,1);
                u(round(sots{j}*st)+1) = 1;
                for i = 1:Ni
                    for k = 1:Nk
                        b = conv(u(:,i),bf(:,k));
                        X0(:,end+1) = b(s);
                    end
                end
            end
        end
        
        Np = size(X,2);
        for j = 1:length(coi)
            Ni = length(sots{coi(j)});  
            Xl = [];
            for k = setdiff(1:length(coi),j) 
                Xl      = [Xl sum(X(:,ri{k}),2)];
            end
            for i = 1:Ni
                Xs = [];
                Xs(:,1) = X(:,ri{j}(i));
                if Ni>1
                    Xs(:,2) = sum(X(:,setdiff(ri{j},ri{j}(i))),2);
                end
                
                Z = [Xs Xl X0 XC K];
                if zflag
                    Z = zscore(Z);   % zscore works on constant column
                end
                %        figure,imagesc(Z);
                
                pX = pinv(Z);  % Needed for LSA below, so use here too
               
                B = pX*d;
        
                Beta{j}(i,:) = B(1,:);
            end
        end
        
    otherwise   
        error('unknown estimation method')
end


return


