function [e,sots,stim,X,df] = fMRI_GLM_efficiency(S);

% Matlab function to calculate efficiency for fMRI GLMs.
%                     rik.henson@mrc-cbu.cam.ac.uk, 2005
%
% See for example:
%
% Henson, R.N. (2006). Efficient experimental design for fMRI. In K. Friston, J. Ashburner, S. Kiebel, T. Nichols, and W. Penny (Eds), Statistical Parametric Mapping: The analysis of functional brain images. Elsevier, London, 2006. pp. 193-210. 
% http://www.mrc-cbu.cam.ac.uk/personal/rik.henson/personal/Henson_SPM_06_preprint.pdf
%
% http://imaging.mrc-cbu.cam.ac.uk/imaging/DesignEfficiency
%
% ...and used for figures in Henson (in press) "Design Efficiency". Brain Mapping: A
% comprehensive reference.
%
% Needs SPM5+ on path
% Written for ease of exposition, not speed!
%
% Function take one of:
%
%  1. cell array of stimulus onset times (SOTS), and possibly durations, eg from SPM
%  2. vector of stimulus codes (ie order of events) (if onset times not calculated)
%  3. a transition matrix in order to generate 1+2 (if no onsets or stimulus train available)
%
% Inputs (fields of structure S):
%
%    S.Ns = Number of scans (REQUIRED)
%    S.CM = Cell array of T or F contrast matrices (NOTE THAT THE L1 NORM OF CONTRASTS SHOULD 
%           MATCH IN ORDER TO COMPARE THEIR EFFICIENCIES, EG CM{1}=[1 -1 1 -1]/2 CAN BE 
%           COMPARED WITH CM{2}=[1 -1 0 0] BECAUSE NORM(Cm{i},1)==2 IN BOTH CASES
%    S.sots    = cell array of onset times for each condition in units of scans (OPTIONAL; SEE ABOVE)
%    S.durs    = cell array of durations for each condition in units of scans (OPTIONAL; SEE ABOVE)
%    S.stim    = vector of events (OPTIONAL; SEE ABOVE)
%    S.TM.prev = Np x Nh transition matrix history of Np possible sequences of the previous Nh events (OPTIONAL; SEE ABOVE)
%    S.TM.next = Np x Nj transition matrix of probabilities of each of next Nj event-types for each of Np possible histories (OPTIONAL; SEE ABOVE)
%    S.SOAmin  = minimal SOA (REQUIRED IF STIMULUS TRAIN OR TRANSITION MATRIX PROVIDED)
%    S.Ni      = Number of stimuli (events) (REQUIRED IF ONLY TRANSITION MATRIX PROVIDED)
%    S.bf = type of HRF basis function (based on spm_get_bf.m) (DEFAULTS TO CANONICAL HRF)
%    S.HC = highpass filter cut-off (s) (DEFAULTS TO 120)
%    S.TR = inter-scan interval (s) (DEFAULTS TO 2)
%    S.t0 = initial transient (s) to ignore (DEFAULTS TO 30)
%
% Outputs:
%
%    e    = efficiency for each contrast in S.CM
%    sots = stimulus onset times (in scans), eg for input into SPM
%    stim = vector of stimulus codes (event-types), just for interest!
%    X    = GLM design matrix
%    df   = degrees of freedom left given model and filter
%
% Examples:
%
%   0. Simple user-specified blocked design
%
% S=[];
% S.TR = 1; 
% S.Ns = 1000;
% S.sots{1} = [0:48:1000];   % Condition 1 onsets
% S.sots{2} = [16:48:1000];  % Condition 2 onsets
% S.durs{1} = 16; S.durs{2} = 16;  % Both 16s epochs
% S.CM{1} = [1 -1]; S.CM{2} = [1 -1]; % Contrast conditions, and conditions vs baseline
% 
% [e,sots,stim,X,df] = fMRI_GLM_efficiency(S);
% 
%
%   1. Code-generated, randomised design with 2 event-types, canonical HRF interest in both differential and common effects
%
% S=[];
% S.Ni = 2000;
% S.TR = 2;
% S.Ns = 1000;          % Must be more than S.Ni*S.SOAmin/TR
% S.CM{1} = [1 -1];     % A-B
% S.CM{2} = [1  1];     % A+B
% S.HC = 120;
% S.TM.prev = [1 2]';
% S.TM.next = [0.5 0.5; 0.5 0.5];
% S.bf = 'hrf';
% 
% soas = [1:20]  % to explore efficiency as function of SOA
%
% S.SOAmin = soas(1); 
% [e(1,:),sots,stim,X,df] = fMRI_GLM_efficiency(S);  % call once to generate stimulus train
% S.stim = stim % ensure same event-train for all SOA
%
% for s = 2:length(soas)
%    S.SOAmin = soas(s);
%    e(s,:)   =  fMRI_GLM_efficiency(S);
% end
% figure,plot(soas,e)
%
%
%    2. Randomised design with FIR and null events (S continued from above)
%
% S.bf = 'Finite Impulse Response'
% S.TM.prev = [1 2]';
% S.TM.next = [1/3 1/3; 1/3 1/3];
% S = rmfield(S,'stim');
% S.SOAmin = 2; e=[];
% [e,sots,stim,X,df] = fMRI_GLM_efficiency(S);
%
%
%    3. Blocked design with canonical HRF
%
% S.bf = 'hrf'
% S = rmfield(S,'stim');
% 
% bls = [1:60]; e=[];
% for s = 1:length(bls)
%     bl = bls(s); 
%     try  S = rmfield(S,'TM'); end
%     for b=1:bl
%         S.TM.prev(b,:)    = [ones(1,bl-(b-1)) 2*ones(1,b-1)];
%         S.TM.prev(b+bl,:) = [2*ones(1,bl-(b-1)) ones(1,b-1)];
%         S.TM.next(b,:)    = [0 1];
%         S.TM.next(b+bl,:) = [1 0];
%     end
%     e(s,:) =  fMRI_GLM_efficiency(S);
% end
% figure,plot(bls,e)

try
    Ns = S.Ns;
catch
    error('Please provide number of scans (volumes) in S.Ns')
end
df = Ns;

try
    CM = S.CM;
catch
    error('Must pass a Nc x Nj matrix of contrasts, where Nc = number of contrasts and Nk = number of evemt-types')
end

try
    TR = S.TR;
catch
    TR = 2;  % Assume 2s TR unless specified otherwise
end

try
    t0 = S.t0;
catch
    t0 = 30;  % Ignore initial 30s transient
end

dt = 0.2;   % Minimum time for convolution space (seconds)
st = round(TR/dt);
t0 = round(t0/TR);

rand('state',0);
try
    sots = S.sots;
    stim = [];  % Can't be bothered to create stim if sots provided!
catch
    try
        SOAmin = S.SOAmin;
        ds     = SOAmin/TR; 
        stim   = S.stim;
        sots   = gensots(stim,ds);
    catch
        try
            TM.prev = S.TM.prev;
            TM.next = S.TM.next;
            Ni      = S.Ni;
            stim    = genstim(TM,Ni);
            sots    = gensots(stim,ds);
        catch
            error('Either provide onsets in S.sots, stimuli in S.stim and SOA in S.SOAmin, or transition matrix in S.TM and number of stimuli in Ni')
        end
    end
end

try
    durs = S.durs;
    for j = 1:length(sots)
        if length(durs{j}) == 1
            durs{j} = repmat(durs{j},1,length(sots{j}));
        elseif length(durs{j}) ~= length(sots{j})
            error('Number of durations (%d) does not match number of onsets (%d) for condition %d',length(durs{j}),length(ons{j}),j)
        end
    end
catch
    for j = 1:length(sots)
        durs{j} = zeros(1,length(sots{j}));  % assume events if durations not specified       
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
else
    % Assume bf is a TxNk matrix of basis functions in units of dt
end

try
    HC = S.HC;
catch
    HC = 120;
end

if HC > 0
    HO = fix(2*(Ns*TR)/HC+1);
    K  = spm_dctmtx(Ns,HO);
    df = df - size(K,2);
    K  = eye(Ns) - K*K';
else
    K  = eye(Ns);
end

Nj = length(sots);
Nk = size(bf,2);
Nt = Ns*st;

%% Create neural timeseries
s = [1:st:Nt];
X = zeros(Ns,Nj*Nk);
for j = 1:Nj

    u = zeros(Nt,1);
    for i = 1:length(sots{j})
        o = [sots{j}(i)*st : st : (sots{j}(i)+durs{j}(i))*st];
        u(round(o)+1) = 1;
    end
    
    for k = 1:Nk
        b = conv(u,bf(:,k));
        X(:,(j-1)*Nk + k) = b(s);
    end
end

%% highpass filter and remove any transients and mean
X = K*X;  
X = X(t0:end,:);  
X = detrend(X,0);
df = df - size(X,2) - 1;

%cc = corrcoef(X*CM');
%max(abs(cc-eye(size(cc,2))))

%% Calculate efficiency for each contrast (row of CM) and average for all contrasts

Nc = length(CM);
bCM = {};
for c = 1:Nc
    bCM{c} = kron(CM{c},eye(Nk)/Nk);
end

iXX = pinv(X'*X);
for c = 1:Nc
	e(1,c) = 1/trace(bCM{c}*iXX*bCM{c}');
    fprintf('Contrast %d: Efficiency = %f\n',c,e(1,c))
end
fprintf('\n')

return





function sots = gensots(stim,ds);

Ni   = length(stim);
Aj   = unique(stim(find(stim>0)));

if Aj(1)~=1 | any(diff(Aj)~=1)
    error('Stimuli must be numbered ordinally from 1 onwards')
end

Nj = length(Aj);
sots = cell(Nj,1);

for i = 1:Ni    
    if stim(i) > 0
        sots{stim(i)}(end+1) = i*ds;
    end
end

return


function stim = genstim(TM,Ni);

Nh = size(TM.prev,2);
Np = size(TM.prev,1);
Nj = size(TM.next,2);

if any(TM.next(:) > 1) | any(TM.next(:) < 0) | any(sum(TM.next,2) > 1)
    error('All values in TM.next must be [0,1], and sum along rows must be [0,1]')
end

hist = ones(1,Nh)*NaN;
stim = [];

for i = 1:Ni     
    h = rfind(hist,TM.prev,Nj);
    r = rand;
    
    if r <= sum(TM.next(h,:))
        j = 1;
        while j < Nj & r > TM.next(h,j)
            r = r-TM.next(h,j);
            j = j+1;
        end
        
        hist = [hist(2:end) j];
        stim(end+1) = j;
    else
        stim(end+1) = 0;   % Null event
    end
end

return


function r = rfind(row,matrix,Nj);

Nr = size(matrix,1);
Nc = size(matrix,2);

h = find(isfinite(row));

%Nj = Nr^(1/Nc);

if isempty(h)
    r = floor(rand*Nj)+1;
else
    for r = 1:Nr
        if all((row(h) - matrix(r,h)) == 0)
            break;
        elseif r==Nr
            error('Should not get here if TM correct!');
        end
    end
end

return




