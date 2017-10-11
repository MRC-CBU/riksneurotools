function [t,df,p,m]=t_matrix(ca,pf,tf,pv,dfc);

%% [t,df,p,m]=t_matrix(ca,pf,tf,pv,dfc);
%%
%% ca = matrix or cell array of data
%% pf = paired flag
%% tf = one-tailed flag
%% pv = pooled variance for unpaired t
%% dfc = Welch-Satterthwaite correction for unequal vars

% leading diagonal of T-matrix is a one-sample t against Ho=0

if(nargin<5)
    dfc=0;
end

if(nargin<4)
    pv=0;
end

if(pv==1 & dfc==1)
    error('Cannot assume pooled variance and correct df')
end

if(nargin<3)
    tf=0;
end

if(nargin<2)
    error('Need to provide 1) data and 2) binary flag as to whether repeated measures');
end


if(~iscell(ca))
    ca = mat2cell(ca,size(ca,1),ones(1,size(ca,2)));
end

n=size(ca);
n=max(n);

if(n==1 & pf~=1)
    'Error - need at least two groups'
    return;
end

t=zeros(n);

if pf		% paired t-tests
    
    for c1=1:n
        m(1,c1)=mean(ca{c1});
        
        for c2=c1:n
            
            if size(ca{c1},1)~=size(ca{c2},1)
                sprintf('Error - unequal sample sizes (conditions %d and %d) for paired t-test',c1,c2)
                return;
            end
            
            if min(size(ca{c1}))~=1
                sprintf('Error - more than one row/column in conditions %d and %d',c1,c2)
                return;
            end
            
            d1=ca{c1}(:);
            
            if(c1==c2)
                d2=0;
            else
                d2=ca{c2}(:);
            end
            
            d=d1-d2;
            t(c1,c2) = mean(d)*sqrt(length(d1))./std(d);
            df(c1,c2)= length(d1)-1;
            
            try
                p(c1,c2) = spm_Tcdf(t(c1,c2),df(c1,c2));
            catch
                p(c1,c2) = tcdf(t(c1,c2),df(c1,c2));;
            end
            
            if(p(c1,c2)>0.5)
                p(c1,c2)=1-p(c1,c2);
            end
            
        end
    end
    
    
else
    
    if n>2 % If more than 2 groups, then can only have one condition per group
        
        for c1=1:n
            m(1,c1)=mean(ca{c1},1);
                
            for c2=c1:n
                
                if min(size(ca{c1}))~=1 | min(size(ca{c2}))~=1
                   error('Error - more than two groups and more than one condition for groups %d and %d',c1,c2)
                end
                
                d1=ca{c1}(:);
                n1=length(d1);
                
                if c1==c2 & size(ca{c1},1) > 1
                    df(c1,c1)= n1-1;
                    t(c1,c1) = mean(d1)/(std(d1)/sqrt(length(d1)));
                else
                    
                    d2=ca{c2}(:);
                    n2=length(d2);
                    
                    df(c1,c2)=n1+n2-2;
                    
                    tpv = pv;
                    if (size(ca{c1},1) == 1 | size(ca{c2},1) == 1) & ~pv
                         warning('Group %d only has one subject, so pooling variance',c1)
                         tpv = 1;
                    end
                    
                    if tpv                     
                        sp = ( (n1-1)*var(d1) + (n2-1)*var(d2) ) / df(c1,c2);
                        
                        t(c1,c2) = (mean(d1) - mean(d2)) / sqrt( sp * (1/n1 + 1/n2) );
                    else
                        
                        t(c1,c2) = (mean(d1) - mean(d2)) / sqrt( var(d1)/n1 + var(d2)/n2 );
                        
                        if(dfc==1)
                            denom=(var(d1)/n1)^2/(n1-1) + (var(d2)/n2)^2/(n2-1);
                            df(c1,c2)=round((var(d1)/n1 + var(d2)/n2)^2/denom);
                        end
                    end
                end
                
                try
                    p(c1,c2) = spm_Tcdf(t(c1,c2),df(c1,c2));
                catch
                    p(c1,c2) = tcdf(t(c1,c2),df(c1,c2));
                end
                
                if(p(c1,c2)>0.5)
                    p(c1,c2)=1-p(c1,c2);
                end
            end
        end
        
    elseif n==2  % Only 2 groups, so compare each condition within groups
        
        if size(ca{1},2)~=size(ca{2},2)
            error('Error - unequal number of conditions for 2 groups')
        end
        n = size(ca{1},2); % overload n!
        t = zeros(n);
        m(1,:)=mean(ca{1},1); m(2,:)=mean(ca{2},1);
        n1=size(ca{1},1); n2=size(ca{2},1);
        
        tpv = pv;
        if (size(ca{1},1) == 1 | size(ca{2},1) == 1) & ~pv
            warning('One of two groups only has one subject, so pooling variance')
            tpv = 1;
        end
        
        for c1=1:n
            for c2=c1:n
                d1=ca{1}(:,c1); d2=ca{2}(:,c2);
                df(c1,c2)=n1+n2-2;
                
                if tpv
                    sp = ( (n1-1)*var(d1) + (n2-1)*var(d2) ) / df(c1,c2);
                    
                    t(c1,c2) = (mean(d1) - mean(d2)) / sqrt( sp * (1/n1 + 1/n2) );
                else
                    
                    t(c1,c2) = (mean(d1) - mean(d2)) / sqrt( var(d1)/n1 + var(d2)/n2 );
                    
                    if(dfc==1)
                        denom=(var(d1)/n1)^2/(n1-1) + (var(d2)/n2)^2/(n2-1);
                        df(c1,c2)=round((var(d1)/n1 + var(d2)/n2)^2/denom);
                    end
                end
                
                try
                    p(c1,c2) = spm_Tcdf(t(c1,c2),df(c1,c2));
                catch
                    p(c1,c2) = tcdf(t(c1,c2),df(c1,c2));
                end
                
                if(p(c1,c2)>0.5)
                    p(c1,c2)=1-p(c1,c2);
                end
            end
        end
    end
end

if ~tf 
    p=p.*2;
end

disp('Means...');  disp(m)
disp('T values...'), disp(t)
disp('dfs...'), disp(df)
disp(sprintf('P-values (%d tailed)...',2-tf)), disp(p)
disp('Bonferonni corrected "P-values" for all pairwise comparisons...'), disp(p*n*(n-1)/2)
