function [y,P,epsilon] = rp_TRP(varargin)
%RP_TRP computes a RP without vertical lines after Ahlstrom et al.2006
%
% Minimum input-arguments : 1
% Maximum input-arguments : 6
%
%    [R,DM,epsilon] = rp_TRP(Y,E,thres_meth,norm,rt_thres,algorithm) 
%    
%    Calculates a recurrence plot 'R' from an phase space vector 'Y', using 
%    the threshold 'E'. 'DM' is an optional output and is the adjacency- or
%    distancematrix. In case you choose the 'var'-fixed threshold selection
%    method optional output 'epsilon' will give the actual calculated
%    value.
%
%
% Input:
%
% 'Y'                       is a N-by-M matrix corresponding to N time points
%                           and M embedding/phase space dimensions.
%
% 'norm' (optional)         norm for distance calculation in phasespace to
%                           'euc' (euclidic) or 'max' (maximum). Default is 
%                           max norm.
%
% 'rt_thres' (optional)     sets the minimum distance (in sampling units),
%                           referred to as recurrence times of first type,
%                           between two consecutive recurrences in a column
%                           of the RP. Default is rt_thres = 1.
%
% 'algorithm' (optional)    specify the way of calculating the distance
%                           matrix here. You can choose from
%                           ['loops','vector','matlabvector']. Default is
%                           'vector'.
%
% 'thres_meth' (optional) specifies how the recurrence threshold epsilon 
% will be calculated. There are three options. Set 'thres_meth' to
%   - 'fix' The RP is computed under a fixed threshold epsilon specified by
%           input parameter 'E'.
%   - 'var' The RP is computed under a fixed threshold epsilon, which
%           corresponds to the lower 'E'-quantile (specified by input parameter
%           'E') of the distance distribution of all points in phasespace.
%   - 'fan' The RP is computed under a variable threshold epsilon using a
%           fixed amount of nearest neighbours in phasespace to compute the
%           epsilon-value for each point of the phasespace trajectory
%           individually.
% Default is 'fix'.  
%
%
%    Example (CRP toolbox needs to be installed):
%      x = sin(linspace(0,5*2*pi,1000));
%      xe = embed(x,2,50);
%      r = rp(xe,.2);
%      [r2,~,~] = rp_TRP(xe,.2,'var','euc',50);
%      figure
%      subplot(1,2,1)
%      imagesc(r), colormap([1 1 1; 0 0 0]), axis xy square
%      title('input RP')
%      subplot(1,2,2)
%      imagesc(r2), colormap([1 1 1; 0 0 0]), axis xy square
%      title('TRP')
%
%
%
% Copyright (c) 2019-
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%% check input
narginchk(1,6)
nargoutchk(0,3)

x = varargin{1};

N = size(x);
if N(1) < N(2)
   warning('Embedding dimension is larger than the length of the vector. Please check!')
end

algoLib={'loops','vector','matlabvector'}; % the possible algorithms
try
    algorithm = varargin{6};
    if ~isa(algorithm,'char') || ~ismember(algorithm,algoLib)
        warning(['Specified algorithm should be one of the following possible values:',...
           10,sprintf('''%s'' ',algoLib{:})])
    end
catch
    algorithm = 'vector';
end


try
    rt_thres = varargin{5};
    if rem(rt_thres,1)~=0 || rt_thres < 1 || rt_thres >= N(1)
        warning('recurrence time threshold needs to be a positive integer with maximum size of the RP. Now set to rt_thres = 1')
        rt_thres = 1;
    end
catch
    rt_thres = 1;
end

methLib={'euc','max'}; % the possible norms
try
    meth = varargin{4};
    if ~isa(meth,'char') || ~ismember(meth,methLib)
       warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
    end
catch
    meth = 'max';
end

thresLib={'fix','var','fan'}; % the possible ways of threshold computation
try
    thres = varargin{3};
    if ~isa(thres,'char') || ~ismember(thres,thresLib)
       warning(['Specified way of calculating threshold should be one of the following possible values:',...
                                10,sprintf('''%s'' ',thresLib{:})])
    end
catch
    thres = 'fix';
end

try
    e = varargin{2};
catch
    e = 1;
end




%% init output

%% bind length of input vector in order to have constant iteration bounds while using parallel computing
M=N(1);
%% calculate distance matrix
switch algorithm 
    case 'loops'
         %% calculation with loops
         y = zeros(N(1),N(1));
         parfor i = 1:M
               for j = 1:M
                switch lower(meth)
                    case 'euc'
                         d = (x(i,:) - x(j,:)).^2;
                         y(i,j) = sqrt(sum(d));
                    otherwise
                         d = abs(x(i,:) - x(j,:));
                         y(i,j) = max(d);
                end
               end
         end
   case 'vector'
    
        %% calculation with vectorisation
        x1 = repmat(x,N(1),1);
        x2 = reshape(repmat(reshape(x,N(1)*N(2),1),1,N(1))',N(1)*N(1),N(2));
        switch lower(meth)
          case 'euc'
              d = (x1 - x2).^2;
              y = sqrt(sum(d,2));
          otherwise
              d = abs(x1 - x2);
              y = max(d,[],2);
        end
        y = reshape(y,N(1), N(1));   
        
    case 'matlabvector'
      %% calculation with matlab's vectorisation
      switch lower(meth)
          case 'euc'
              y = squareform(pdist(x));
          otherwise
              y = squareform(pdist(x,'chebychev'));
      end
end
P = y;

if strcmp(thres,'fix')
    % apply threshold
    y = double(y < e);
    epsilon = e;
    
elseif strcmp(thres,'var')    
    % get lower (e*100)%-quantile of distance-distribution
    epsilon = quantile(P(:),e);
    y = double(y < epsilon);        

elseif strcmp(thres,'fan')
    % compute variable threshold for each point in order to get fixed
    % number of nearest neighbours
    q = quantile(P,e); % distance that corresponds to the fraction e of rec. points per column
    thresholds = repmat(q,N(1),1); % q has to be applied for each row in d
    % apply individual thresholds
    epsilon = e;
    % apply threshold(s)
    y=double(y<thresholds);
end

% Ahlstrom correction

% loop over the columns of the RP
for j = 1:size(y,2)
    % loop over the lines left in the upper triangle
    thres_flag = false;
    cnt = 1;
    for k = 1 : j-1
        if rem(cnt,rt_thres+1) == 0
            thres_flag = true;
        end        
        if ~y(j-k,j) && ~thres_flag
            cnt = cnt + 1;
        elseif ~y(j-k,j) && thres_flag
            cnt = 1;
            thres_flag = true;
        elseif y(j-k,j) && ~thres_flag
            y(j-k,j) = 0;
            cnt = cnt + 1;           
        elseif y(j-k,j) && thres_flag
            cnt = 1;
            thres_flag = false;
        end
        
    end
    
    % loop over the lines left in the lower triangle
    thres_flag = false;
    cnt = 1;
    for k = 1 :size(y,1) - j
        if rem(cnt,rt_thres+1) == 0
            thres_flag = true;
        end        
        if ~y(j+k,j) && ~thres_flag
            cnt = cnt + 1;
        elseif ~y(j+k,j) && thres_flag
            cnt = 1;
            thres_flag = true;
        elseif y(j+k,j) && ~thres_flag
            y(j+k,j) = 0;
            cnt = cnt + 1;           
        elseif y(j+k,j) && thres_flag
            cnt = 1;
            thres_flag = false;
        end
        
    end
end

        
        

end


