function [y,P,epsilon] = rp_LM2P(varargin)
%RP_LM2P local minima thresholded recurrence plot with 2 parameters after 
%Schultz et al., Int.J.of Bifurc. Vol 21, 2011 and Wendi & Marwan, Chaos 28, 2018
%
% Minimum input-arguments : 1
% Maximum input-arguments : 6
%
%    [R,DM,epsilon] = rp_LM2P(Y,E,thres_meth,tau,norm,algorithm) 
%    
%    Calculates a recurrence plot R from a phase space vector Y using 
%    the recurrence threshold 'E'. It is constructed from the adjacency- or
%    distancematrix 'DM' (optional output) by looking for the local minima 
%    in each column which fall under the threshold 'E'. In case you choose 
%    the 'var'-fixed threshold selection method optional output 'epsilon' 
%    will give the actual calculated value. The additional parameter 'tau'
%    gives the minimum distance between consecutive minima in a column.
%
% Input:
%
% 'Y'                       is a N-by-M matrix corresponding to N time points
%                           and M embedding dimensions.
%
% 'tau'                     a positive integer, which determines the
%                           minimum tolerable distance (in sampling
%                           units) between consecutive minima in each
%                           column of the distance matrix. Default is tau =
%                           2.
%
% 'norm' (optional)         norm for distance calculation in phasespace to
%                           'euc' (euclidic) or 'max' (maximum). Default is 
%                           max norm.
%
% 'algorithm' (optional)    specify the way of calculating the distance
%                           matrix here. You can choose from
%                           ['loops','vector','matlabvector']. Default is
%                           'vector'.
%
% 'thres_meth' (optional) specifies how the threshold epsilon will
% be calculated. There are three options. Set 'thres_meth' to
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
%    Example (CRP toolbox needs to be installed and MATLABs Signal Processing
%    toolbox is required):
%      x = sin(linspace(0,5*2*pi,1000));
%      xe = embed(x,2,50);
%      r = rp(xe,.2);
%      [r2,~,~] = rp_LM2P(xe,.2,'fix',2);
%      figure
%      subplot(1,2,1)
%      imagesc(r), colormap([1 1 1; 0 0 0]), axis xy square
%      title('input RP')
%      subplot(1,2,2)
%      imagesc(r2), colormap([1 1 1; 0 0 0]), axis xy square
%      title('LM2P')
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
       algorithm = 'vector';
    end
catch
    algorithm = 'vector';
end


try
    tau = varargin{4};
    if rem(tau,1)~=0 || tau < 1 || tau >= N(1)
        warning('recurrence time threshold needs to be a positive integer with maximum size of the RP. Now set to tau = 2')
        tau = 2;
    end
catch
    tau = 2;
end

methLib={'euc','max'}; % the possible norms
try
    meth = varargin{5};
    if ~isa(meth,'char') || ~ismember(meth,methLib)
       warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
       meth = 'max';
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
        thres = 'fix';
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
    % reverse distances in distance matrix
    PP = -P;
    PPP = 10^10*ones(size(PP));
    % look for localmaxima with a peak-to-peak-distance of tau
    for j = 1:size(P,2)
        [~,locs] = findpeaks(PP(:,j),'MinPeakDistance',tau);
        % store reversed distances in final local minima matrix PPP and
        % threshold it
        PPP(locs,j) = P(locs,j);
        PPP(:,j) = double(PPP(:,j)<e);
    end    
    y = PPP;
    epsilon = e;
    
elseif strcmp(thres,'var')
    % get lower (e*100)%-quantile of distance-distribution
    epsilon = quantile(P(:),e);
    % reverse distances in distance matrix
    PP = -P;
    PPP = 10^10*ones(size(PP));
    % look for localmaxima with a peak-to-peak-distance of tau
    for j = 1:size(P,2)
        [~,locs] = findpeaks(PP(:,j),'MinPeakDistance',tau);
        % store reversed distances in final local minima matrix PPP and
        % threshold it        
        PPP(locs,j) = P(locs,j);
        PPP(:,j) = double(PPP(:,j)<epsilon);
    end  
    y = PPP;
    
elseif strcmp(thres,'fan')
    % compute variable threshold for each point in order to get fixed
    % number of nearest neighbours
    q = quantile(P,e); % distance that corresponds to the fraction e of rec. points per column
    thresholds = repmat(q,N(1),1); % q has to be applied for each row in d
    % apply individual thresholds
    epsilon = e;
    % reverse distances in distance matrix
    PP = -P;
    PPP = 10^10*ones(size(PP));
    % look for localmaxima with a peak-to-peak-distance of tau
    for j = 1:size(P,2)
        [~,locs] = findpeaks(PP(:,j),'MinPeakDistance',tau);
        % store reversed distances in final local minima matrix PPP and
        % threshold it
        PPP(locs,j) = P(locs,j);
        PPP(:,j) = double(PPP(:,j)<thresholds(:,j));
    end 
    
    y = PPP;

end 

end
