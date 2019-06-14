function [iso_r,r,sp,epsilon] = rp_iso(varargin)
% RP_ISO    Calculates the isodirectional recurrence plot after Horai et
% al. 2002
%
% Minimum input-arguments : 1
% Maximum input-arguments : 5
%
%    [RP_iso,RP_normal,sp,epsilon]=rp_iso(Y,E1,E2,thres_meth,tau,norm)
%
%    computes the isodirectional recurrence plot 'RP_iso' (also the 
%    conventional recurrence plot 'RP_normal') from an embedding/phase space
%    vector 'Y' using the threshold selection method 'thres_meth' under the 
%    recurrence threshold 'E1'. Distance computations are carried out using 
%    a specified 'norm'. In case you choose the 'var'-fixed threshold selection
%    method, the optional output 'epsilon' will give the actual calculated
%    value. Output variable 'sp' is a matrix, that contains all the pairwise 
%    distances of the lines set at each point in phase space, formed looking
%    'tau' time steps ahead from the reference point (Default is tau = 2). 
%    The second recurrence threshold 'E2' is applied to 'sp'.
%
%
% Input:
%
% 'Y'                       is a N-by-M matrix corresponding to N time points
%                           and M embedding dimensions.
%
% 'tau'                     determines the length and direction of the
%                           diagonal lines at each reference point, by
%                           looking tau-timesteps ahead and building the
%                           difference vector. Default is tau = 2.
%
% 'E2'                      is the recurrence threshold applied to all the 
%                           pairwise distances of the lines set at each point
%                           in phase space, formed looking 'tau' time steps 
%                           ahead from the reference point.
%
% 'norm' (optional)         norm for distance calculation in phasespace to
%                           'euc' (euclidic) or 'max' (maximum). Default is 
%                           max norm.
%
% 'thres_meth' (optional) specifies how the threshold epsilon will
% be calculated. There are three options. Set 'thres_meth' to
%   - 'fix' The RP is computed under a fixed threshold epsilon specified by
%           input parameter 'E1'.
%   - 'var' The RP is computed under a fixed threshold epsilon, which
%           corresponds to the lower 'E1'-quantile (specified by input parameter
%           'E1') of the distance distribution of all points in phasespace.
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
%      [r2,r,~,~] = rp_iso(xe,.3,.01,'fix');
%      figure
%      subplot(1,2,1)
%      imagesc(r), colormap([1 1 1; 0 0 0]), axis xy square
%      title('input RP')
%      subplot(1,2,2)
%      imagesc(r2), colormap([1 1 1; 0 0 0]), axis xy square
%      title('isodirectional RP')
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
nargoutchk(1,4)

%% set default values for input parameters

tau = 2; % distance between points at phase space vector for constructing the tangential trajectory

%% get input arguments
% embedding vector
x = varargin{1};
N = size(x); % size of the embedding vector
if N(1) < N(2)
   error('Embedding dimension is larger than the length of the vector. Please check!')
end

% set threshold value 1
if nargin > 1
    if isa(varargin{2},'double')
        e = varargin{2};
    else
        warning('Threshold has to be numeric.')
    end
end

% set threshold value 2
if nargin > 2
    if isa(varargin{3},'double')
        e2 = varargin{3};
    else
        warning('Threshold has to be numeric.')
    end
end

% set threshold selection method
thresLib={'fix','var','fan'}; % the possible ways of threshold computation
try
    thres = varargin{4};
    if ~isa(thres,'char') || ~ismember(thres,thresLib)
       warning(['Specified way of calculating threshold should be one of the following possible values:',...
                                10,sprintf('''%s'' ',thresLib{:})])
    end
catch
    thres = 'fix';
end

% set time delay 
if nargin > 4
    if isa(varargin{5},'double')
        tau = varargin{5};
    else
        warning('Time delay has to be numeric.')
    end
    if rem(tau,1)~=0 || tau < 1
        warning('Time delay has to be a positive integer and is now set to the Default tau=2.')
    end
end

% set norm
methLib={'euc','max'}; % the possible norms
try
    meth = varargin{6};
    if ~isa(meth,'char') || ~ismember(meth,methLib)
       warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
    end
catch
    meth = 'max';
end




%% calculate distance matrix D and RP
% distance matrix using MATLAB's pdist function
switch lower(meth)
  case 'euc'
      d = squareform(pdist(x));
  otherwise
      d = squareform(pdist(x,'chebychev'));
end

% apply threshold and get the RP

if strcmp(thres,'fix')
    % apply threshold
    r = double(d < e);
    epsilon = e;
    
elseif strcmp(thres,'var')    
    % get lower (e*100)%-quantile of distance-distribution
    epsilon = quantile(d(:),e);
    r = double(d < epsilon);    

elseif strcmp(thres,'fan')
    % compute variable threshold for each point in order to get fixed
    % number of nearest neighbours
    q = quantile(d,e); % distance that corresponds to the fraction e of rec. points per column
    thresholds = repmat(q,N(1),1); % q has to be applied for each row in d
    % apply individual thresholds
    epsilon = e;
    % apply threshold(s)
    r=double(d<thresholds);
end


%% estimate the tangential vector
% estimation variant:
% use reference point and another point in the future, tau time steps ahead
% (after Horai et al, 2002)

tangential_vec = zeros(size(x));
tangential_vec(1:end-tau,:) = x(1+tau:end,:) - x(1:end-tau,:);


%% calculate distance between the "tangential" vectors
sp = zeros(N(1),N(1));
for i = 1:N   
    % distances between tangential vectors
    switch lower(meth)
        case 'euc'
            sp(i,:) = sqrt(sum((tangential_vec - repmat(tangential_vec(i,:),N(1),1)).^ 2,2));
        otherwise
            dis = abs(tangential_vec - repmat(tangential_vec(i,:),N(1),1));
            sp(i,:) = max(dis,[],2);
    end   
end

% apply second threshold to sp matrix
iso_r = r .* (sp < e2);

   
      
end