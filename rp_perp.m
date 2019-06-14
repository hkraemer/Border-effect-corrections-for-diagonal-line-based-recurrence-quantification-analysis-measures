function [perp_r,sp,r,epsilon] = rp_perp(varargin)
% RP_PERP    Calculates the perpendicular recurrence plot after Choi et al.
% 1999 and modified by the authors of this script.
%
% Minimum input-arguments : 1
% Maximum input-arguments : 6
%
%    [RP_perp,dot_prod_matrix,RP_normal,epsilon]=rp_perp(Y,E,thres_meth,w,tau,norm)
%
%    calculates the perpendicular recurrence plot 'RP_perp' (also the 
%    conventional recurrence plot 'RP_normal') from an phase space vector 
%    'Y' using the threshold selection method 'thres_meth' under the recurrence
%    threshold 'E'. Distance computations carried out using a specified
%    'norm'. In case you choose the 'var'-fixed threshold selection
%    method, the optional output 'epsilon' will give the actual calculated
%    value. Output variable 'dot_prod_matrix' is matrix, that contains all
%    normalized dot products of the tangential of each reference point 
%    (formed by using 'tau'/2 proceding and subsequent phase space points;
%    Default is tau = 2) to all distance vectors. 
%
%
% Input:
%
% 'Y'                       is a N-by-M matrix corresponding to N time points
%                           and M embedding dimensions.
% 'w' (optional)            is a threshold parameter which allows two points
%                           in phase space to be considered perpendicular,
%                           when the angle between then tangential of the
%                           reference point and the difference vector is
%                           cos(90)+-w. Floating number in the interval [0,1]
%                           in radians. Default is w = 0.025.
% 'norm' (optional)         norm for distance calculation in phasespace to
%                           'euc' (euclidic) or 'max' (maximum). Default is 
%                           max norm.
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
%    Example (CRP toolbox needs to be installed):
%      x = sin(linspace(0,5*2*pi,1000));
%      xe = embed(x,2,50);
%      [r2,~,r,~] = rp_perp(xe,.2,'var',0.95);
%      figure
%      subplot(1,2,1)
%      imagesc(r), colormap([1 1 1; 0 0 0]), axis xy square
%      title('input RP')
%      subplot(1,2,2)
%      imagesc(r2), colormap([1 1 1; 0 0 0]), axis xy square
%      title('perpendicular RP')
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
e = 1; % recurrence threshold
w = .025; % angle threshold
tau = 2; % distance between points at phase space vector for constructing the tangential trajectory

%% get input arguments
% embedding vector
x = varargin{1};
N = size(x); % size of the embedding vector
if N(1) < N(2)
   error('Embedding dimension is larger than the length of the vector. Please check!')
end

% set threshold value
if nargin > 1
    if isa(varargin{2},'double')
        e = varargin{2};
    else
        warning('Threshold has to be numeric.')
    end
end

% set threshold selection method
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

% set angle threshold value
if nargin > 3
    if isa(varargin{4},'double')
        w = varargin{4};
    else
        warning('Threshold has to be numeric.')
    end
end

% set time delay
if nargin > 4
    if isa(varargin{5},'double')
        tau = varargin{5};
    else
        warning('Time delay has to be numeric.')
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
% two estimation variants:
% use reference point and another point in the future, tau time steps ahead
% (after Horai et al, 2002)
tangential_vec = zeros(size(x));
tangential_vec(1+floor(tau/2):end-ceil(tau/2),:) = x(1+tau:end,:) - x(1:end-tau,:);


%% calculate dot product
sp = zeros(N(1),N(1));
for i = 1:N
   % create distance vector between reference point and potential neighbours
   dist_vect = x - repmat(x(i,:),N(1),1);
   % dot product
   sp(i,:) = abs(dot(dist_vect, repmat(tangential_vec(i,:),N(1),1),2))';
   % normalize
   sp(i,:) = sp(i,:) ./ (vecnorm(dist_vect,2,2) .* vecnorm(repmat(tangential_vec(i,:),N(1),1),2,2))';
end

% apply threshold to dot product matrix to check whether perpendicular
perp_r = r .* (sp < w);

% add LOI
perp_r = perp_r + eye(size(perp_r));

end
      

