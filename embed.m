function y = embed(varargin)
% EMBED   Create embedding vector.
%     Y = EMBED(X, M, T) creates embedding vector Y from
%     time series X using a time delay embedding with 
%     dimension M and time delay T. The resulting embedding 
%     vector has length N-T*(M-1), where N is the length 
%     of the original time series.
%
%     Example: x = sin(0:0.1:10*2*pi);
%              y = embed(x,2,16);
%              plot(y(:,1),y(:,2))

% Copyright (c) 2021-
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% $Date:  $
% $Revision: $
%
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% check input arguments
error(nargchk(1,3,nargin))
error(nargoutchk(0,1,nargout))

if nargin < 3
    tau = 1;
else
    tau = varargin{3};
end
if nargin < 2
    m = 2;
else
    m = varargin{2};
end

x = varargin{1}(:);

N = length(x) - (m-1)*tau;
y = buffer(x,N,N-tau,'nodelay');
