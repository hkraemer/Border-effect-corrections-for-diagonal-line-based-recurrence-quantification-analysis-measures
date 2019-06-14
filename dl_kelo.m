function [a_out, b_out]=dl_kelo(varargin)
% DL_KELO   Mean of the diagonal line lengths and their distribution - 
% Correction for border lines: KEep LOngest diagonal line (kelo)
%    A=dl_kelo(X) computes the mean of the length of the diagonal 
%    line structures in a recurrence plot X. In this correction, just the 
%    longest border line (in each triangle) of the RP is counted. All other 
%    border lines are discarded.
%
%    A=dl_kelo(X,'semi') computes the mean of the length of the diagonal 
%    line structures in a recurrence plot X using the mentionded correction. 
%    Not only lines starting AND ending at a border of the RP, but also semi
%    border lines - lines, that start OR end at a border of the RP - are 
%    denoted as border lines. The longest of these count.
%
%    [A B]=dl_kelo(X,'semi') computes the mean A and the lengths of the
%    found diagonal lines of the recurrence plot X, stored in B, using the 
%    correction mentioned above and also accounts for semi-border diagonals.
%    In order to get the histogramme of the line lengths, simply call 
%    HIST(B,[1 MAX(B)]).
%
%    Examples (CRP toolbox needs to be installed):
%       x = sin(linspace(0,5*2*pi,1050));
%       xe = embed(x,2,50);
%       r = rp(xe,.2);
%       [l l_dist] = dl_kelo(r);
%       subplot(1,2,1)
%       imagesc(r), colormap([1 1 1;0 0 0]), axis xy square
%       title('underlying RP')
%       subplot(1,2,2)
%       histogram(l_dist,1000)
%       xlim([0 1000])
%       xlabel('diagonal line length')
%       ylabel('counts')
%       title('diagonal line length histogram - kelo correction')
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

X = varargin{1};
styleLib={'normal','semi'}; % the possible borderline-style to look for
try
    type = varargin{2};
    if ~isa(type,'char') || ~ismember(type,styleLib)
        warning(['Specified RP type should be one of the following possible values:',...
           10,sprintf('''%s'' ',styleLib{:})])
    end
catch
    type = 'normal';
end

[Y,~] = size(X);
if issymmetric(X)
    [lines(1),borderlines(1)] = getLinesOnDiag(X,-Y+1,type); % init with first (scalar) diagonal
    for j=-Y+2:-1
        [ll,bl] = getLinesOnDiag(X,j,type);
        lines = horzcat(lines,ll);
        borderlines = horzcat(borderlines,bl);
    end
    % append lines for second triangle 
    lines = horzcat(lines,lines);
    % append longest border lines for second triangle (but exclude LOI)
    lines = horzcat(lines,max(borderlines),max(borderlines));
else
    [lines(1),borderlines(1)] = getLinesOnDiag(X,-Y+1,type); % init with first (scalar) diagonal
    for j=-Y+2:Y-1
        [ll,bl] = getLinesOnDiag(X,j,type);
        lines = horzcat(lines,ll);
        borderlines = horzcat(borderlines,bl);
    end
    borderlines = sort(borderlines,'descend');
    % add longest borderlines to lines
    lines = horzcat(lines,borderlines(2:3));
end

% remove lines of length zero (=no line)
zero_lines = lines(:)==0;
lines(zero_lines) = []; 

b_out= sort(lines,'descend')';
a_out = mean(b_out);
end

function [lines, borderline] = getLinesOnDiag(M,j,type)
    d = diag(M,j);
    border_line_length = length(d);
    if ~any(d)
        lines = 0;
        borderline = 0;
        return
    end
    starts = find(diff([0; d],1)==1);
    ends = find(diff([d; 0],1)==-1);

    lines = zeros(1,numel(starts));
    borderline = zeros(1,numel(starts));
    
    if strcmp(type,'normal')
        for n=1:numel(starts)
            if ends(n) - starts(n) + 1 < border_line_length
                lines(n) = ends(n) - starts(n) +1;
            elseif ends(n) - starts(n) + 1 == border_line_length
                borderline = ends(n) - starts(n) +1;
            end
        end
    elseif strcmp(type,'semi')
        for n=1:numel(starts)
            if ends(n) ~= border_line_length && starts(n) ~=1               
                lines(n) = ends(n) - starts(n) +1;
            else
                borderline(n) = ends(n) - starts(n) +1;
            end
        end    
    end
    
end
