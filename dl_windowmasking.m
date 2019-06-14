function [a_out, b_out] = dl_windowmasking(x)
% DL_WINDOWMASKING   Mean of the diagonal line lengths and their distribution
% additionally with the corresponding indices of the lines of a 45 degree 
% rotated RP.
%    A=dl_windowmasking(X) computes the mean of the length of the diagonal 
%    line structures in a 45deg-rotated recurrence plot 'X'.
%
%    [A B]=dl_windowmasking(X) computes the mean A and the lengths of the
%    found diagonal lines of a 45deg-roated RP 'X', stored in the first line 
%    of B. B is a 3 line matrix storing the found diagonal line lengths in 
%    its columns. Line 2 and 3 store the indices i, j of the startpoint of
%    the diagonal line stored in the same column in the first line.
%    In order to get the histogramme of the line lengths, simply call 
%    HIST(B(1,:),[1 MAX(B(1,:))]).
%
%    Examples (CRP toolbox needs to be installed):
%       x = sin(linspace(0,5*2*pi,1050));
%       xe = embed(x,2,50);
%       r = rp(xe,.2);
%       [l l_dist] = dl_windowmasking(r);
%       subplot(1,2,1)
%       imagesc(r), colormap([1 1 1;0 0 0]), axis xy square
%       title('underlying RP')
%       subplot(1,2,2)
%       histogram(l_dist(1,:),1000)
%       xlim([0 1000])
%       xlabel('diagonal line length')
%       ylabel('counts')
%       title('diagonal line length histogram - window masking')
%
% Donath, J. 2016 - Bachelor thesis, Untersuchung alternativer Zeitfenster-
% formen fuer die quantitative Rekurrenzanalyse anhand von Modellsystemen 
% und Klimadaten, HU Berlin
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


narginchk(1,1)
nargoutchk(0,2)

% transform RP into a different shape
X = shape3(x);


[a_out, b_out] = dl_e(X);



function [a_out, b_out] = dl_e(X)
% DL_E   Mean of the diagonal line lengths and their distribution
% additionally with the corresponding indices of the lines.
%    A=DL_E(X) computes the mean of the length of the diagonal 
%    line structures in a recurrence plot.
%
%    [A B]=DL_E(X) computes the mean A and the lengths of the
%    found diagonal lines, stored in the first line of B. B is a 3 line
%    matrix storing the found diagonal line lengths in its columns. Line 2
%    and 3 store the indices i, j of the startpoint of the diagonal line
%    stored in the same column in the first line.
%    In order to get the 
%    histogramme of the line lengths, simply call 
%    HIST(B(1,:),[1 MAX(B(1,:))]).
%
% Copyright (c) 2019-
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
%
% $Date: 2018/09/19 $
% $Revision:  $
%
% $Log: dl.m,v $
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.
[Y,~] = size(X);
lines(:,1) = getLinesOnDiag(X,-Y+1); % init with first (scalar) diagonal
for j=-Y+2:Y-1
    lines = horzcat(lines,getLinesOnDiag(X,j)); 
end

% remove lines of length zero (=no line)
zero_lines = lines(1,:)==0;
lines(:,zero_lines) = []; 

b_out= sortrows(lines','descend')';
a_out = mean(b_out(1,:));


function lines = getLinesOnDiag(M,j)
    d = diag(M,j);
    if ~any(d)
        lines = [0;0;0];
        return
    end
    starts = find(diff([0; d],1)==1);
    ends = find(diff([d; 0],1)==-1);

    lines = zeros(3,numel(starts));
    for n=1:numel(starts)
        ll = get_indices(starts(n),j);
        for k = 1:length(ll)
            lll = ll{k};
            lines(2,n) = lll(1);
            lines(3,n) = lll(2);
            lines(1,n) = ends(n) - starts(n) +1;
        end
    end
    


function tuples = get_indices(indices,position_from_LOI)
% GET_INDICES gets true indices from the 'sourceRP' of the lines
% represented in the column vector 'indices' and its position in the
% 'sourceRP' is determined by 'position_from_LOI'.
%
% tuples = get_indices(indices,position_from_LOI,sourceRP)
%
% Input:
% 'indices' is a column vector, 'position_from_LOI' a integer
%
% Output: a cell array of tuples conating the true line and column index in
% the sourceRP.
%
% Copyright (c) 2019
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% check input
narginchk(2,2)
nargoutchk(1,1)

if size(indices,1)<size(indices,2)
    indices = indices';
end
if rem(position_from_LOI,1)~=0
    error('position_from_LOI needs to be a integer')
end



%%
tuples = cell(1,length(indices));

for i = 1:length(indices)
    
    if position_from_LOI < 0
        start_line = indices(i) + abs(position_from_LOI);
        start_column = indices(i);
    elseif position_from_LOI > 0
        start_line = indices(i);
        start_column = indices(i) + abs(position_from_LOI); 
    elseif position_from_LOI == 0
        start_line = indices(i);
        start_column = indices(i);  
    end    
   
    tuples{i}=[start_line start_column];
end



function [Y] = shape3(X)
%==========================================================================
%Creates a new window Y with lenght s based on X. All border diagonals have
%the lenght s in Y.
%==========================================================================
s = floor(size(X,1)/2);
sc = floor(size(X,1)/4);
Y = zeros(size(X,1));
for i = 1:s
    for j = 0:sc-1
        Y(sc-j+i,sc+j+i) = X(sc-j+i,sc+j+i);
        Y(sc-j+i,sc+j+i+1) = X(sc-j+i,sc+j+i+1);
        Y(sc+j+i,sc-j+i) = X(sc+j+i,sc-j+i);
        Y(sc+j+i,sc-j+i+1) = X(sc+j+i,sc-j+i+1);
    end
end


function [Y] = shape2(X)
%==========================================================================
%Creates a new window Y with lenght s based on X. All border diagonals have
%the lenght s in Y.
%==========================================================================
s = floor(size(X,1)/2);
sc = floor(size(X,1)/4);
Y1 = tril(X(:,sc:(s+sc)));
Y2 = fliplr(flipud(tril(fliplr(flipud((X(:,sc:(s+sc))))))));
Y = Y1 .* Y2;
