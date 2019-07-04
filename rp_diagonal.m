function [X_new,dl_new] = rp_diagonal(varargin)
% RP_DIAGONAL  RP with corrected lines ...
% 
%    [RP_new, dl_new] = rp_diagonal(RP) 
%    computes a new recurrence plot 'RP_new' by altering diagonal line structures 
%    in the input recurrence plot 'RP': Slubs, but also block structures, 
%    are deleted in favour of the longest diagonal lines (skeletonization). 
%    Whenever a diagonal line (starting with the longest lines contained in 
%    the diagonal line length histogram) encounters an adjacent diagonal 
%    line, this adjacent line and - recursively - all its consecutive 
%    adjacent lines, get deleted. 
%
%    Output:
%    You receive the corrected recurrence plot 'RP_new' INCLUDING the line 
%    of identity. Optional you also receive a matrix containing all lines 
%    of this new RP. This matrix has three lines and as many columns as 
%    there are lines in the RP. In the first line the total length of each 
%    line is stored. In the second and third line, the corresponding line 
%    and column indices are stored.
%
%    Example (CRP toolbox needs to be installed):
%      x = sin(linspace(0,5*2*pi,1000));
%      xe = embed(x,2,50);
%      r = rp(xe,.2);
%      [r2, ~] = rp_diagonal(r);
%      figure
%      subplot(1,2,1)
%      imagesc(r), colormap([1 1 1; 0 0 0]), axis xy square
%      title('input RP')
%      subplot(1,2,2)
%      imagesc(r2), colormap([1 1 1; 0 0 0]), axis xy square
%      title('diagonal RP')
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
narginchk(1,3)
nargoutchk(1,2)

X = varargin{1};

% size of the RP
[N_org,M_org] = size(X);
if N_org~=M_org
    error('Input needs to be a squared, binary recurrence matrix')
end
if sum(sum(X>1))~=0 || sum(sum(X<0))~=0 || sum(sum(rem(X,1)))~=0
    error('Input needs to be a squared, binary recurrence matrix')
end

% check whether input RP is symmetric
if issymmetric(X)
    symm = true;
    % if yes, just take the lower triangle
    X2 = tril(X);
    
    % convert this RP into a close returns map and just use the upper half
    X_cl = convertRP(X2);
 
    % get line distributions
    [~, lines_1] = dl_h(X_cl); % black horizontal lines
        
    % make a copy of the line matrix
    lines_1_copy = lines_1;
    
    Nlines = size(lines_1,2); % number of found lines
    
    % create a close returns map with horizontal lines represented by 
    % numbers, equal to its lengths
    X_hori = zeros(size(X_cl));
    for i = 1:size(lines_1,2)
        line_ind = lines_1(2,i);
        column_ind = lines_1(3,i);
        for j = 0:lines_1(1,i)-1
            X_hori(line_ind,column_ind+j) = lines_1(1,i);
        end
    end
      
    
else
    symm = false;
    % if not, store lower triangle in X2 and upper triangle transposed in
    % X3
    X2 = tril(X);
    X3 = triu(X)';
    
    % convert these RPs into close returns maps and just use the upper half
    X_cl = convertRP(X2);
    X_cl2 = convertRP(X3);
    
    % get line distributions
    [~, lines_1] = dl_h(X_cl); % black horizontal lines
    [~, lines_2] = dl_h(X_cl2); % black horizontal lines
    
    % make a copy of the line matrices
    lines_1_copy = lines_1;
    lines_2_copy = lines_2;
    
    Nlines = size(lines_1,2); % number of found lines in lower triangle
    Nlines2 = size(lines_2,2); % number of found lines in upper triangle
    
    % create a close returns map with horizontal lines represented by 
    % numbers, equal to its lengths
    X_hori = zeros(size(X_cl));
    for i = 1:size(lines_1,2)
        line_ind = lines_1(2,i);
        column_ind = lines_1(3,i);
        for j = 0:lines_1(1,i)-1
            X_hori(line_ind,column_ind+j) = lines_1(1,i);
        end
    end
    X_hori2 = zeros(size(X_cl2));
    for i = 1:size(lines_2,2)
        line_ind = lines_2(2,i);
        column_ind = lines_2(3,i);
        for j = 0:lines_2(1,i)-1
            X_hori2(line_ind,column_ind+j) = lines_2(1,i);
        end
    end
end


% scan the lines, start with the longest one and discard all adjacent lines

% initialize final line matrix
line_matrix_final = zeros(3,1);

% go through all lines stored in the sorted line matrix
[N,M] = size(X_hori);
for l_ind = 1:Nlines

    % check if line is still in the rendered line matrix
    if ~ismember(lines_1(:,l_ind)',lines_1_copy','rows')
        continue
    end
    
    % get index pair for start of the line
    linei = lines_1(2,l_ind); 
    columni = lines_1(3,l_ind);
    
    % copy this line in the final line matrix
    line_matrix_final = horzcat(line_matrix_final,lines_1(:,l_ind));
    
    % delete this line from the RP
    X_hori = delete_line_from_RP(X_hori,lines_1(:,l_ind));
 
    % go along each point of the line and check for neighbours
    l_max = lines_1(1,l_ind);
    for l = 1:l_max
        
        % scan each line twice - above and underneth
        for index = -1:2:1
            
            % make sure not to exceed RP-boundaries
            if linei+index > N | linei+index == 0
                break
            end
                           
            % if there is a neighbouring point, call recursive scan-function
            if X_hori(linei+index,columni+l-1)

                [X_hori,lines_1_copy]=scan_lines(X_hori,lines_1_copy,...
                    linei+index,columni+l-1);

            end

        end
        
    end

end

% if not symmetric input RP, than compute for the upper triangle as well
if ~symm
    % initialize final line matrix
    line_matrix_final2 = zeros(3,1);
    for l_ind = 1:Nlines2

        % check if line is still in the rendered line matrix
        if ~ismember(lines_2(:,l_ind)',lines_2_copy','rows')
            continue
        end

        % get index pair for start of the line
        linei = lines_2(2,l_ind); 
        columni = lines_2(3,l_ind);

        % copy this line in the final line matrix
        line_matrix_final2 = horzcat(line_matrix_final2,lines_2(:,l_ind));

        % delete this line from the RP
        X_hori2 = delete_line_from_RP(X_hori2,lines_2(:,l_ind));

        % go along each point of the line and check for neighbours
        l_max = lines_2(1,l_ind);
        for l = 1:l_max

            % scan each line twice - above and underneth
            for scan = 1:2

                if scan == 1
                    index = 1;
                    % make sure not to exceed RP-boundaries
                    if linei+index > N
                        break
                    end
                else
                    index = -1;
                    % make sure not to exceed RP-boundaries
                    if linei+index == 0
                        break
                    end 
                end

                % if there is a neighbouring point, call recursive scan-function
                if X_hori2(linei+index,columni+l-1)

                    [X_hori2,lines_2_copy]=scan_lines(X_hori2,lines_2_copy,...
                        linei+index,columni+l-1);

                end

            end

        end      

    end
end

% build RP based on the histogramm of the reduced lines

X_cl_new = zeros(N,M);
if ~symm
   X_cl2_new = zeros(N,M); 
end

% fill up close returns map with lines stored in the new line matrix
for i = 1:size(line_matrix_final,2)
    l_max = line_matrix_final(1,i);
    linei = line_matrix_final(2,i);
    columni = line_matrix_final(3,i);
    for j = 1:l_max
        X_cl_new(linei,columni+j-1) = 1;
    end
end

if symm
    % revert this close returns map into a legal RP
    XX = revertRP(X_cl_new);    
    X_new = XX + (XX-eye(size(XX)))';
else 
    % fill up close returns map with lines stored in the new line matrix
    for i = 1:size(line_matrix_final2,2)
        l_max = line_matrix_final2(1,i);
        linei = line_matrix_final2(2,i);
        columni = line_matrix_final2(3,i);
        for j = 1:l_max
            X_cl2_new(linei,columni+j-1) = 1;
        end
    end
    % revert this close returns map into a legal RP
    XX = revertRP(X_cl_new);
    XXX= revertRP(X_cl2_new);
    X_new = XX + (XXX-eye(size(XXX)))';
    
end

% bind optional output

% get line distributions of new RP
if nargout > 1
    [~, lines3] = dl_e(X_new);
    dl_new = sortrows(lines3','descend')';
end

end

%%%%% Helper functions %%%%%%
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
%    Examples: X = crp(rand(200,1),1,1,.3,'fan','silent');
%              [l l_dist] = dl_e(X);
%              hist(l_dist(1,:),200)
%
%    See also CRQA, TT, DL.

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
end

function lines = getLinesOnDiag(M,j)
% getLinesOnDiag computes lines on a diagonal 'j' of a matrix 'M'. The
% output is a 3-by-number_of_lines_in_the_digonal-matrix. In the first line
% the length of the line is stored and in line 2 and 3 the respective line
% and column index of the starting point of that diagonal.
%
% Copyright (c) 2019
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
% Nikolaus Koopmann
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.
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
    
end


function tuples = get_indices(indices,position_from_LOI)
% GET_INDICES gets "true" indices from of the lines represented in the column 
% vector 'indices' and its diagonal determined by 'position_from_LOI'.
% "True" indices in this case means the indices in the corresponding RP,
% not the close returns map.
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


end

function [a_out, b_out]=dl_h(x)
% DL_H   Mean of the horizontal line lengths and their distribution in a
% close returns map, additionally with the corresponding indices of the lines.
%    A=DL_H(X) computes the mean of the length of the horizontal 
%    line structures in a close returns map of a recurrence plot.
%
%    [A B]=DL_H(X) computes the mean A and the lengths of the
%    found horizontal lines, stored in the first line of B. B is a 3 line
%    matrix storing the found horizontal line lengths in its columns. Line 2
%    and 3 store the indices i, j of the startpoint of the horizontal line
%    stored in the same column in the first line.
%    In order to get the 
%    histogramme of the line lengths, simply call 
%    HIST(B(1,:),[1 MAX(B(1,:))]).
%
%    Examples: X = crp(rand(200,1),1,1,.3,'fan','silent');
%              [l l_dist] = dl_h(convertRP(X));
%              hist(l_dist(1,:),200)
%
%    See also CRQA, TT, DL, convertRP, revertRP

% Copyright (c) 2018-
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


narginchk(1,1)
nargoutchk(0,2)

[N,~] = size(x);

liness = zeros(3,1);
for j = 1:N
    d = x(j,:)';
    starts = find(diff([0; d],1)==1);
    ends = find(diff([d; 0],1)==-1);
    
    if ~isempty(starts)
        lines = zeros(3,numel(starts));
        for n=1:numel(starts)
            lines(2,n) = j;
            lines(3,n) = starts(n);
            lines(1,n) = ends(n) - starts(n) +1;        
        end
    else
        lines = zeros(3,1);
    end
    liness = horzcat(liness,lines);
end

% remove lines of length zero (=no line)
zero_lines = liness(1,:)==0;
liness(:,zero_lines) = []; 

b_out= sortrows(liness','descend')';
a_out = mean(b_out(1,:));
end

function Y = convertRP(X)
% CONVERTRP   Transforms the standard RP to a close returns map
%    Y = convertRP(X)
%
%    Example:
%      x = sin(linspace(0,5*2*pi,200));
%      xe = embed(x,2,5);
%      r = rp(xe,.2);
%      imagesc(r)
%      c = convertRP(r);
%      imagesc(c)

%% size of RP
N = size(X);

%% initialize new matrix
Y = zeros(2*N(1)+1,N(1));

%% fill rows of Y by the diagonals of X
% upper triangle
for i = 0:N(1)-1
   Y(N(1)+i+1,(1:(N(1)-i))) = diag(X,i);
end
   
% lower triangle
for i = 0:N(1)-1
   Y(N(1)-i+1,(1:(N(1)-i))+i) = diag(X,-i);
end

end

function X = revertRP(Y)
% REVERTRP   Transforms a close returns map to the standard RP
%    X = revertRP(Y)
%
%    Example:
%      x = sin(linspace(0,5*2*pi,200));
%      xe = embed(x,2,5);
%      r = rp(xe,.2);
%      c = convertRP(r);
%      r2 = revertRP(c);
%      imagesc(r2)
%      imagesc(r-r2) % shows the difference between original and reverted

%% size of close returns map
N = size(Y);

%% initialize new matrix
X = zeros(N(2),N(2));

%% make Y to a square matrix, fill the new part with zeros
Z = [Y zeros(N(1),N(2)+1)];
Z = flipud(Z); % flip upside down

%% fill columns of  by the diagonals of Z (but only the first N points) 
for i = 1:N(2)
    di = diag(Z,-i);
    X(:,N(2)-i+1) = di(1:N(2));
end
end


function RP = delete_line_from_RP(RP,l_vec)
% deletes a line, specified in 'l_vec' (line vector, with first line being
% the total line length, the second line the line-index of the starting point
% and the third line the column-index of the starting pint) from the 'RP'.
    
    RP(l_vec(2),l_vec(3)+(1:l_vec(1))-1) = 0;
%    X = RP;
end


% function [XX,YY]= scan_lines(X,l_vec,line,column)
%     
%     % for the input index tuple look for the start indices
%     index = 0;
%     while true        
%         % check whether the input index tuple is a listed index for starting
%         % points of line lengths in the line matrix 
%         loc_line = find(line==l_vec(2,:));
%         del_ind = loc_line(column+index==l_vec(3,loc_line));
%         
%         if del_ind
%             break
%         else
%             index = index - 1;
%         end
%     end   
%     
%     % delete the line from RP
%     XX = delete_line_from_RP(X,l_vec(:,del_ind));
%     
%     % bind line length, line & column starting index
%     len = l_vec(1,del_ind);
%     li = l_vec(2,del_ind);
%     co = l_vec(3,del_ind);
%     
%     % delete the line from the line matix
%     l_vec(:,del_ind) = [];
%     
%     [N,~] = size(XX);
%     
%     %%%%%% check for borders of the RP %%%%%%
%     if li-1 < 1
%         flag1 = false;
%     else
%         flag1 = true;
%     end
%     if li+1 > N
%         flag2 = false;
%     else
%         flag2 = true;
%     end
%     
%     for i = 1:len
%     
%         % check above the line for a neighbour
%         if flag1 && XX(li-1,co+i-1)
%             % call function itself
%             [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i-1);
%             
%         % check underneeth the line for a neighbour    
%         elseif flag2 && XX(li+1,co+i-1)
%             % call function itself
%             [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i-1);
%         end
%     end
%     
%     YY = l_vec;
% 
% end

function [XX,YY]= scan_lines(XX,l_vec,line,column)
    
    % for the input index tuple look for the start indices
    index = 0;
    while true        
        % check whether the input index tuple is a listed index for starting
        % points of line lengths in the line matrix 
        loc_line = find(line==l_vec(2,:));
        del_ind = loc_line(column+index==l_vec(3,loc_line));
        
        if del_ind
            break
        else
            index = index - 1;
        end
    end   
    
    % delete the line from RP
    %XX = delete_line_from_RP(XX,l_vec(:,del_ind));
    XX(l_vec(2,del_ind),l_vec(3,del_ind)+(1:l_vec(1,del_ind))-1) = 0;
    
    % bind line length, line & column starting index
    len = l_vec(1,del_ind);
    li = l_vec(2,del_ind);
    co = l_vec(3,del_ind);
    
    % delete the line from the line matix
    l_vec(:,del_ind) = [];
    
    [N,M] = size(XX);
    
    %%%%%% check for borders of the RP %%%%%%
    if li-1 < 1
        flag1 = false;
    else
        flag1 = true;
    end
    if li+1 > N
        flag2 = false;
    else
        flag2 = true;
    end
    
    for i = 1:len
        
        % check for borders of the RP
        if li-1 < 1 || co+i-2 == 0
            flag1b = false;
        else
            flag1b = true;
        end
        if li-1 < 1 || co+i > M
            flag1c = false;
        else
            flag1c = true;
        end
        if li+1 > N || co+i-2 == 0
            flag2b = false;
        else
            flag2b = true;
        end
        if li+1 > N || co+i > M
            flag2c = false;
        else
            flag2c = true;
        end
        
        % check above left for a neighbour
        if flag1b && XX(li-1,co+i-2)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i-2);
            
        % check above the line for a neighbour
        elseif flag1 && XX(li-1,co+i-1)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i-1);
            
        % check above right for a neighbour
        elseif flag1c && XX(li-1,co+i)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i);
            
        % check underneeth left for a neighbour    
        elseif flag2b && XX(li+1,co+i-2)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i-2);
            
        % check underneeth the line for a neighbour    
        elseif flag2 && XX(li+1,co+i-1)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i-1);
            
        % check underneeth right for a neighbour    
        elseif flag2c && XX(li+1,co+i)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i);
        end
    end
    
    YY = l_vec;

end

function y = embed(varargin)
%
%   version 1.0
%
%
%     Create embedding vector using time delay embedding
%     Y=EMBED(X,M,T) create the embedding vector Y from the time
%     series X using a time delay embedding with dimension M and
%     delay T. The resulting embedding vector has length N-T*(M-1),
%     where N is the length of the original time series.
%
%     Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         y = embed(x,2,17);
%         plot(y(:,1),y(:,2))
%
% Copyright (c) 2012
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

narginchk(1,3)
nargoutchk(0,1)

try
    t = varargin{3};
catch
        t = 1;
end
        
try
    m = varargin{2};
catch
        m = 1;
end    

x = varargin{1}(:);

% length of time series
Nx = length(x);
NX = Nx-t*(m-1);
if t*(m-1) > Nx
   warning('embedding timespan exceeding length of time series')
end

%% create embeeding vector using Matlab's vectorization
% index series considering the time delay and dimension
for mi = 1:m
    jx(1+NX*(mi-1):NX+NX*(mi-1)) = 1+t*(mi-1):NX+t*(mi-1);
end

% the final embedding vector
y = reshape(x(jx),NX,m);
end