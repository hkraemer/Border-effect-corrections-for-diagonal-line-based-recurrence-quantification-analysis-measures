function [y,P,epsilon] = rp(varargin)
%
% Minimum input-arguments : 2
% Maximum input-arguments : 6
%
%    [R,DM,epsilon] = RP(Y,E,thres_calc,norm,type,algorithm) 
%    
%    Calculates a recurrence plot R from an embedding vector Y and using 
%    the threshold 'E'. 'DM' is an optional output and is the adjacency- or
%    distancematrix. In case you choose the 'var'-fixed threshold selection
%    method optional output 'epsilon' will give the actual calculated
%    value.
%
%
% Input:
%
% 'Y'                       is a N-by-M matrix corresponding to N time points
%                           and M embedding dimensions.
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
% 'threshold-calc' (optional) specifies how the threshold epsilon will
% be calculated. There are three options. Set 'threshold-calc' to
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
% 'type' (optional) specifies the type of the RP.
%   - 'normal'      The RP is computed after the definition of Eckmann et
%                   al.1987
%   - 'diagonal'    The RP is computed after the definition of Eckmann et
%                   al.1987 and then line corrected after Kraemer and
%                   Marwan 2019
%   - 'shape'       The RP is computed the definition of Eckmann et
%                   al. 1987 and then shape-converted after J.Donath 2016
%                   (windowshape 3).
%
%    Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         xVec = embed(x,2,17);
%         R = rp(xVec,.1);
%         imagesc(R)

% Copyright (c) 2019
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Modified by Hauke Kr�mer,Potsdam Institute for Climate Impact Research, 
% Germany http://www.pik-potsdam.de
%
% Contact: hkraemer@pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%% check input
narginchk(1,6)
nargoutchk(0,3)

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

typeLib={'normal','diagonal','shape'}; % the possible types
try
    type = varargin{5};
    if ~isa(type,'char') || ~ismember(type,typeLib)
        warning(['Specified RP type should be one of the following possible values:',...
           10,sprintf('''%s'' ',typeLib{:})])
    end
catch
    type = 'normal';
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

x = varargin{1};

N = size(x);
if N(1) < N(2)
   warning('Embedding dimension is larger than the length of the vector. Please check!')
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

if strcmp(type,'diagonal')
    [y, ~] = rp_line_correct(y);
elseif strcmp(type,'shape')
    y = shape3(y);
end

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

end

function [X_new,dl_new] = rp_line_correct(varargin)
% 
%    [RP_new, dl_new] = rp_line_correct(RP) 
%    computes a new recurrence plot by altering diagonal line structures in
%    the input recurrence plot 'RP': Slubs, but also block structures, are 
%    deleted in favour of the longest diagonal lines. Whenever a diagonal 
%    line encounters an adjacent diagonal line of lower lengths, this 
%    adjacent line and all its adjacent smaller lines, get deleted.
%
%    Output:
%    You receive the corrected recurrence plot 'RP_new' INCLUDING the LOI.
%    Optional you also receive a matrix containing all lines of this new
%    RP. This matrix has three lines and as many columns as there are lines
%    in the RP. In the first line the total length of each line is stored.
%    In the second and third line, the corresponding line and column index
%    is stored.
%    
% Copyright (c) 2019-
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
narginchk(1,3)
nargoutchk(1,2)

X = varargin{1};

% size of the RP
[N,M] = size(X);
if N~=M
    error('Input needs to be a squared, binary recurrence plot')
end

% check if LOI is excluded. If not, exclude LOI
if any(diag(X))
    X = double(X) - eye(size(X));
end

% get line distributions
[~, line_matrix] = dl_e(X); % black diagonal lines conventional

% sort the line length / index matrix

[lines_1, b1] = sort(line_matrix,2,'descend');
lines_1(2,:) = line_matrix(2,b1(1,:));
lines_1(3,:) = line_matrix(3,b1(1,:));


% make a copy of the line matrix
lines_1_copy = lines_1;

% scan the lines, start with the longest one and discard adjacent lines
% which are smaller

% go through all lines stored in the sorted line matrix

for l_ind = 1:size(lines_1,2)
    
    % get index pair for start of the line
    line = lines_1(2,l_ind);
    column = lines_1(3,l_ind);
    
    % check if line is still in the rendered line matrix
    if ~ismember(lines_1(:,l_ind)',lines_1_copy','rows')
        continue
    end
    
    % check whether line is in upper or lower triangle of the RP
    if line > column
        upper = 0;
    else
        upper = 1;
    end
 
    % go along each point of the line and check for neighbours
    l_max = lines_1(1,l_ind);
    cnt = 1;
    for l = 1:l_max
        
        if l == 1
            first = true;
        else
            first = false;
        end
        
        if l == l_max
            last = true;
        else
            last = false;
        end
        
        % check underneeth the line
        flag = true;
        index = 1;
        while flag
            
            if ~upper
                
                % LOWER TRIANGLE
                
                % make sure not to exceed RP-boundaries
                if line+index+l-1 > N
                    break
                end

                if X(line+index+l-1,column+l-1)
                    % check whether this index is a listed index for starting
                    % points of line lengths in the line matrix             
                    loc_line = find(line+index+l-1==lines_1_copy(2,:));
                    loc_column = find(column+l-1==lines_1_copy(3,:));

                    % if this is a starting point of a line stored in the line
                    % matrix then check whether that line is longer than l_max
                    % - cnt + 1 than delete this line from the line
                    % matrix
                    if any(ismember(loc_line,loc_column))
                        del_ind = loc_line(ismember(loc_line,loc_column));
                        if lines_1_copy(1,del_ind) <= l_max 
                            lines_1_copy(:,del_ind) = [];
                        end

                    % if this is not the case, then backtrack it until you
                    % reach the starting point
                    else
                        if first
                            index3 = 1;
                            flag3 = true;
                            while flag3
                               if line + index - index3 == 0 || column - index3 == 0
                                   break
                               end
                               if X(line+index-index3,column-index3)
                                   % localize the point in the line matrix                               
                                   loc_line = find(line+index-index3==lines_1_copy(2,:));
                                   loc_column = find(column-index3==lines_1_copy(3,:));
                                   % check if line and column index share the
                                   % same column in the line matrix. 
                                   % If they do, delete this line when it is 
                                   % smaller then the actual line
                                   if any(ismember(loc_line,loc_column))                                     
                                       del_ind = loc_line(ismember(loc_line,loc_column));
                                       if lines_1_copy(1,del_ind) <= l_max
                                           lines_1_copy(:,del_ind) = [];
                                       end
                                   end                  
                               else
                                   flag3 = false;
                               end
                               index3 = index3 + 1;
                            end
                        end
                        
                        
                    end

                else
                    % stop looking for neighbouring lines
                    flag = false;
                end
                index = index + 1;
                
                
            else
                % UPPER TRIANGLE
                
                % make sure not to exceed RP-boundaries
                if column+index+l-1 > N
                    break
                end
                

                if X(line+l-1,column+index+l-1)
                    % check whether this index is a listed index for starting
                    % points of line lengths in the line matrix             
                    loc_line = find(line+l-1==lines_1_copy(2,:));
                    loc_column = find(column+index+l-1==lines_1_copy(3,:));

                    % if this is a starting point of a line stored in the line
                    % matrix then check whether that line is longer than l_max
                    % - cnt + 1 than delete this line from the line
                    % matrix
                    if any(ismember(loc_line,loc_column))
                        del_ind = loc_line(ismember(loc_line,loc_column));
                        if lines_1_copy(1,del_ind) <= l_max 
                            lines_1_copy(:,del_ind) = [];
                        end
                    % if this is not the case, then backtrack it until you
                    % reach the starting point
                    else
                        if first
                            index3 = 1;
                            flag3 = true;
                            while flag3
                               if line - index3 == 0 || column + index - index3 == 0
                                   break
                               end
                               if X(line-index3,column+index-index3)
                                   % localize the point in the line matrix                               
                                   loc_line = find(line-index3==lines_1_copy(2,:));
                                   loc_column = find(column+index-index3==lines_1_copy(3,:));
                                   % check if line and column index share the
                                   % same column in the line matrix. 
                                   % If they do, delete this line when it is 
                                   % smaller then the actual line
                                   if any(ismember(loc_line,loc_column))                                     
                                       del_ind = loc_line(ismember(loc_line,loc_column));
                                       if lines_1_copy(1,del_ind) <= l_max
                                           lines_1_copy(:,del_ind) = [];
                                       end
                                   end                  
                               else
                                   flag3 = false;
                               end
                               index3 = index3 + 1;
                            end
                        end
                        
                        
                    end
                    
                else
                    % stop looking for neighbouring lines
                    flag = false;
                end
                index = index + 1;                
            end   
             
        end
         

        
        % check above the line
        flag = true;
        index = 1;
        while flag
            
            if ~upper
                
                % LOWER TRIANGLE
                
                % make sure not to exceed RP-boundaries
                if line+l-1-index == 0
                    break
                end              

                if X(line-index+l-1,column+l-1)

                    % check whether this index is a listed index for starting
                    % points of line lengths in the line matrix
                    loc_line = find(line-index+l-1==lines_1_copy(2,:));
                    loc_column = find(column+l-1==lines_1_copy(3,:));

                    % if this is a starting point of a line stored in the line
                    % matrix then check whether that line is longer than l_max
                    % - cnt + 1 than delete this line from the line
                    % matrix                
                    if any(ismember(loc_line,loc_column))                    
                        del_ind = loc_line(ismember(loc_line,loc_column));                                                
                        if lines_1_copy(1,del_ind) <= l_max 
                            lines_1_copy(:,del_ind) = [];
                        end                    
                    end
                
                elseif X(line+l-1,column+l-1+1) && last
                    
                    flag2 = true;
                    index2 = 1;
                    while flag2
                        if X(line+l-1,column+l-1+index2)
                            % delete that point
                            loc_line = find(line+l-1==lines_1_copy(2,:));
                            loc_column = find(column+l-1+index==lines_1_copy(3,:));
                            
                            if any(ismember(loc_line,loc_column))
                                del_ind = loc_line(ismember(loc_line,loc_column));
                                lines_1_copy(:,del_ind) = [];
                            end
                        else
                            flag2 = false;
                        end
                        index2 = index2 + 1;
                        if index2>N
                            break
                        end
                    end
                    
                else
                    % stop looking for neighbouring lines
                    flag = false;
                end
                index = index + 1;
            else
                
                % UPPER TRIANGLE
                
                % make sure not to exceed RP-boundaries
                if column+l-1-index == 0
                    break
                end

                if X(line+l-1,column-index+l-1)

                    % check whether this index is a listed index for starting
                    % points of line lengths in the line matrix
                    
                    loc_line = find(line+l-1==lines_1_copy(2,:));
                    loc_column = find(column-index+l-1==lines_1_copy(3,:));

                    % if this is a starting point of a line stored in the line
                    % matrix then check whether that line is longer than l_max
                    % - cnt + 1 than delete this line from the line
                    % matrix                
                    if any(ismember(loc_line,loc_column))  
       
                        del_ind = loc_line(ismember(loc_line,loc_column)); 
                        if lines_1_copy(1,del_ind) <= l_max                           
                            lines_1_copy(:,del_ind) = [];
                        end                    
                    end

                elseif X(line+l-1+1,column+l-1) && last
                    
                    flag2 = true;
                    index2 = 1;
                    while flag2
                        if X(line+l-1+index2,column+l-1)
                            % delete that point
                            loc_line = find(line+l-1+index2==lines_1_copy(2,:));
                            loc_column = find(column+l-1==lines_1_copy(3,:));
                            
                            if any(ismember(loc_line,loc_column))
                                del_ind = loc_line(ismember(loc_line,loc_column));
                                lines_1_copy(:,del_ind) = [];
                            end
                        else
                            flag2 = false;
                        end
                        index2 = index2 + 1;
                        if index2>N
                            break
                        end
                    end
                    
                else
                    % stop looking for neighbouring lines
                    flag = false;
                end
                index = index + 1;
            
            end
            
        end
               
        cnt = cnt + 1;
    end

end

% build RP based on the histogramm of the reduced lines

X_old = zeros(N,M);
% add LOI
X_old = double(X_old) + eye(size(X_old));
% fill up RP with lines stored in the new line matrix
for i = 1:size(lines_1_copy,2)
    l_max = lines_1_copy(1,i);
    line = lines_1_copy(2,i);
    column = lines_1_copy(3,i);
    for j = 1:l_max
        X_old(line+j-1,column+j-1) = 1;
    end
end

X_new = X_old;

% bind optional output
dl_new = lines_1_copy;


end
