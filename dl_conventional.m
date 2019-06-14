function [a_out, b_out]=dl_conventional(x)
% DL_CONVENTIONAL   Mean of the diagonal line lengths and their distribution.
%    A=dl_conventional(X) computes the mean of the length of the diagonal 
%    line structures in a recurrence plot.
%
%    [A B]=dl_conventional(X) computes the mean A and the lengths of the
%    found diagonal lines, stored in B. In order to get the 
%    histogramme of the line lengths, simply call 
%    HIST(B,[1 MAX(B)]).
%
%    Examples (CRP toolbox needs to be installed):
%       x = sin(linspace(0,5*2*pi,1050));
%       xe = embed(x,2,50);
%       r = rp(xe,.2);
%       [l l_dist] = dl_conventional(r);
%       subplot(1,2,1)
%       imagesc(r), colormap([1 1 1;0 0 0]), axis xy square
%       title('underlying RP')
%       subplot(1,2,2)
%       histogram(l_dist,1000)
%       xlabel('diagonal line length')
%       ylabel('counts')
%       title('diagonal line length histogram - conventional counting')
%
%    See also CRQA, TT.
%
% Copyright (c) 2008-
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

warning off
if any(x(:))

  if min(size(x))>100000       % this should speed up the routine; the value
                             % depends on the available memory
    x2=uint8(x);
    N=size(x2);
    x3=zeros(2*N(2)+N(1),N(2));
    x3(N(2)+1:N(2)+N(1),1:N(2))=x2;
    N3=size(x3);
    
    i2=repmat(((1:1+N(2))+N(1)+N(2))',1,N(2));
    i4=i2+repmat((2*N(2)+N(1)+1)*[0:N(2)-1],size(i2,1),1);
    i4(:,end)=[];
    i4=reshape(i4,size(i4,1)*size(i4,2),1);
    x3(i4)=[];
    x3(end)=[];
    x2=(reshape(x3,N(1)+N(2),N(2)))';
  
    x2(end+1,:)=0;
    x=reshape(x2,size(x2,1)*size(x2,2),1);
    x2=x(2:end);x(end)=[];
    z0=find(x==0&x2==1);
    z1=find(x2==0&x==1);
  
  else
  
    N=size(x);
  %   x3=zeros(2*N(2)+N(1),N(2));
  %   x3(N(2)+1:N(2)+N(1),1:N(2))=x;
  %   N3=size(x3);
  %   
  %   i2=repmat(((1:1+N(2))+N(1)+N(2))',1,N(2));
  %   i4=i2+repmat((2*N(2)+N(1)+1)*[0:N(2)-1],size(i2,1),1);
  %   i4(:,end)=[];
  %   i4=reshape(i4,size(i4,1)*size(i4,2),1);
  %   x3(i4)=[];
  %   x3(end)=[];
  %   x=(reshape(x3,N(1)+N(2),N(2)))';
  %  
  %   x(end+1,:)=0;
    
  %  for i1=-ceil(N(2)/2):ceil(N(2)/2); temp=diag(x,i1); X(1:length(temp),1+i1+ceil(N(2)/2))=temp;
  %  end, x=double(X);
    x1=spdiags(double(x));
    z=reshape(x1,size(x1,1)*size(x1,2),1);
    z2(2:length(z)+1)=z;z2(1)=0;z2(end+1)=0;
    z=diff(z2);
    z0=find(z==1);
    z1=find(z==-1);
  
  end
  if length(z0)>length(z1), z0(end)=[]; end
  if length(z1)>length(z0), z1(end)=[]; end
  
  if isempty(z0), z0=0; end
  if isempty(z1), z1=0; end
  
  if z0(1)>z1(1)
    z0(2:end+1)=z0(1:end);z0(1)=0; 
    if length(z0)>length(z1) 
       z0(end)=[];
    end
  end

  l=sort(z1-z0); %l(end)=[];
  l1=l(find(l-1));
  
  if nargout==2
     b_out=zeros(length(l),1);
     b_out=l';
  end
  
  if nargout>0
     a_out=mean(l1);
  else
     mean(l1)
  end
  
else

  if nargout==2
     b_out=NaN;
  end

  if nargout>0
     a_out=NaN;
  else
     NaN
  end

end

warning on
