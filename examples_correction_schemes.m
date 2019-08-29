%% Applying correction schemes for border effects and tangential motion

% Here we exemplary show the funcionality of the different approaches in 
% this toolbox. 

% 1)    we generate the synthetic data
% 2)    we compute the RP's and the according line length histograms with
%       and without applying the correction schemes for border effects.
% 3)    we apply correction schemes for tangential motion on the RP's and
%       look at the resulting diagonal line length entropy with and without
%       applying the correction schemes for border effects.


clear, clc, close all

%% 0) Parameter setting and preallocation

% set length of the time series output
N = 2000; 

% select the threshold selection method and threshold:

% For the map data exampel:

% Here we apply a fixed recurrece threshold as 1/10 of the range of the
% underlying signal (this is done in section 1), where the map data will be
% generated)
thres_sel_1 = 'fix';

% For the flow data example:

% Here we apply a fixed recurrence threshold corresponding to 4% global 
% recurrence rate (we take the 4th percentile of all pairwise distances,
% cf. Kraemer et al., Chaos 28, 085720 (2018))
thres_sel_2 = 'var';
e2 = 0.04;

% select noise level (times the mean std of the multivariate signal)
sigma = 0; % 0% of the mean std of the signal, i.e. no noise contamination
% sigma = 0.1; % 10% of the mean std of the signal
% sigma = 0.3; % 30% of the mean std of the signal

% Choose minimal line length for computation of the entropy in Sect. 3.2)
% This is important, if you have chosen a non-zero noiselevel
l_min = 2;

% select the norm
norm = 'euc';
% norm = 'max';


% select type of border diagonal counting 

%'normal': lines, which start AND end at a RP boundary are denoted as
%          border diagonals
%'semi':   lines, which start OR end at a RP boundary are also denoted as
%          border diagonals

% type = 'normal';
type = 'semi';

% select the number of data points, which will be removed due to transient
% behavior
transients = 2500;

% preallocate storing matrices for
t = zeros(2,N); % time vectors for for regular and chaotic flow data

Y1 = cell(1,2); % phase space vectors for regular and chaotic map data
Y2 = cell(1,2); % phase space vectors for regular and chaotic flow data

RP1 = cell(6,2); % recurrence plots for regular and chaotic map data
RP2 = cell(6,2); % recurrence plots for regular and chaotic flow data

hl1 = cell(5,2); % all diagonal lines for regular and chaotic map data
hl2 = cell(5,2); % all diagonal lines for regular and chaotic flow data

pl1 = cell(5,2); % histogram of diagonal lines for regular and chaotic map data
pl2 = cell(5,2); % histogram of diagonal lines for regular and chaotic flow data

ent_Sum1 = cell(2,6); % entropies for regular and chaotic map data
ent_Sum2 = cell(2,6); % entropies for regular and chaotic flow data


%% 1) generate synthetic data

% set random number generator (for reproducibility)
rng(1)


%% 1.1) Map data 
% Logistic map in two parameter settings, leading to regular and chaotic 
% motion

% set parameters
a(1) = 3.5;
a(2) = 3.8;


x0 = rand; % set initial condition
xx(1,:) = ones(1,N+transients) * x0; % preallocate storing vector for regular motion
xx(2,:) = ones(1,N+transients) * x0; % preallocate storing vector for chaotic motion

% generate data
for i = 2:N+transients
    xx(1,i) = a(1) .* xx(1,i-1) .* ( 1 - xx(1,i-1));
    xx(2,i) = a(2) .* xx(2,i-1) .* ( 1 - xx(2,i-1));
end
x(1,:) = xx(1,transients+1:end);
x(2,:) = xx(2,transients+1:end);

% set recurrence thresholds as 10% of the range of the time series
e1(1) = range(x(1,:))/10;
e1(2) = range(x(2,:))/10;


% embed time series

% embedding dimension
m = 2;
% time delay
tau = 1;

% generate embedding vectors
Y1{1}= embed(x(1,:),m,tau);
Y1{2}= embed(x(2,:),m,tau);

% additive white noise
Y1{1} = Y1{1} + sigma*mean(std(Y1{1}))*randn(size(Y1{1}));
Y1{2} = Y1{2} + sigma*mean(std(Y1{2}))*randn(size(Y1{2}));

% Plot Phase Space

% Plot sampled Phase Space
figure('Units','normalized',...
    'Position',[.2 .2 .6 .5])
% markersize
sz = 20;

subplot(1,2,1)
scatter(Y1{1}(:,1),Y1{1}(:,2),sz,'*')
xlabel('x')
ylabel('y')
title('Logistic map (regular)')
set(gca,'Linewidth',2)
set(gca,'Fontsize',14)
grid on

subplot(1,2,2)
scatter(Y1{2}(:,1),Y1{2}(:,2),sz,'*')
xlabel('x')
ylabel('y')
title('Logistic map (chaotic)')
set(gca,'Linewidth',2)
set(gca,'Fontsize',14)
grid on


%% 1.2) Flow data 
% Roessler system in two parameter settings, leading to regular and chaotic 
% motion

% initial condition
x0 = [2*rand();2*rand();2*rand()];

% parameters
a = 0.15;
b(1) = 10; % parameter leading to regular motion
b(2) = 0.2; % parameter leading to chaotic motion
c = 10;

% integrate Roessler system with given parameters
[t(1,:),Y2{1}] = Rossler(a,b(1),c,x0,N,transients);
[t(2,:),Y2{2}] = Rossler(a,b(2),c,x0,N,transients);

Y2{1} = Y2{1}';
Y2{2} = Y2{2}';

Y2{1} = Y2{1} + sigma*mean(std(Y2{1}))*randn(size(Y2{1}));
Y2{2} = Y2{2} + sigma*mean(std(Y2{2}))*randn(size(Y2{2}));
    

% Plot Phase Space

% Plot sampled Phase Space
figure('Units','normalized',...
    'Position',[.2 .2 .9 .8])

subplot(1,2,1)
plot3(Y2{1}(:,1),Y2{1}(:,2),Y2{1}(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
title('numerically integrated Roessler Attractor (regular)')
set(gca,'Linewidth',2)
set(gca,'Fontsize',14)
grid on

subplot(1,2,2)
plot3(Y2{2}(:,1),Y2{2}(:,2),Y2{2}(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
title('numerically integrated Roessler Attractor (chaotic)')
set(gca,'Linewidth',2)
set(gca,'Fontsize',14)
grid on


%% 2) compute RPs and diagonal line length histograms of border correction schemes

% loop over the two different regimes
for i = 1:2
    
    % compute normal RPs for the Logistic map
    [RP1{1,i},~,~] = rp(Y1{i},e1(i),thres_sel_1,norm); 
    
    % compute normal RPs for the Roessler system
    [RP2{1,i},~,epsilon] = rp(Y2{i},e2,thres_sel_2,norm);     
    
    X1=RP1{1,i};
    X2=RP2{1,i};
    
    % get line distributions    
    [~, hl1{1,i}] = dl_conventional(X1); % black diagonal lines conventional
    [~, hl2{1,i}] = dl_conventional(X2); 
    
    [~, hl1{2,i}] = dl_dibo(X1,type); % dibo correction
    [~, hl2{2,i}] = dl_dibo(X2,type); 
    
    [~, hl1{3,i}] = dl_censi(X1,type); % Censi correction 
    [~, hl2{3,i}] = dl_censi(X2,type);
    
    [~, hl1{4,i}] = dl_kelo(X1,type); % kelo correction
    [~, hl2{4,i}] = dl_kelo(X2,type); 
        
    [~, hll] = dl_windowmasking(X1); % window masking
    hl1{5,i} = hll(1,:)';
    [~, hll2] = dl_windowmasking(X2); % window masking
    hl2{5,i} = hll2(1,:)';
     
    
    % loop over the different border effect correction approaches and make
    % a proper histogram
    for l = 1:5       
        % make proper histogram
        pl1{l,i} = hist(hl1{l,i},1:size(RP1{1,i},1)); 
        pl2{l,i} = hist(hl2{l,i},1:size(RP2{1,i},1)); 
    end
    
    
    % create string array for storing title names of the subplot, which we
    % will plot for all the correction schemes
    titlestr = {'conventional','dibo','Censi','kelo','window masking'};
    
    % Plot results for the Logistic map
    figure('Units','normalized',...
        'Position',[.1 .1 .9 .8])

    subplot(2,3,1)
    imagesc(RP1{1,i}(400:600,400:600)), colormap([1 1 1; 0 0 0]), axis xy square
    if i == 1
        title('RP of Logistic map (regular)')
    else
        title('RP of Logistic map (chaotic)')
    end
    
    for l = 1:5
        subplot(2,3,l+1)
        h = bar(1:size(RP1{1,i},1),pl1{l,i});
        col = get(h,'facecolor'); % color of the bar plot
        set(h,'edgecolor',col,'linewidth',1) % use the color also for surrounding line
        axis square
        xlabel('Length'),ylabel('Frequency')
        title(titlestr{l})
        set(gca,'LineWidth',1)
        set(gca,'FontSize',10) 
        yl = ylim; xl = xlim;
        text(xl(1) - .3*range(xl),yl(2) + .18*range(yl),char(64+i),'FontWeight','bold','FontSize',14)      
        
        % set x- and y- limits in the plots
        if i == 1
            ylim([0 3])
            xlim([0 N+10])
        else
            xlim([0 N/40])            
        end
    end
    
    % Plot results for the Roessler system
    figure('Units','normalized',...
        'Position',[.1 .1 .9 .8])

    subplot(2,3,1)
    imagesc(t(i,:),t(i,:),RP2{1,i}), colormap([1 1 1; 0 0 0]), axis xy square
    if i == 1
        title('RP of Roessler system (regular)')
    else
        title('RP of Roessler system (chaotic)')
    end
    
    for l = 1:5
        subplot(2,3,l+1)
        h = bar(1:size(RP2{1,i},1),pl2{l,i});
        col = get(h,'facecolor'); % color of the bar plot
        set(h,'edgecolor',col,'linewidth',1) % use the color also for surrounding line
        axis square
        xlabel('Length'),ylabel('Frequency')
        title(titlestr{l})
        set(gca,'LineWidth',1)
        set(gca,'FontSize',10) 
        yl = ylim; xl = xlim;
        text(xl(1) - .3*range(xl),yl(2) + .18*range(yl),char(64+i),'FontWeight','bold','FontSize',14)      
        
        % set x- and y- limits in the plots
        if i == 1
            ylim([0 3])
            xlim([0 N+10])
        else
            xlim([0 N/10])            
        end
    
    end
    
    
end


%% 3) tangential motion correction schemes

% set parameters for correction schemes

% perpendicular RP: angle threshold
w = 0.258;

% isodirectional RP: translation parameter T
T = 8;
% isodirectional RP: choice of second recurrence threshold. Give this as a
% percentage of the first recurrence threshold, i.e. if set to 0.8, the
% second threshold will be e2 = 0.8 * e1
iso_threshold_factor = 0.8;

% TRP threshold
TRP_threshold = 8;

% LM2P RP threshold
LM2P_threshold = 8;

%% 3.1) Compute and display RPs from the tangential motion correction schemes

% loop over the two different regimes
for i = 1:2
    
    % apply correction schemes for tangential motion
   
    % Logistic map data
    RP1{2,i} = rp_perp(Y1{i},e1(i),thres_sel_1,w,norm);
    RP1{3,i} = rp_iso(Y1{i},e1(i),e1(i)*iso_threshold_factor,thres_sel_1,T,norm);
    RP1{4,i} = rp_TRP(Y1{i},e1(i),thres_sel_1,TRP_threshold,norm); 
    RP1{5,i} = rp_LM2P(Y1{i},e1(i),thres_sel_1,LM2P_threshold,norm);
    RP1{6,i} = rp_diagonal(RP1_normal{i});
    
    % Roessler system data
    RP2{2,i} = rp_perp(Y2{i},e2,thres_sel_2,w,norm);
    RP2{3,i} = rp_iso(Y2{i},e2,epsilon*iso_threshold_factor,thres_sel_2,T,norm);
    RP2{4,i} = rp_TRP(Y2{i},e2,thres_sel_2,TRP_threshold,norm); 
    RP2{5,i} = rp_LM2P(Y2{i},e2,thres_sel_2,LM2P_threshold,norm);
    RP2{6,i} = rp_diagonal(RP2_normal{i});
   
   
    % plot results 
    
    methods = {'normal','perpendicular','isodirectional','True','LM2P',...
        'diagonal'};
    
    % Logistic map
    
    figure('Units','normalized',...
         'Position',[.1 .1 .9 .8])
    % loop over the correction schemes
    for j = 1:6
        
        subplot(2,3,j)
        imagesc(RP1{j,i}), colormap([1 1 1; 0 0 0]), axis xy square
        if i == 1
            titlestring = strcat(methods{j},'{ }','RP (Logistic map regular)');
            title(titlestring)
        else
            titlestring = strcat(methods{j},'{ }','RP (Logistic map chaotic)');
            title(titlestring)
        end
        xlabel('datapoint index')
        ylabel('datapoint index')
        set(gca,'LineWidth',2)
        set(gca,'FontSize',14)     
    end
    
    % Roessler system
    
    figure('Units','normalized',...
         'Position',[.1 .1 .9 .8])
    % loop over the correction schemes
    for j = 1:6
        
        subplot(2,3,j)
        imagesc(RP2{j,i}), colormap([1 1 1; 0 0 0]), axis xy square
        if i == 1
            titlestring = strcat(methods{j},'{ }','RP (Roessler system regular)');
            title(titlestring)
        else
            titlestring = strcat(methods{j},'{ }','RP (Roessler system chaotic)');
            title(titlestring)
        end
        xlabel('datapoint index')
        ylabel('datapoint index')
        set(gca,'LineWidth',2)
        set(gca,'FontSize',14)     
    end
      
end
    

%% 3.2) Compute and display entropy from the tangential motion correction schemes
    
% loop over the two different regimes
for i = 1:2
    
    % loop over the different tangential motion correction schemes
    
    for j = 1:6
        
        X1 = RP1{j,i};
        X2 = RP2{j,i};
        
        % get line distributions    
        [~, hl1{1,i}] = dl_conventional(X1); % black diagonal lines conventional
        [~, hl2{1,i}] = dl_conventional(X2); 

        [~, hl1{2,i}] = dl_dibo(X1,type); % dibo correction
        [~, hl2{2,i}] = dl_dibo(X2,type); 

        [~, hl1{3,i}] = dl_censi(X1,type); % Censi correction 
        [~, hl2{3,i}] = dl_censi(X2,type);

        [~, hl1{4,i}] = dl_kelo(X1,type); % kelo correction
        [~, hl2{4,i}] = dl_kelo(X2,type); 

        [~, hll] = dl_windowmasking(X1); % window masking
        hl1{5,i} = hll(1,:)';
        [~, hll2] = dl_windowmasking(X2); % window masking
        hl2{5,i} = hll2(1,:)';
        
        % loop over the different border effect correction approaches and make
        % a proper histogram
        for l = 1:5       
            % make proper histogram
            pl1{l,i} = hist(hl1{l,i},1:size(RP1{1,i},1)); 
            pl2{l,i} = hist(hl2{l,i},1:size(RP2{1,i},1)); 
            
            % adjust histograms to chosen minimal line length
            pl1{l,i} = pl1{l,i}(lmin:end);
            pl2{l,i} = pl2{l,i}(lmin:end);
            
            % compute diagonal line length entropy

            % get probability distribution from histogram
            l_prob1 = pl1{l,i}/sum(pl1{l,i});
            l_prob2 = pl2{l,i}/sum(pl2{l,i});

            ent_Sum1{i,j}(l) = l_prob1 .* log(l_prob1)/log(length(pl1{l,i}));
            ent_Sum2{i,j}(l) = l_prob1 .* log(l_prob1)/log(length(pl2{l,i}));
        end
        
    end  
    
    % Plot entropies
    
    figure('Units','normalized',...
         'Position',[.1 .1 .9 .8])
    % loop over the correction schemes
    for j = 1:6
        
        subplot(2,3,j)       
        bar(1:5,ent_Sum1{i,j}(:))
        set(gca,'xticklabel',{titlestr{1},titlestr{2},titlestr{3},...
            titlestr{4},titlestr{5}})
        set(gca,'LineWidth',2)
        set(gca,'FontSize',14) 
        ylabel('normalized entropy')
        
        if i == 1
            titlestr = strcat('Entropy of ',methods{j},'{ }','RP (Logistic map regular)');
            title(titlestr)
        else
            titlestr = strcat('Entropy of ',methods{j},'{ }','RP (Logistic map chaotic)');
            title(titlestr)
        end

    end
    
    figure('Units','normalized',...
         'Position',[.1 .1 .9 .8])
    % loop over the correction schemes
    for j = 1:6
        
        subplot(2,3,j)       
        bar(1:5,ent_Sum2{i,j}(:))
        set(gca,'xticklabel',{titlestr{1},titlestr{2},titlestr{3},...
            titlestr{4},titlestr{5}})
        set(gca,'LineWidth',2)
        set(gca,'FontSize',14) 
        ylabel('normalized entropy')
        
        if i == 1
            titlestr = strcat('Entropy of ',methods{j},'{ }','RP (Roessler system regular)');
            title(titlestr)
        else
            titlestr = strcat('Entropy of ',methods{j},'{ }','RP (Roessler system chaotic)');
            title(titlestr)
        end

    end
    
    
end
 
    

%% Helper functions

function [t2,x2] = Rossler(varargin)
% This function numerically integrates the Roessler-System with input
% paramters a, b and c starting with input inital condition 
% x0 = [x(t0); y(t0); z(t0)]. 
%
% [t,x] = Rossler(a,b,c,x0,N,transients)
%
% t is the time vector containing just the times, where the samping took
% place. The corresponding values to these t2-times are stored in x.
%
% Here the sampling rate is 0.2 and the output is a vector-series of 
% length N (input parameter) with 'transients' removed samples.


% bind input
a = varargin{1};
b = varargin{2};
c = varargin{3};
x0= varargin{4};
N = varargin{5};
transients = varargin{6};


% Define your sampling time
time_sample = 0.2;

time_interval = (N+transients)* time_sample;     % k+1 is the number of columns in the solution matrix M

% This results in a corresponding time vector
t2 = 0:time_sample:time_interval;

 
f = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)];  % Rossler system
SOL = ode45(f,[t2(1) t2(end)],x0);
    
x2 = deval(SOL,t2);

% Store sampled time series    
x2 = x2(:,transients+1:end-1);
t2 = t2(1:end-transients-1);

end

    
