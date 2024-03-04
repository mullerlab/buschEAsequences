% Code to process data from an example session
% Inputs: example session neural and behavioural data (found in ./data)
% Outputs: processed data file for one example session -- equivalent to one
% row in the data structures contained in ./processed_data/

%%%%%%%%%%%% Step 1: Process Spiking Data to get Sequences %%%%%%%%%%%%%%%%

clearvars; clc;
addpath( './helper_functions')

%%%%% CHOOSE PARAMETERS %%%%%

filetype = 'wm';
array = 2; % 0 for NSP0, 1 for NSP1, 2 for both
pval = 0.05; % for determining tuned cells
s = 0.1; % kernel size for SDF function

% set file name
if array == 0, arr = 'NSP0'; elseif array == 1, arr = 'NSP1'; else, arr = 'NSP0-1'; end
SaveFile = strcat( 'fr1_', filetype, '_', arr, '_spiketimes.mat');

%task info
timeOnset=1;
timeOffset=7000;
targetName = [ 6, 9, 8, ...
    7, 1, 3, ...
    12,11,10];
% targetNames = [1,3,6,7,8,9,10,11,12]; % original naming
% want to remap to the numbers 1-9 as follows:
% 1, 2, 3 = left column front, middle, back
% 4, 5, 6 = center column, front, middle, back
% 7, 8, 9 = right column, front, middle, back

% store
session_info = struct();
spiketimes = struct();

% process one example session
load('./data/example_session_neural.mat')

ff=1; % file number
session_info(ff).outcome = data.WM.trialOutcome; % save trial outcomes

% remap target names to target numbers 1:9
target_names = data.WM.cond;
targets = zeros( size(target_names));
for ii = 1: 9
    targets(target_names == targetName(ii)) = ii;
end
session_info(ff).targets = targets; % save trial conditions

tuned = 2; % both tuned and untuned cells
[raster_correct, raster_incorrect, neuronCond] = get_raster_info(data, array, tuned, pval, timeOnset, timeOffset);
session_info(ff).all_cells.neuronCond = neuronCond;

% correct trials
spiketimes.correct  = get_spiketimes(raster_correct, neuronCond, timeOnset, timeOffset, s);

% incorrect trials
spiketimes.incorrect = get_spiketimes(raster_incorrect, neuronCond, timeOnset, timeOffset, s);

% save
session_info(ff).all_cells.spiketimes = spiketimes;
%save(SaveFile, 'session_info')

%%

%%%%%%%%%%%%%%%%%%%% Step 2: Get Neural Boundaries %%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except ff session_info arr filetype array 
clc

db = 250; % bin width for setting neural boundaries
bins = 0:db:7000;

% load data
%load( strcat( './fr1_', filetype, '_', array, '_spiketimes.mat'))
SaveFile = strcat( 'fr1_', filetype, '_', arr, '_cell_selectivity.mat');

% thresholds for selectivity to scan over
height_threshold = 0.5;
width_threshold = 1000;

S = struct;
spiketimes = session_info(ff).all_cells.spiketimes.correct; % get spike times
Nc = size(spiketimes, 1); % number of cells
time_selective = zeros( Nc, 1);

for c = 1: Nc % for each cell
    
    tmp = spiketimes(c,:) ; % spike times of that cell across all correct trials of all 9 conditions
    [f,x] = ksdensity( tmp ); % get distribution of spike times
    m = max(f); % peak of distribution
    y = m*height_threshold; % height of dist to consider
    v = x(f>y); v = [v(1), v(end)]; % points where dist crosses that height
    width = abs( diff(v)); % observed width of dist at that threshold
    if width < width_threshold(ff) % if observed width is below width threshold, that cell is "time selective"
        time_selective(c) = 1;
    end
end

S(ff).times = spiketimes;
S(ff).time_selective = logical( time_selective );
S(ff).prct_selective = sum(time_selective)./Nc;

%save( SaveFile, 'S')

% Note: the neural boundaries defined in the paper are determined by
% repeating this process for each session, then pooling across all sessions
% of each subject to determine the peak times for that subject


%%
%%%%%%%%%%%%%%%%%%%% Step 3: Get Correlations %%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except ff session_info S arr filetype array 

% choose parameters 
boundary_type = 'neural'; % 'neural' or 'task'
epoch = 'delay'; % 'cue', 'delay', 'nav', 'delnav', 'all'
ndim = 1:3; % for projection
iterations = 100; % for shuffle
type = 'Spearman';  % for correlations

% task info
timeOnset=1;
delayOnset = 3000;
delayOffset = 5000;
timeOffset=7000;

load('./processed_data/fr1_wm_boundary_times.mat') % load neural boundaries (as computed in Step 2) 
load( './processed_data/fr1_wm_trajectories.mat') % load behavioural trajectories 
SaveFile = strcat( 'fr1_', filetype, '_', arr, '_neuralboundary_', epoch, '.mat');

goal_loc = session_data( 9 ).traj_dm; % load behavioural data corresponding to the example session of neural data 
subject = 'b';
spiketimes = session_info(ff).all_cells.spiketimes;
correct = struct; incorrect = struct;

% subset sequences by epoch, using desired boundaries
sequences_correct = get_epoch_sequences( spiketimes.correct, epoch, boundaries.(subject)(1), boundaries.(subject)(2));
sequences_incorrect = get_epoch_sequences( spiketimes.incorrect, epoch, boundaries.(subject)(1), boundaries.(subject)(2));
session_info(ff).all_cells.spiketimes.correct = sequences_correct;
session_info(ff).all_cells.spiketimes.incorrect = sequences_incorrect;

% get projection, correlations : correct trials
targets = session_info(ff).targets( session_info(ff).outcome );
[rho, D, V, prj] = get_projection_linear(sequences_correct, ndim);
[centroids, dm, cc, cc_shuff] = get_correlations( prj, targets, goal_loc, iterations, type);

% the "session_info" structure created here corresponds to one row of the
% "sesison_info" structures in the ./processed_data folder, specifically 
% row 15 of fr1_wm_NSPO-1_neuralboundary_delay.mat, corresponding to
% session date 2017-12-07










