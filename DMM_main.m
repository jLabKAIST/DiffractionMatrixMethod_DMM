%% Introduction
% Require Statistics & Machine Learning Toolbox
% Require RETICOLO RCWA by Jean-Paul Hugonin & Philippe Lalanne
% Hugonin, Jean Paul, and Philippe Lalanne. "Reticolo software for grating analysis." 
% arXiv preprint arXiv:2101.00901 (2021).
addpath('reticolo_allege');

clear
close all
fix(clock)
tic

%% Parameter Define
% All variables are in SI unit (m)
% 
param.wavelengths = 520*1e-9;     % Wavelength of light
param.period = 340*1e-9;          % Period of grating

% refractive indices at wavelength 520 nm (Set your refractive indicies)
n_top = 1.0;
n_org = 1.8;
n_Ag = 0.129807+3.09889i;
n_al= 0.83901+6.32423i;

Grating_thickness = 20*1e-9;      % Thickness of grating
Ag_thickness = 15*1e-9;           % Thickness of top Ag layer

param.nns = 20;                   % Fourier orders kept in RCWA (-nns:nns)
param.resolution = 70;            % # of points between 0th and 1st order in k-space

param.k0=2*pi/(param.wavelengths);
param.dx = param.wavelengths/param.period/n_org;  % spacing in u(k/nk0)-space region
param.uxsize = [-3.0001*param.resolution/param.dx:3*param.resolution/param.dx]*param.dx/param.resolution;
param.uysize = param.uxsize;


% define layers below/upper organic layer
% If periodic squared grating layer exists, set down.grating = 1
% If the layer is planar, down.grating = 0
% duty : the ratio of 2nd material (1: filled with 2nd material, 0: filled
% with 1st material)
% The grating is composed of only two diffeernt matrials

% starts from organic layer -> set down.n(1) to be n_org
% set thickness of 1st and last layer to be 0

down.grating = [0,1,0];
down.duty = [0,0.5,0];
down.n = [n_org,[n_org,n_al],n_al];
down.thickness = [0,Grating_thickness,0];

up.grating = [0,0,0];
up.duty = [0,0,0];
up.n = [n_org,n_Ag,n_top];
up.thickness = [0,Ag_thickness,0];

% [Logic Switch] Detect Grating Location
% Identify grating location and set flags 
% Error if both (up & down) have gratings
hasDownGrating = any(down.grating);  % True if any layer in 'down' has a grating
hasUpGrating   = any(up.grating);    % True if any layer in 'up' has a grating

% if RCWA simulation data exists, load file, or perform RCWA and save file
% Append structure info (Down/Up) to filename to prevent overwriting
if hasDownGrating
    struct_type = 'DownGrat';
else
    struct_type = 'UpGrat';
end

filename = ['RCWA\Orgair_' struct_type '_nn' num2str(param.nns) '_res' num2str(param.resolution) ...
            '_period_' num2str(param.period*1e9) '_G',num2str(Grating_thickness*10^9), 'nm.mat'];

if isfile(filename)

    load(filename)
    timeLoad = toc

else

    % RCWA/TMM process
    if hasDownGrating
        % Case 1: Down Grating
        % Down -> RCWA, Up -> TMM
        [rd,td,sd] = RCWA_DMM(param,down);
        [ru,tu,su] = TMM_DMM(param,up);
        
    elseif hasUpGrating
        % Case 2: Up Grating
        % Down -> TMM, Up -> RCWA
        [rd,td,sd] = TMM_DMM(param,down);
        [ru,tu,su] = RCWA_DMM(param,up);
    end

    % Record RCWA/TMM time
    timeRCWA = toc

    % Create folder if it doesn't exist
    if ~exist('RCWA', 'dir')
       mkdir('RCWA')
    end
    save(filename,'-v7.3')

end

%% Define DMM parameters
param.Org_thickness = 215*1e-9;         % thickness of Organic layer (from the top of down.grating to bottom of up.grating)
param.Dipole_pos = 55*1e-9;             % height of the dipole (from the top of down.grating to dipole)

param.Dox = 1;                          % 1: Calculate x-oriented dipole
param.Doy = 1;                          % 1: Calculate y-oriented dipole
param.Doz = 1;                          % 1: Calculate z-oriented dipole (if Dox,Doy,Doz ==1, calculate incoherently randomly oriented dipole)

param.sweep = 0;                            % 1: perform sweep
param.sweeprange.Org = [200:10:300]*1e-9;   % if sweep == 1, replace Org_thickness, give sweep range of thickness of org layer
param.sweeprange.Dip = [50:10:150]*1e-9;    % if sweep == 1, replace Dipole_pos, give sweep range of height of dipole

param.singledip = 1;                    % 0: incoherent planar emission, 1: single dipolar emission
param.x_dipole = 3/4*param.period;      % If singledip == 1, determine the position of dipole (from left side of slab)

param.uk = 1;                           % 1: Analyze power dissipation by its wavevectors
param.showstructure = 1;                % 1: show structure of corrugated OLED

%% DMM
if param.sweep ~= 1
    % Case: Single Run (Not a sweep)
    if hasDownGrating
        % If Grating is in Down structure -> Run DMM_process_down
        result = DMM_process_down(param,up,rd,ru,tu);
    else
        % If Grating is in Up structure -> Run DMM_process_up
        result = DMM_process_up(param,up,rd,ru,tu);
    end
else
    % Case: Sweep 
    if hasDownGrating
        % If Grating is in Down structure -> Run DMM_process_down
        result = DMM_sweep_down(param,up,rd,ru,tu);
    else
        % If Grating is in Up structure -> Run DMM_process_up
        result = DMM_sweep_up(param,up,rd,ru,tu);
    end
end

% Record calculation time
time = toc

%% Utility
% show structure real refractive index image
if param.showstructure == 1
    struc = showimage(param,up,down,su,sd);
end

% show far-field image
% if sweep == 1, show far-field image of 1st data
if param.Dox*param.Doy*param.Doz == 1

    figure
    if param.sweep ~= 1
        imagesc(result.Ifin(result.d:result.b,result.d:result.b))
    else
        imagesc(result.Ifin(1,1).data(result.d:result.b,result.d:result.b))
    end
    colorbar
    title(['Far-field Pattern (' struct_type ')']);

    % Analyze power dissipation by its wavevector
    if param.uk == 1
        [u,uK,K_prime,uK_u,uK_uK,uK_K_prime] = Get_uK(param,up,down,result);
    end

end