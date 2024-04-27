%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             simulate 3D UTE data for multiple TEs   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Original Code written by Peter J. Niedbalski and edited by Abdullah to include noise and
% multiple TEs and calcualte T2* for undersampling comparison
% Ian R. Stecker made the rest of the changes 

% For any question, please contact:
%   Author: Abdullah Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
 

clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Imaging Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Things needed to set-up 
% IS NOTE: Basically I want to set up a range of values from 0.3 to 0.5 ms
% lets say I want to do this randomly and apply onlmatlaby to the lungs. I've
% added code to do that but for all simulated tissues

T2s = 0.0004; %In Seconds - Mouse lung in vivo T2* is ~0.4 ms, % Not using this anymore!!!

ImSize = 64; % 128 is what we usually run at for mouse imaging 
GenSize = ImSize*2; % THIS is to help reduce radial artifacts

SampPct_all = 100; % Desired Sampling Percentage will be used to calculate number of projections

% Gaussian noise variable is included for any future simulation works
noisesd = 50; % SD of the added noise-- 

Save_SimulationResult = 'no'; % if you choose "yes", a window will pop up to choose where you want to save the simulation workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things you might want to leave alone as they are (defult values)
TR = 0.007; %In Seconds
TE= [0.080; 0.250; 0.500; 1.250; 2.500];  %In mSeconds
nTE=length(TE);
Dwell = 1/277777; %In Seconds (We usually run mouse scans at BW = 277777, hence this number for the dwell
RampComp = 'Yes'; %Being Rigorous about simulation, and using a "real" gradient shape with a ramp up
      %Do we want to compensate for points acquired on the ramp
      %(Usually we do this in practice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arbitrarily create some general gradient trapezoid ~1/3 ramp, the rest 
GradShape = [linspace(0,1,33) ones(1,67)];

%Get Trajectory Shape from Arbitrary Grad Shape by integrating
TrajShape = cumtrapz(GradShape);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SampPct = SampPct_all;
TE = TE/1000;  %In Seconds  

% Based on Imaging Variables, calculate projections 
ptsPerPro = ImSize/2;  %Theoretical Number of points along a projection is 1/2 image size
NPro_FS = ceil(4*pi*(ptsPerPro^2)); %Number of Projections for Fully sampled image - Nyquist Criterion
NPro_Act = ceil(SampPct/100*NPro_FS); %Number of Projections for Scan, based on SampPct;

%In order to have the same resolution as initially, we want to compensate
%for points acquired on the ramp
if(strcmp(RampComp,'Yes'))
    nPts = length(find(TrajShape<=ptsPerPro));
else
    nPts = ptsPerPro;
end
TrajShapeFin = TrajShape(1:nPts); %Final Trajectory Shape to use for Sampling

%Golden Means EigenValues (From Jinbang Guo) and gradient multipliers
phi1 = 0.46557123;
phi2 = 0.6823278;
gs = 1;
gr = 1;
gp = 1;

%Preallocate memory
r = zeros(1,NPro_Act);
p = zeros(1,NPro_Act);
s = zeros(1,NPro_Act);
alpha = zeros(1,NPro_Act);
kz = zeros(1,NPro_Act);

%Rotation code from Jinbang's UTE sequence (Calculate Read, Phase, and
%Slice Components of each trajectory rotation using Golden Means calculation)
for i = 0:(NPro_Act-1)
    kz(i+1) = (i*phi1-floor(i*phi1))*2-1;
    ts = kz(i+1)*gs;
    alpha(i+1) = (i*phi2-floor(i*phi2))*2*pi;
    tr = sqrt(1-kz(i+1)*kz(i+1))*cos(alpha(i+1))*gr;
    tp = sqrt(1-kz(i+1)*kz(i+1))*sin(alpha(i+1))*gp;
    r(i+1) = tr;
    p(i+1) = tp;
    s(i+1) = ts;
end
zang = asin(kz);

%Preallocate memory
trajx = zeros(length(TrajShapeFin),NPro_Act);
trajy = trajx;
trajz = trajx;

%Calculate all x, y, and z trajectories
for i = 1:NPro_Act
    trajx(:,i) = r(i)*TrajShapeFin';
    trajy(:,i) = p(i)*TrajShapeFin';
    trajz(:,i) = s(i)*TrajShapeFin';
end

%Combine everything into a single trajectory matrix
traj = cat(3,trajx,trajy,trajz);

%Permute trajectory matrix to a format more pleasing to my eye
traj = permute(traj,[3,1,2]);

% % End of Trajectory Calculations % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phantom and noise Creation and Decay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create 3-Dimensional Shepp-Logan Phantom using phantom3d function from
%internet: We could make some other 3D phantom, but this is easy.

% 07/08/2020 IRS Note: I created another phantom within phantom3d.m, I
% changed the phantom size and added additional ellipsoids to make phantom
% appear more lung-like to illustrate pulmonary vascualture voxels. The
% phantom3d.m function now defaults to the phantom model:
% irs_undersampling. Please note, that the modified shepp-logan phantom
% still exists within this function, an additional input arg is required
% if you wish to switch to this model instead (e.g., phtm =
% phantom3d(128,'modified shepp-logan');) [this is from the original
% writer's code]

MyPhantom = flip(phantom3d('irs-fakemouse', GenSize));
disp('Creating Phantom Completed')

%We could do a lot of really rigorous simulation, but for this version,
%just do T2* simulation --> i.e start with main phantom and each successive
%point on a ray should be decayed according to T2* - First phantom should
%also be decayed according to T2star

%Leave out Flip angle and T1 considerations (For Proton, I think it's
%largely irrelevant... not too hard to add if that is not the case.

%Decay Phantom Images according to T2star, and FFT into k-space for each
%decay step - This takes a bit of time for the FFT

% Implementing Varying T2* to the Lung, Heart, and Surrounding Tissues
% First, let's just change the value within the lungs, don't like
% it when signal is less than 0

% 07/08/2020 IRS Note: I have change the signal intensity of the phantom to
% more closely match the signal seen commonly for the in-vivo images
MyPhantom(MyPhantom < 0) = 500000; % original value: 0; new set value is close to in-vivo images
MyPhantom(MyPhantom > 0.1 & MyPhantom < 0.2) = 1000000; % original value: 0.2;
MyPhantom(MyPhantom > 0.2 & MyPhantom < 0.3) = 900000; % original value: 0.3;
MyPhantom(MyPhantom > 0.4 & MyPhantom < 1)   = 1050000; % original value: 0.5;
MyPhantom(MyPhantom == 1) = 1000000; % original value: 1; No outer rings present within in-vivo images, setting value to mouse tissue instead

% Now let's make some binary masks
% 07/08/2020 IRS Note: I am adding an additional mask for the simulated
% lung vasculature
ParenMask  = double(MyPhantom == 500000);
TissueMask = double(MyPhantom == 1000000);
HeartMask  = double(MyPhantom == 900000);
VascMask   = double(MyPhantom == 1050000);
NoiseMask  = double(MyPhantom == 0);

% Now, let's add some varying T2* based on size of phantom
% The size will change if using larger matrix sizes

% 07/08/2020 IRS Note: Here, I am only assigning a single T2* for each
% masked section. 

% Assigning T2* of 0.4ms to the ParenMask (i.e., lung parenchyma voxels)
T2Paren = ones(size(ParenMask(ParenMask == 1)));
T2Paren(1:round(size(T2Paren,1)/4,0)) = 0.0004;
T2Paren(round(size(T2Paren,1)/4,0)+1:round(size(T2Paren,1)*2/4,0)) = 0.0004;
T2Paren(round(size(T2Paren,1)*2/4,0)+1:round(size(T2Paren,1)*3/4,0)) = 0.0004;
T2Paren(round(size(T2Paren,1)*3/4,0)+1:end) = 0.0004;

% I actually don't know what the T2 of the heart or surrounding tissue
% would be at 7T. Defaulting to parenchyma T2* 
T2Heart = ones(size(HeartMask(HeartMask == 1)));
T2Heart(1:round(size(T2Heart,1)/4,0)) = 0.0004;
T2Heart(round(size(T2Heart,1)/4,0)+1:round(size(T2Heart,1)*2/4,0)) = 0.0004;
T2Heart(round(size(T2Heart,1)*2/4,0)+1:round(size(T2Heart,1)*3/4,0)) = 0.0004;
T2Heart(round(size(T2Heart,1)*3/4,0)+1:end) = 0.0004;

T2Tissue = ones(size(TissueMask(TissueMask == 1)));
T2Tissue(1:round(size(T2Tissue,1)/4,0)) = 0.0004;
T2Tissue(round(size(T2Tissue,1)/4,0)+1:round(size(T2Tissue,1)*2/4,0)) = 0.0004;
T2Tissue(round(size(T2Tissue,1)*2/4,0)+1:round(size(T2Tissue,1)*3/4,0)) = 0.0004;
T2Tissue(round(size(T2Tissue,1)*3/4,0)+1:end) = 0.0004;

% Here I do know what the T2* of the vasculature is: ~3.16ms
T2Vasc = ones(size(VascMask(VascMask == 1)));
T2Vasc(1:round(size(T2Vasc,1)/4,0)) = 0.00316;
T2Vasc(round(size(T2Vasc,1)/4,0)+1:round(size(T2Vasc,1)*2/4,0)) = 0.00316;
T2Vasc(round(size(T2Vasc,1)*2/4,0)+1:round(size(T2Vasc,1)*3/4,0)) = 0.00316;
T2Vasc(round(size(T2Vasc,1)*3/4,0)+1:end) = 0.00316;

% Preallocate memory for T2* variables of different tissues
addT2Lung = zeros(GenSize, GenSize, GenSize);
addT2Heart = zeros(GenSize, GenSize, GenSize);
addT2Tissue = zeros(GenSize, GenSize, GenSize);
addT2Vasc = zeros(GenSize, GenSize, GenSize);

% Find the locations where the binary masks == 1 for each tissue
LungLoc = find(ParenMask==1);
HeartLoc = find(HeartMask==1);
TissueLoc = find(TissueMask==1);
VascLoc = find(VascMask==1);

% Lastly, add the varying T2* for each tissue in the binary locations
% Changing 0's to 1's to not delete any information about the phantom
addT2Lung(LungLoc) = T2Paren; %addT2Lung(addT2Lung == 0) = Inf;
addT2Heart(HeartLoc) = T2Heart; %addT2Heart(addT2Heart == 0) = Inf;
addT2Tissue(TissueLoc) = T2Tissue; %addT2Tissue(addT2Tissue == 0) = Inf;
addT2Vasc(VascLoc) = T2Vasc; %addT2OutEdge(addT2OutEdge == 0) = Inf;
addT2Combined = addT2Lung + addT2Heart + addT2Tissue + addT2Vasc;

% Clear no longer needed variables
clear addT2Lung LungLoc addT2Heart HeartLoc addT2Tissue TissueLoc addT2Vasc
clear VascLoc LungMask HeartMask TissueMask VascMask NoiseMask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare for K-space Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For sampling, trajectories need to be centered at ImSize/2 and go between
% 0 and ImSize - Do that here
%         SamplingTraj = traj + ImSize / 2;
SamplingTraj = traj + GenSize / 2;

%For Recon, trajectories need to be centered at 0 and go between -.5 and .5
radius = sqrt(traj(1,:,:).^2 + traj(2,:,:).^2 + traj(3,:,:).^2);
ReconTraj = traj / max(radius(:))/2;

% Round sampling trajectories to the nearest grid point so that we don't
% have to interpolate
SamplingTraj_round = round(SamplingTraj);

%Since we're rounding to the nearest grid point, we need to edit 
ReconTraj_round = SamplingTraj_round - GenSize / 2;
rad = sqrt(ReconTraj_round(1,:,:).^2 + ReconTraj_round(2,:,:).^2 + ReconTraj_round(3,:,:).^2);

ReconTraj2 = ReconTraj_round / max(rad(:)) / 2;

%Pre-Allocate Memory
KSpaceData2 = zeros(nPts, NPro_Act, nTE);

%% Sanity Check 1 - Look at first 100 generated trajectories
% figure('Name','Sanity Check 1 - Generated Trajectories')
% hold on
% for i = 1:250
% %     plot3(ReconTraj_round(1,:,i),ReconTraj_round(2,:,i),ReconTraj_round(3,:,i));
%     plot3(traj(1,:,i),traj(2,:,i),traj(3,:,i));
% end
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clearing Out Variables no Longer Needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear addT2Heart addT2Lung addT2OutEdge addT2Tissue
clear alpha GenSize gp gr gs GradShape HeartLoc HeartMask i kz LungLoc
clear LungMask NoiseMask OutEdgeLoc OutEdgeMask phi1 phi2 r radius RampComp
clear T2Heart T2Lung T2Vasc T2s T2Tissue TissueLoc TissueMask tp tr TrajShape
clear TrajShapeFin trajx trajy trajz ts traj p s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample K-Space - This is probably one of the longest steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method using rounding of trajectory points
% MyPhantom_T2s = zeros(GenSize,GenSize,GenSize,nPts,nTE);
% MyPhantom_KSpace = zeros(GenSize,GenSize,GenSize,nPts,nTE);
tic
for k = 1:nTE
    parfor i = 1:nPts
        MyPhantom_T2s = MyPhantom.*exp(-(TE(k)+(i-1)*Dwell)./addT2Combined);
        MyPhantom_KSpace = fftshift(ifftn(ifftshift(MyPhantom_T2s)));
        for j = 1:NPro_Act
            KSpaceData2(i,j,k) = MyPhantom_KSpace(SamplingTraj_round(1,i,j),SamplingTraj_round(2,i,j),SamplingTraj_round(3,i,j));
        end
    end
end
Round_Dur = toc/60
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear no longer needed variables
clear MyPhantom_T2s MyPhantom_KSpace Noise_phantom_KSpace MyPhantom_KSpaceRecon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add noise
N_real= normrnd(0,noisesd,size(KSpaceData2)); %real
N_imag= normrnd(0,noisesd,size(KSpaceData2)); %imaginary  
Noise_phantom = complex(N_real,N_imag); 
KSpaceData2 = KSpaceData2 + Noise_phantom;

% Make sure we don't have NaNs in the data:
KSpaceData2(isnan(KSpaceData2)) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare Data for Recon - Data and Trajectories both need to be in columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for k = 1:nTE
    FID_Col2(:,k) = reshape(KSpaceData2(:,:,k),1,[])';
end
Recon_Traj2x = reshape(ReconTraj2(1,:,:),1,[])';
Recon_Traj2y = reshape(ReconTraj2(2,:,:),1,[])';
Recon_Traj2z = reshape(ReconTraj2(3,:,:),1,[])';

Recon_Traj2 = [Recon_Traj2x Recon_Traj2y Recon_Traj2z];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recon Original Image - Uses Scott Robertson's Recon Code placed within a function for ease of use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear no longer needed variables
clear KSpaceDataAddNoise KSpaceOriginal
Recon2 = zeros(ImSize, ImSize, ImSize, nTE);
for k = 1:nTE
    Recon2(:,:,:,k) = ScottRecon3D_deap(ImSize,FID_Col2(:,k),Recon_Traj2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sanity Check 5 - Look at Reconstructed Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotate and flip reconstruction to be in correct orientation
FinalRecon = Recon2;

clear Recon2 Recon_Traj2 Recon_Traj2x Recon_Traj2y Recon_Traj2z ReconTraj ReconTraj2
clear ReconTraj_round SamplingTraj SamplingTraj_round FID_Col2 KSpaceData2
clear rad TR zang

% Final_Recon = rot90(abs(FinalRecon),1);
% Final_Recon = flip(FinalRecon,1);

% % First we need the Final Recon Phantom to have prescribed binary masks
PrePhantom = flip(phantom3d(ImSize));
PrePhantom(PrePhantom < 0) = 0.1;
PrePhantom(PrePhantom > 0.1 & PrePhantom < 0.2) = 0.2;
PrePhantom(PrePhantom > 0.2 & PrePhantom < 0.3) = 0.3;
PrePhantom(PrePhantom > 0.4 & PrePhantom < 1)   = 0.5;
PrePhantom(PrePhantom == 1)                     = 0.2;

clear Final_Recon Final_Recon_with_bgr_all
clear Noise_phantom_KSpace MyPhantom_KSpace KSpaceData KSpaceNoise Recon_Traj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the final binary masks of prescribed matrix size
FinParenMask = double(PrePhantom == 0.1);
FinHeartMask = double(PrePhantom == 0.3);
FinTissueMask = double(PrePhantom == 0.2);
FinVascMask = double(PrePhantom == 0.5);

FinLungMask = FinParenMask + FinVascMask;
[x, y, z] = ndgrid(-1:1);
hood = sqrt(x.^2 + y.^2 + z.^2) <= 1.0;
FinLungMask = imerode(FinLungMask, hood);

% This ensures erosion along outer edge of mask
FinParenMask = FinParenMask .* FinLungMask;

FinNoiseMask = double(PrePhantom == 0);
[xx,yy,zz] = ndgrid(-5:5);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 5.0;
FinNoiseMask = imerode(FinNoiseMask, nhood);

signal_n_avg = zeros(1,nTE);
mean_noise_n = zeros(1,nTE);
std_noise_n = zeros(1,nTE);
SNR_vec_n = zeros(1,nTE);

for n = 1:nTE
    %Calculating the siganl vector
    signal_vec = abs(FinalRecon(:,:,:,n));
    signal_vec(FinParenMask == 0) = [];
    signal_n_avg(n) = mean(signal_vec);

    %Calculating the noise vector
    noise_vec = abs(FinalRecon(:,:,:,n));
    noise_vec(FinNoiseMask == 0) = [];
    mean_noise_n(n) = mean(noise_vec);               % the MEAN of the noise    
    std_noise_n(n) = std(noise_vec);                 % the standard deviation of the noise
    SNR_vec_n(n) = signal_n_avg(n) / std_noise_n(n); % signal to noise ratio
end

SNR_vec_n_all = SNR_vec_n;
disp(num2str(SNR_vec_n))
disp('Calcuating SNR Completed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating a T2* Mask of Image (NonLinear Fit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making an empty array for T2
T2array_nonlinear = zeros(ImSize^3, 1);
SoArray_nonlinear = zeros(ImSize^3, 1);
arrayLength = ImSize^3;

% Using parallel looping for nonlinear modeling
disp('Calculating T2*...Please wait...')
modeloptions = fitoptions('Method', 'NonlinearLeastSquares',...
           'Lower', [0, 0, 0],...
           'Upper', [Inf, Inf, Inf],...
           'StartPoint', [1e6 0.35 1e3]);
% Echo Times used in Experiment (in ms)
TE = [0.080; 0.250; 0.500; 1.250; 2.500]; %%3/7/19 - IS avoiding UTE use

% Equation to Fit Data of Lung Images
t2fit = fittype('a * exp(-x / b) + c', 'coeff', {'a', 'b', 'c'}, 'options', modeloptions);

Im_UTE_TE1 = abs(FinalRecon(:,:,:,1));
Im_UTE_TE2 = abs(FinalRecon(:,:,:,2));
Im_UTE_TE3 = abs(FinalRecon(:,:,:,3));
Im_UTE_TE4 = abs(FinalRecon(:,:,:,4));
Im_UTE_TE5 = abs(FinalRecon(:,:,:,5));

parfor index = 1:arrayLength
    if FinParenMask(index) == 1
        y = [Im_UTE_TE1(index); Im_UTE_TE2(index); Im_UTE_TE3(index); Im_UTE_TE4(index); Im_UTE_TE5(index)];
        f = fit(TE, y, t2fit);
        coeff_values = coeffvalues(f);
        SoArray_nonlinear(index) = coeff_values(1);
        T2array_nonlinear(index) = coeff_values(2);
    end
end
toc

T2Image = reshape(T2array_nonlinear, [ImSize ImSize ImSize]);
SoImage = reshape(SoArray_nonlinear, [ImSize ImSize ImSize]);
Mean_T2 = T2Image;
Mean_T2(Mean_T2<0) = 0;
Mean_T2(Mean_T2>100)=0;

Mean_T2(FinParenMask == 0) = [];
Mean_T2 = mean(Mean_T2);

% To view your resulting image in real time
% figure; imslice(T2Image); colormap jet, caxis([0 2]); colorbar;
%figure; imslice(So_image3ParameterFit);
% cmap = colormap('jet'); cmap(:,1) = 0;
% figure; imslice(addT2Combined, 'addT2Combined Map'); colormap(cmap); clim([0 2]); 
% figure; imslice(T2Image, 'T2Star Map'); colormap(cmap); clim([0 2]); 
% figure; imslice(abs(FinalRecon),'Reconstructed Images, No Gaussian Noise Added');

disp('Calculating T2* Completed')

clear FID_Col
clear FID_Col_noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% True T2 That was Prescribed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is to check to and check with final phantom
T2Paren = ones(size(FinParenMask(FinParenMask == 1)));
T2Paren(1:round(size(T2Paren,1)/4,0)) = 0.0004;
T2Paren(round(size(T2Paren,1)/4,0)+1:round(size(T2Paren,1)*2/4,0)) = 0.0004;
T2Paren(round(size(T2Paren,1)*2/4,0)+1:round(size(T2Paren,1)*3/4,0)) = 0.0004;
T2Paren(round(size(T2Paren,1)*3/4,0)+1:end) = 0.0004;

T2Heart = ones(size(FinHeartMask(FinHeartMask == 1)));
T2Heart(1:round(size(T2Heart,1)/4,0)) = 0.0004;
T2Heart(round(size(T2Heart,1)/4,0)+1:round(size(T2Heart,1)*2/4,0)) = 0.0004;
T2Heart(round(size(T2Heart,1)*2/4,0)+1:round(size(T2Heart,1)*3/4,0)) = 0.0004;
T2Heart(round(size(T2Heart,1)*3/4,0)+1:end) = 0.0004;

T2Tissue = ones(size(FinTissueMask(FinTissueMask == 1)));
T2Tissue(1:round(size(T2Tissue,1)/4,0)) = 0.0004;
T2Tissue(round(size(T2Tissue,1)/4,0)+1:round(size(T2Tissue,1)*2/4,0)) = 0.0004;
T2Tissue(round(size(T2Tissue,1)*2/4,0)+1:round(size(T2Tissue,1)*3/4,0)) = 0.0004;
T2Tissue(round(size(T2Tissue,1)*3/4,0)+1:end) = 0.0004;

T2Vasc = ones(size(FinVascMask(FinVascMask == 1)));
T2Vasc(1:round(size(T2Vasc,1)/4,0)) = 0.00316;
T2Vasc(round(size(T2Vasc,1)/4,0)+1:round(size(T2Vasc,1)*2/4,0)) = 0.00316;
T2Vasc(round(size(T2Vasc,1)*2/4,0)+1:round(size(T2Vasc,1)*3/4,0)) = 0.00316;
T2Vasc(round(size(T2Vasc,1)*3/4,0)+1:end) = 0.00316;

% Preallocate memory for T2* variables of different tissues
addT2Lung = zeros(ImSize, ImSize, ImSize);
addT2Heart = zeros(ImSize, ImSize, ImSize);
addT2Tissue = zeros(ImSize, ImSize, ImSize);
addT2Vasc = zeros(ImSize, ImSize, ImSize);

% Find the locations where the binary masks == 1 for each tissue
LungLoc = find(FinParenMask==1);
HeartLoc = find(FinHeartMask==1);
TissueLoc = find(FinTissueMask==1);
VascLoc = find(FinVascMask==1);

% Lastly, add the varying T2* for each tissue in the binary locations
% Changing 0's to 1's to not delete any information about the phantom
addT2Lung(LungLoc) = T2Paren; %addT2Lung(addT2Lung == 0) = Inf;
addT2Heart(HeartLoc) = T2Heart; %addT2Heart(addT2Heart == 0) = Inf;
addT2Tissue(TissueLoc) = T2Tissue; %addT2Tissue(addT2Tissue == 0) = Inf;
addT2Vasc(VascLoc) = T2Vasc; %addT2OutEdge(addT2OutEdge == 0) = Inf;
addT2Combined = addT2Lung + addT2Heart + addT2Tissue + addT2Vasc;

clear addT2Heart addT2Lung addT2Tissue addT2OutEdge arrayLength
clear Im_UTE_TE1 Im_UTE_TE2 Im_UTE_TE3 Im_UTE_TE4 Im_UTE_TE5
clear LungLoc HeartLoc Tissueloc OutEdgeLoc So_image3ParameterFit
clear k m n T2_image3parameterfit T2array_nonlinear SoArray_nonlinear
clear nTE nPts NTE p_smp SampPct sd_N ptsPerPro Dwell
clear difference_percent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Save_SimulationResult,'yes') == 1
    path = uigetdir;
    save_path=[path, '\UTEundersampling_Result1.mat'];
    save(save_path);
    % save('E:\Lab\Ian UTE\3D Simulation_undersampling\undersampling_Result1.mat');
    disp('Workspace has been saved')
else 
    disp('Workspace has not been saved')
end
disp('Simulation completed')
Simulation_Run_Time = toc/60

SimulationOutput.AddedT2.Lung = T2Paren;
SimulationOutput.AddedT2.Heart = T2Heart;
SimulationOutput.AddedT2.Tissue = T2Tissue;
SimulationOutput.AddedT2.OuterEdge = T2Vasc;

SimulationOutput.TrueT2Mask = addT2Combined;

SimulationOutput.BinaryMask.Paren = FinParenMask;
SimulationOutput.BinaryMask.Heart = FinHeartMask;
SimulationOutput.BinaryMask.Tissue = FinTissueMask;
SimulationOutput.BinaryMask.Vasc = FinVascMask;
SimulationOutput.BinaryMask.Noise = FinNoiseMask;

SimulationOutput.SimulationParams.MatrixSize = ImSize;
SimulationOutput.SimulationParams.SamplingPercent = SampPct_all;
SimulationOutput.SimulationParams.PrescribedProj = NPro_FS;
SimulationOutput.SimulationParams.ActualProj = NPro_Act;
SimulationOutput.SimulationParams.GaussianNoiseSD = noisesd;
SimulationOutput.SimulationParams.OrignalPhantom = MyPhantom;
SimulationOutput.SimulationParams.TruePhantom = PrePhantom;
% SimulationOutput.SimulationParams.DecayPhantom = MyPhantom_T2s;

SimulationOutput.ImgRecon = FinalRecon;
SimulationOutput.AllMeanT2 = Mean_T2;
SimulationOutput.AllMeanSNR = SNR_vec_n_all;

SimulationOutput.T2Fit.T2Img = T2Image;
SimulationOutput.T2Fit.S0Img = SoImage;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SNR_vec_n signal_vec noise_vec
clear T2Heart T2Lung T2Tissue T2OutEdge TissueLoc Simulation_Run_Time
clear full_sampl
