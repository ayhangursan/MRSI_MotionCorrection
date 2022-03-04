%% Simulation enviroment for testing methods

% This is a 2D simulation
% A body slice shaped digital phantom is generated using Shepp-Logan Phantom.
% Stomach, liver and subcutaeous fat compartments are generated on this
% image.

% Deuterium compounds are used in the simulation. Linewidths and chemical
% shifts are used for these metabolites:
% 1-HDO
% 2-Glucose
% 3-Glx
% 4-Lipid


% Ayhan Gursan, UMC Utrecht 2021
% a.gursan@umcutrecht.nl


%% Load Shepp-Logan Phantom
clear;clc;
% close all
PhantomSize=64;
SLPhantom=phantom(PhantomSize);
SLPhantom=imrotate(SLPhantom,-90);
SLPhantom=round(SLPhantom,2);
%% Show phantom image
figure('WindowState','maximized')
imagesc(SLPhantom)
daspect([1 1 1]); colorbar
%% Generate masks for body compartements
% Subcutaneous Lipid
Mask.Lipid=zeros(size(SLPhantom));
Mask.Lipid(SLPhantom==1)=1;
% Stomach
Mask.Stomach=zeros(size(SLPhantom));
Mask.Stomach(SLPhantom==0.3)=1;
Mask.Stomach(:,1:PhantomSize/2)=0;
% Liver
Mask.Liver=zeros(size(SLPhantom));
Mask.Liver(SLPhantom==0)=1;
InCirc=round(PhantomSize*0.273);
Mask.Liver(1:InCirc,:)=0;
Mask.Liver(PhantomSize-InCirc:end,:)=0;
Mask.Liver(:,1:InCirc)=0;
Mask.Liver(:,PhantomSize-InCirc:end)=0;
Mask.Liver=circshift(Mask.Liver,round(-PhantomSize/6),2);
MaskMatrix(:,:,1)=Mask.Lipid;
MaskMatrix(:,:,2)=Mask.Stomach;
MaskMatrix(:,:,3)=Mask.Liver;

% Whole Body for water signal
Mask.Water=imfill(Mask.Lipid,'holes');
MaskMatrix(:,:,4)=Mask.Water;
% %% Show the mask
% figure('WindowState','maximized')
% subplot(2,3,[1,4])
% imagesc(SLPhantom)
% daspect([1 1 1]); colorbar;
% title('Phantom image')
%
% subplot(2,3,2)
% imagesc(squeeze(MaskMatrix(:,:,1)))
% daspect([1 1 1]); colorbar;
% title('Lipid mask image')
%
% subplot(2,3,3)
% imagesc(squeeze(MaskMatrix(:,:,2)))
% daspect([1 1 1]); colorbar;
% title('Stomach mask image')
%
% subplot(2,3,5)
% imagesc(squeeze(MaskMatrix(:,:,3)))
% daspect([1 1 1]); colorbar;
% title('Liver mask image')
%
% subplot(2,3,6)
% imagesc(squeeze(MaskMatrix(:,:,4)))
% daspect([1 1 1]); colorbar;
% title('Water on full body')

%% Move objects
TR=0.333;
seconds=60;
time=round((0:1*(PhantomSize*PhantomSize)-1));
Movingphantom=zeros([size(SLPhantom) numel(time)]);
omeg=20;
livermotionrange=0.07*PhantomSize;
stomachmotionrange=0;
LiverMotionVector=round(livermotionrange*sind(omeg*time(1:numel(time))));
StomachMotionVector=round(stomachmotionrange*sind(omeg*time(1:numel(time))));

for frame=1:numel(time)
    Movingphantom(:,:,frame)=Mask.Lipid+circshift(Mask.Liver,LiverMotionVector(frame),1)+circshift(Mask.Stomach,StomachMotionVector(frame),1);
end

implay(Movingphantom,24)
%% Reconstruct moving data from acquisition
Acquiredkspace=zeros(PhantomSize,PhantomSize);
kspaceduringAcq=zeros(size(Movingphantom));
for frame=1:numel(time)
    kspaceduringAcq(:,:,frame)=fftshift(fft(fftshift(fft(Movingphantom(:,:,frame),[],1),1),[],2),2);
    kpace2D=kspaceduringAcq(:,:,frame);
    Acquiredkspace(frame)=kpace2D(frame);
end

%% Effect of motion
figure('WindowState','maximized')
subplot(1,2,1)
imagesc(real(Movingphantom(:,:,1)));
daspect([1 1 1])

subplot(1,2,2)
imagesc(real(ifft(ifftshift(ifft(ifftshift(Acquiredkspace,1),[],1),2),[],2)));
daspect([1 1 1]);title('Pulse- acq')

%% Motion correction
% Apply phase shift in k-space to correct image
kx1 = (PhantomSize-1)/PhantomSize;
ky1 = (PhantomSize-1)/PhantomSize;

[kx,ky] = ndgrid(linspace(-kx1,kx1,PhantomSize),linspace(-ky1,ky1,PhantomSize));

Correctedkspace=zeros(PhantomSize,PhantomSize);

for frame=1:numel(time)
    Correctedlines=kspaceduringAcq(:,:,frame).*exp(pi*1i*kx.*LiverMotionVector(frame) + pi*1i*ky.*LiverMotionVector(frame)*0);
    Correctedkspace(frame)=Correctedlines(frame);
end

%% Effect of motion correction
figure('WindowState','maximized')
% imagesc(real(fftshift(fft(fftshift(fft(Acquiredkspace(:,:,1),[],1),1),[],1),1)))
subplot(1,3,1)
imagesc(real(Movingphantom(:,:,1)));
daspect([1 1 1])

subplot(1,3,2)
imagesc(real(ifft(ifftshift(ifft(ifftshift(Acquiredkspace,1),[],1),2),[],2)));
daspect([1 1 1]);title('Pulse-acq')

subplot(1,3,3)
imagesc(real(ifft(ifftshift(ifft(ifftshift(Correctedkspace,1),[],1),2),[],2)));
daspect([1 1 1]);title('Corrected')

%% Apply motion correction on voxel base
% A reference and dynamic motion fields are required to estimate voxel
% motion during acquisition.
Mask.Moving=zeros(size(Mask.Liver));
Mask.Moving(20:43,11:33)=1;
% ones(size(Mask.Liver))
% MotionDynamics=permute(repmat(Mask.Moving,[1 1 numel(time)]),[3 1 2]).*repmat(LiverMotionVector.',[1 PhantomSize PhantomSize]);
% MotionDynamics=permute(MotionDynamics,[2 3 1]);
% PhaseCorr=exp(pi*1i*kx.*MotionDynamics + pi*1i*ky.*MotionDynamics*0); %
% PhaseCorr=ones([size(SLPhantom) numel(time)]); % DEBUG
Correctedkspace_2=zeros([PhantomSize,PhantomSize, numel(time)]);
%
tic
for pixel=1:numel(time)
    [x, y]=ind2sub(size(SLPhantom),pixel);
    if Mask.Moving(x,y)==0
        for frame=1:numel(time)
                [kxind, kyind]=ind2sub(size(SLPhantom),frame);
            Correctedkspace_2(kxind,kyind,pixel)=kspaceduringAcq(kxind,kyind,frame);
        end
    else
        for frame=1:numel(time)
            [kxind, kyind]=ind2sub(size(SLPhantom),frame);
            Correctedlines=kspaceduringAcq(:,:,frame).*exp(pi*1i*kx.*LiverMotionVector(frame) + pi*1i*ky.*LiverMotionVector(frame)*0);
            Correctedkspace_2(kxind,kyind,pixel)=Correctedlines(kxind,kyind);
        end
    end
end
toc
%%
% Voxel based recon
VoxCorrImag=zeros(size(Correctedkspace_2));
VoxCorrImag=ifft(ifftshift(ifft(ifftshift(Correctedkspace_2,1),[],1),2),[],2);
CombImag=zeros(size(SLPhantom));
% figure('WindowState','maximized')

for frame=1:numel(time)
    [x, y]=ind2sub(size(SLPhantom),frame);
    CombImag(x,y)=VoxCorrImag(x,y,frame);
end
figure
imagesc(real(squeeze(CombImag)))

%% Show all results
figure('WindowState','maximized')
subplot(2,2,1)
imagesc(real(Movingphantom(:,:,1)));
daspect([1 1 1])

subplot(2,2,2)
imagesc(real(ifft(ifftshift(ifft(ifftshift(Acquiredkspace,1),[],1),2),[],2)));
daspect([1 1 1]);title('Pulse-acq')

subplot(2,2,3)
imagesc(real(ifft(ifftshift(ifft(ifftshift(Correctedkspace,1),[],1),2),[],2)));
daspect([1 1 1]);title('Corrected')

subplot(2,2,4)
imagesc(real(squeeze(CombImag)))
daspect([1 1 1]);title('Voxel based recon')

%% Spectroscopy related lines
% %% Introduce chemical shifts (ppm)
% % HDO signal is on-resonance
% FreqOffset=4.7;
% ChemicalShift.Lipid=1.3-FreqOffset;
% ChemicalShift.Glx=2.25-FreqOffset;
% ChemicalShift.Glucose=3.8-FreqOffset;
% ChemicalShift.HDO=4.7-FreqOffset;
% %% Introduce lorentzian lineshape (Hz)
% Linewidth.Lipid=30;
% Linewidth.Glx=30;
% Linewidth.Glucose=30;
% Linewidth.HDO=30;
% %% Introduce amplitude
% % Change metabolite amplutides in here
% Amplitude.Lipid=10;
% Amplitude.Glx=1;
% Amplitude.Glucose=8;
% Amplitude.HDO=10;
% %% Simulation parameters
% SimParam.BW=10000;
% SimParam.NP=2048;
% SimParam.B0=7;
% SimParam.GyromagRatio=6.53569*10^6;
% SimParam.Freq=SimParam.GyromagRatio*SimParam.B0;
% SimParam.time=single([0:SimParam.NP-1]*(1/SimParam.BW));
% SimParam.ppmwindow=((SimParam.Freq/SimParam.BW)^-1)*10^6;
% SimParam.xaxis=linspace(-SimParam.ppmwindow/2,SimParam.ppmwindow/2,SimParam.NP)+FreqOffset;
%
% %% Generate lineshape map
% tic
% LineShape.Lipid=single(permute(repmat(Mask.Lipid,[1 1 SimParam.NP]),[3 1 2]).*exp(-abs(SimParam.time) * Linewidth.Lipid.' * pi).');
% LineShape.Glx=single(permute(repmat(Mask.Liver,[1 1 SimParam.NP]),[3 1 2]).*exp(-abs(SimParam.time) * Linewidth.Glx.' * pi).');
% LineShape.Glucose=single(permute(repmat(Mask.Liver+Mask.Stomach,[1 1 SimParam.NP]),[3 1 2]).*exp(-abs(SimParam.time) * Linewidth.Glucose.' * pi).');
% LineShape.HDO=single(permute(repmat(Mask.Water,[1 1 SimParam.NP]),[3 1 2]).*exp(-abs(SimParam.time) * Linewidth.HDO.' * pi).');
% toc
% disp(['Lineshapes are generated with BW=',num2str(SimParam.BW),'  Number of points=',num2str(SimParam.NP),'.'])
% %% Generate metabolite maps
% tic
% Metabmap.Lipid=single(Amplitude.Lipid.*LineShape.Lipid.*exp(2*pi*1i*(ChemicalShift.Lipid)*(SimParam.Freq/(10^6))*SimParam.time).');
% Metabmap.Glx=single(Amplitude.Glx.*LineShape.Glx.*exp(2*pi*1i*(ChemicalShift.Glx)*(SimParam.Freq/(10^6))*SimParam.time).');
% Metabmap.Glucose=single(Amplitude.Glucose.*LineShape.Glucose.*exp(2*pi*1i*(ChemicalShift.Glucose)*(SimParam.Freq/(10^6))*SimParam.time).');
% Metabmap.HDO=single(Amplitude.HDO.*LineShape.HDO.*exp(2*pi*1i*(ChemicalShift.HDO)*(SimParam.Freq/(10^6))*SimParam.time).');
%
% Metabmap.All= Metabmap.Lipid + Metabmap.Glx + Metabmap.Glucose + Metabmap.HDO;
% Metabmap.All_Spec=fftshift(fft(Metabmap.All,[],1),1);
% toc
% disp(['Metabolite maps are generated at high spatial resolution.'])
% %% Check signal from specific voxel at High spatial resolution scan
% figure;
% AP=90;
% RL=90;
% plot(SimParam.xaxis,real(Metabmap.All_Spec(:,AP,RL)))
% xlim([0 10])
% set(gca,'XDir','reverse')
% %% Downscale matrix for low spatial resolution
% SimulatedDataSet.HighRes_kspace=fftshift(ifftn(ifftshift(Metabmap.All)));
% SimulatedDataSet.HighRes_kspace=fftshift(fft(SimulatedDataSet.HighRes_kspace,[],1),1);
%
% K_X=16;
% K_Y=16;
%
% SimulatedDataSet.LowRes_kspace=SimulatedDataSet.HighRes_kspace(:,(PhantomSize-K_X)/2+1:(PhantomSize+K_X)/2,(PhantomSize-K_Y)/2+1:(PhantomSize+K_Y)/2);
% SimulatedDataSet.LowRes_image=ifftshift(fftshift(fftn(SimulatedDataSet.LowRes_kspace)),1);
% SimulatedDataSet.LowRes_image=fftshift(fft(Phasecorrection(ifft(ifftshift(SimulatedDataSet.LowRes_image,1),[],1)),[],1),1); % Phase correction on time domain
%
% [idealhammingwindow,Correctionmask,SimulatedDataSet.LowRes_kspaceHamming]=Hammingfilter_2D(SimulatedDataSet.LowRes_kspace,ones(K_X,K_Y));
% SimulatedDataSet.LowRes_imageHamming=ifftshift(fftshift(fftn(SimulatedDataSet.LowRes_kspaceHamming)),1);
% SimulatedDataSet.LowRes_imageHamming=fftshift(fft(Phasecorrection(ifft(ifftshift(SimulatedDataSet.LowRes_imageHamming,1),[],1)),[],1),1); % Phase correction on time domain
%
% %%
% figure
% subplot(2,2,1)
% imagesc(squeeze(sum(abs(SimulatedDataSet.HighRes_kspace),1)));colormap gray
% daspect([1 1 1]);title('K-space of high resolution image')
% subplot(2,2,2)
% imagesc(squeeze(sum(abs(SimulatedDataSet.LowRes_kspace),1)));colormap gray
% daspect([1 1 1]);title('K-space of low resolution image')
%
% subplot(2,3,4)
% imagesc(squeeze(sum(abs(Metabmap.All),1)));colormap gray
% daspect([1 1 1]);title('High resolution image')
% subplot(2,3,5)
% imagesc(squeeze(sum(abs(SimulatedDataSet.LowRes_image),1)));colormap gray
% daspect([1 1 1]);title('Low resolution image')
% subplot(2,3,6)
% imagesc(squeeze(sum(abs(SimulatedDataSet.LowRes_imageHamming),1)));colormap gray
% daspect([1 1 1]);title('Low resolution image after Hamming filter')
%
% %% Check signal from specific voxel - Low spatial resolution
% figure;
% AP=7;
% RL=7;
% plot(SimParam.xaxis,real(SimulatedDataSet.LowRes_image(:,AP,RL)))
% hold on
% plot(SimParam.xaxis,real(SimulatedDataSet.LowRes_imageHamming(:,AP,RL)))
% hold off
% xlim([0 10])
% set(gca,'XDir','reverse')
% legend('Low Spatial Resolution','Low Res + Hamming')
% title(['AP:',num2str(AP),'  RL:',num2str(RL)])
%
% %%
% function [idealhammingwindow,Correctionmask,filteredCSI]=Hammingfilter_2D(fiddata,acqpattern)
% dims=size(fiddata);
% grid=dims(2:end);
% for k=1:prod(grid)
%     [x,y]= ind2sub(grid,k);
%     filterfuncx=0.54+0.46*cos(2*pi*x/grid(1));
%     filterfuncy=0.54+0.46*cos(2*pi*y/grid(2));
%     hammingwindow(k)=filterfuncx*filterfuncy;
% end
%
% hammingwindow=fftshift(reshape(hammingwindow,grid));
%
% idealhammingwindow=circshift(hammingwindow,1,1);
% idealhammingwindow=circshift(idealhammingwindow,1,2);
%
% Correctionmask=idealhammingwindow./acqpattern;
% % In datasets acquired with old k-space weighting patch shutter leaves some
% % points with no line acquired(acqpattern=0). This dividing by zero may lead to inf points in correction mask.
% % To avoid that check if acqpattern has zeros, if so correctionmask(correctionmask==inf)=0
% if numel(find(~acqpattern))>0
%     disp('Non-acquired k-space points in the dataset.')
%     Correctionmask(isinf(Correctionmask))=0;
% end
%
% dimmatchedmask=repmat(Correctionmask,[1 1 1 dims(1)]);
% dimmatchedmask=permute(dimmatchedmask,[4 1 2 3]);
%
% filteredCSI=fiddata.*dimmatchedmask;
%
% end
