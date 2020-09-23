%Zijian Dong
%CBF calculation

%Read data.
perfusion_8 = niftiread('/Users/zijiandong/Desktop/research/B50032/8_1_pCASL_EPIv8_pCASL_EPIv8_perf_0offset.nii.gz');
perfusion_9 = niftiread('/Users/zijiandong/Desktop/research/B50032/9_1_pCASL_EPIv8_pCASL_EPIv8_perf_0offset.nii.gz');

%Interleave perfusion image.
perfusion = zeros([100,100,30,60]);
for i = 1:size(perfusion_8,4)
    for j = 1:size(perfusion_8,3)
        perfusion(:,:,2*j-1,i)=perfusion_8(:,:,j,i);
        perfusion(:,:,2*j,i)=perfusion_9(:,:,j,i);
    end
end

%Separate control and label.
con = perfusion(:,:,:,1:2:59);
lab = perfusion(:,:,:,2:2:60);
deltaM = sum(abs(con-lab),4)./30;

% Parameters
T1vals = T1vals.*1000;
TR = 15703;
M0 = (sum(con,4)./30)./(1-exp(-TR./T1vals));
lambda= 0.9;
w= 0.2;
T1b = 2.3;
IE = 0.663;
tau = 3000;

%Calculation
CBF_vals=(6000*lambda.*deltaM.*exp(w./T1b))./(2.*IE.*T1vals.*M0.*(1.-exp(-tau./T1vals)))*1000;
CBF_vals(CBF_vals>500)=0;

%Smooth
mask = CBF_vals;
mask(CBF_vals>0)=1;
BW = bwareaopen(mask,16,8);
CBF_vals = CBF_vals.*BW;
gf = fspecial('gaussian',[3 3],1);
CBF_vals = imfilter(CBF_vals,gf);

