%Zijian Dong
%CBF calculation

% %Read data.
% perfusion_8 = niftiread('/Users/zijiandong/Desktop/research/B50032/8_1_pCASL_EPIv8_pCASL_EPIv8_perf_0offset.nii.gz');
% perfusion_9 = niftiread('/Users/zijiandong/Desktop/research/B50032/9_1_pCASL_EPIv8_pCASL_EPIv8_perf_0offset.nii.gz');
% 
% %Interleave perfusion image.
% perfusion = zeros([100,100,30,60]);
% for i = 1:size(perfusion_8,4)
%     for j = 1:size(perfusion_8,3)
%         perfusion(:,:,2*j-1,i)=perfusion_8(:,:,j,i);
%         perfusion(:,:,2*j,i)=perfusion_9(:,:,j,i);
%     end
% end
cd('/Users/zijiandong/Desktop/research')
mkdir CBF

maps_list = dir(['/Users/zijiandong/Desktop/research/warped_T1/','*.nii.gz']);
for i = 1:length(maps_list)
    T1vals = niftiread(['/Users/zijiandong/Desktop/research/warped_T1/',maps_list(i).name]);
    [index_3,index_4] = xlsread('/Users/zijiandong/Desktop/research/Inversion_Efficiencies.xlsx',1,'A2:M58');
    for j = 1:size(index_4,1)
        if strcmp(extractBefore(index_4{j,1},':'),extractBefore(maps_list(i).name,'W'))
            ine = index_3(j);
        end
    end
    file_name = ['/Users/zijiandong/Desktop/research/interleaved/Animal_',extractBefore(maps_list(i).name,'W'),'_perfusion.nii.gz'];
    perfusion = niftiread(file_name);
    info = niftiinfo(file_name);
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
    tau = 3000;

%Calculation
    CBF_vals=(6000*lambda.*deltaM.*exp(w./T1b))./(2.*ine.*T1vals.*M0.*(1.-exp(-tau./T1vals)))*1000;
    CBF_vals(isnan(CBF_vals))=0;
    CBF_vals(CBF_vals==Inf)=0;
    CBF_vals(CBF_vals>1000)=0;
%     ori = CBF_vals;
%     ori(ori>500) = 0;
%     CBF_vals(CBF_vals>400)=400;
%     BW = imbinarize(ori,0);
%     for k = 1:30
%         fill = imfill(BW(:,:,k),'hole');
%         CBF_vals(:,:,k) = CBF_vals(:,:,k).*fill;
%     end
    
%     CBF_vals(CBF_vals == Inf)=0;
%     CBF_vals(CBF_vals>500)=500;
%Smooth
%     mask = niftiread('/Users/zijiandong/Desktop/research/reg_result/mask.nii.gz');
%     mask(mask ~= 0) = 1;
%     CBF_vals_m = CBF_vals.*mask;
%     mask = CBF_vals;
%     mask(CBF_vals>0)=1;
%     BW = bwareaopen(mask,16,8);
%     CBF_vals = CBF_vals.*BW;
    gf = fspecial('gaussian',[3 3],1);
    CBF_vals = imfilter(CBF_vals,gf);
    
    cd('/Users/zijiandong/Desktop/research/CBF/')
    nii = load_untouch_nii(['/Users/zijiandong/Desktop/research/interleaved/Animal_',extractBefore(maps_list(i).name,'W'),'_perfusion.nii.gz']);
    nii.img = CBF_vals;
    new_dims=ones([1 length(nii.hdr.dime.dim)]);
    num_dims=length(size(nii.img));
    new_dims(1)=num_dims;
    new_dims(2:(num_dims+1))=size(nii.img);
    nii.hdr.dime.dim=new_dims;
    save_untouch_nii(nii,[extractBefore(maps_list(i).name,'W'),'.nii.gz'])
    display(extractBefore(maps_list(i).name,'W'))
end            
