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

maps_list = dir(['/Users/zijiandong/Desktop/research/T1maps/','*mat']);
for i = 1:length(maps_list)
    mapfile = load(['/Users/zijiandong/Desktop/research/T1maps/',maps_list(i).name]);
    T1vals = mapfile.T1vals;
    [index_3,index_4] = xlsread('/Users/zijiandong/Desktop/research/Inversion_Efficiencies.xlsx',1,'A2:M49');
    for j = 1:size(index_4,1)
        if strcmp(extractBefore(index_4{j,1},':'),extractBefore(maps_list(i).name,':'))
            IE = index_3(j);
        end
    end
    file_name = ['/Users/zijiandong/Desktop/research/interleaved/Animal_',extractBefore(maps_list(i).name,':'),'_perfusion.nii.gz'];
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
    CBF_vals=(6000*lambda.*deltaM.*exp(w./T1b))./(2.*IE.*T1vals.*M0.*(1.-exp(-tau./T1vals)))*1000;
    CBF_vals(CBF_vals>500)=0;

%Smooth
%     mask = niftiread('/Users/zijiandong/Desktop/research/reg_result/mask.nii.gz');
%     mask(mask ~= 0) = 1;
%     CBF_vals_m = CBF_vals.*mask;
    gf = fspecial('gaussian',[3 3],1);
    CBF_vals = imfilter(CBF_vals,gf);
    
    cd('/Users/zijiandong/Desktop/research/CBF/')
    save(extractBefore(maps_list(i).name,':'),'CBF_vals')
%     file = make_nii(CBF_vals);
%     save_nii(file,[extractBefore(maps_list(i).name,':'),'.nii.gz'])
    display(extractBefore(maps_list(i).name,'.'))
end            
