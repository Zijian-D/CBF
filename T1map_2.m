%Zijian Dong
%T1map generation for B50032(Interleaving and resampling for T1image included).
%Running time 218s(2.6 GHz 6-Core Intel Core i7)

%Read data
T1_11 = niftiread(['/Users/zijiandong/Desktop/research/B50032/11_1_Pe'...
'rfusion_FAIR_EPI.nii.gz']);
T1_12 = niftiread(['/Users/zijiandong/Desktop/research/B50032/12_1_Pe'...
'rfusion_FAIR_EPI.nii.gz']);
TI = [30.000; 45.400; 68.704; 103.972; 157.342; 238.110; 360.337; 545.306; 825.224; 1248.83; 1889.882; 2859.993; 4328.1; 6549.808; 9911.969; 15000]./1000;

%Interleave T1image,combine two 100x100x20x16 into one 100x100x40x16.
T1_inL = zeros([100,100,40,16]);
for i = 1:size(T1_11,4)
    for j = 1:size(T1_11,3)
        T1_inL(:,:,2*j-1,i)=T1_11(:,:,j,i);
        T1_inL(:,:,2*j,i)=T1_12(:,:,j,i);
    end
end

%Resample combined T1image to 100x100x30x16.(Interpolation).
T1 = zeros([100,100,30,16]);
for i = 1:size(T1_inL,4)
    T1(:,:,:,i)=imresize3(T1_inL(:,:,:,i),[100 100 30]);
end
    
%Remove background(Thresholding),medium filtering.
T1(T1<10)=0;
for m = 1:16
    for n = 1:30
    T1(:,:,n,m)=medfilt2(T1(:,:,n,m));
    end
end

% flip
M_z = reshape(T1,100*100*30,16);
[min,index]=min(M_z,[],2);
for k = 1:100*100*30
    for j = 1:index(k)-1
        M_z(k,j) = -M_z(k,j);
    end
end

%fitting
parfor i = 1:100*100*30
%     xdata=TI.';
%     ydata=double(M_z(i,:));
%     fun = @(x,xdata)x(1)*(1-2*0.7*exp(-xdata/x(2)));
%     x0 = [0,0];
%     x = lsqcurvefit(fun,x0,xdata,ydata);
%     T1vals(i,1)=x(2);
%     M0vals(i,1)=x(1);
%     disp(i)
%     disp(x(2))
%     figure,clf
%     plot(TI,M_z(i,:),'ko',TI,fun(x,TI),'b-')
%     pause
% end
    fun=@(x)sseval(x, TI, M_z(i,:).'); 
    x0=[0 0]; 
    options=optimset('MaxFunEvals',10000); 
    curveFit=fminsearch(fun, x0, options); 
    T1vals(i,1)=curveFit(2); 
    M0vals(i,1)=curveFit(1); 
    disp(i)
end

%Remove extreme values,apply Gaussian filtering to smooth.
T1vals = reshape(T1vals,100,100,30);
T1vals(T1vals>4)=0;
gf = fspecial('gaussian',[3 3],1);
T1vals = imfilter(T1vals,gf);
M0vals = reshape(M0vals,100,100,30);

%Cost function
function [sse]= sseval(x, TIdata, MZdata)
T1=x(2);
M0=x(1);
LMmod=M0.*(1-2.*0.663.*exp(-TIdata./T1)); 
sse=sum((MZdata-LMmod).^2); 
end
