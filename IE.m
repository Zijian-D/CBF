%Zijian DOng
%IE calculation(manully define the roi of carotid)
FcFlash = niftiread(['/Users/zijiandong/Desktop/research/B50032/10_1_'...
'pCASL_FcFLASHv2_pCASL_FcFLASHv2_effficiencym_no4']);

realLab=FcFlash(47:48,37:38,1,1);
imgLab=FcFlash(47:48,37:38,1,2);
realCon=FcFlash(49:50,36:37,1,3);
imgCon=FcFlash(49:50,36:37,1,4);
IE1 = (sqrt((realLab-realCon).^2+(imgLab-imgCon).^2))./(2.*sqrt(realCon.^2+imgCon.^2));

realLab=FcFlash(69:70,34:35,1,4);
imgLab=FcFlash(69:70,34:35,1,3);
realCon=FcFlash(68:69,33:34,1,2);
imgCon=FcFlash(68:69,33:34,1,1);
IE2 = (sqrt((realLab-realCon).^2+(imgLab-imgCon).^2))./(2.*sqrt(realCon.^2+imgCon.^2));

ie=(sum(sum(IE1))+sum(sum(sum(IE2))))/8;