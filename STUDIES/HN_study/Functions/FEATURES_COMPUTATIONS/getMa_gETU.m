function [Ma_gETU]= getMa_gETU(ROIonlyPET,a)

nVoxels = sum(~isnan(ROIonlyPET(:)));
ROIonlyPET = ROIonlyPET .^a;

Ma_gETU = (sum(ROIonlyPET(~isnan(ROIonlyPET)))/nVoxels)^(1/a);
end