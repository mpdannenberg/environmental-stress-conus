% Calculate PET and VPD from PRISM data following FAO Penman-Monteith
% method (Allen et al. 1998)

cd prism; % Folder with PRISM data
f1 = matfile('PRISM_TDMEAN');
f2 = matfile('PRISM_TMAX');
f3 = matfile('PRISM_TMIN');

load prism_latlon;
year = f1.year;

cd ../data;
[ELEV, R] = geotiffread('us_dem_4km.tif');
ELEV(ELEV<-1000) = NaN;

COAST = geotiffread('us_CoastalBoundary_4km.tif');
COAST(COAST<-1000) = NaN;

nx = length(PRISMlon);
ny = length(PRISMlat);

[ny,nx,nt,nm] = size(f1,'TDMEAN');

cd ../prism; % Folder with PRISM data

f4 = matfile('PRISM_PET.mat', 'Writable',true);
f5 = matfile('PRISM_VPD.mat', 'Writable',true);

tic
idx = 1;
temp = NaN(1,1,nt,nm);
for i = 1:ny
    for j = 1:nx
        
        tmin = squeeze(f3.TMIN(i, j, :, :));
        tmax = squeeze(f2.TMAX(i, j, :, :));
        tdmean = squeeze(f1.TDMEAN(i, j, :, :));
        
        if sum(sum(~isnan(tmin))) > 0
            [pet, ~, vpd, ~, ~] = fao_pm(tmax', tmin', tdmean', PRISMlat(i), ELEV(i, j), COAST(i, j), year);
            temp(1,1,:,:) = pet';
            if isreal(pet)
                f4.PET(i, j, :, :) = temp;
            else
                f4.PET(i, j, :, :) = real(temp);
            end
            temp(1,1,:,:) = vpd';
            f5.VPD(i, j, :, :) = temp;
        end
        
        
        
        if toc > 60
            str = strcat(num2str(100*idx/(nx*ny)), '% Complete');
            disp(str)
            tic
        end
        
        idx = idx+1;
        
    end
end

cd ..;

