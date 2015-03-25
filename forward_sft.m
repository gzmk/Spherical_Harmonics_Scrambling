% Read in a light field, sample it and forward transform using spherical
% harmonics using S2kit (http://rodgers.org.uk/software/s2kit)

% author @ gizem
% March 2015


lf = imread('envmap.bmp','bmp'); % make sure it read in uint8 format
% convert image to grayscale
lf_bw = squeeze(lf(:,:,1)*0.2989 + lf(:,:,2)*0.5870 + lf(:,:,3)*0.1140);
% downsample image by 8
lf_bw_sm = imresize(lf_bw, [128 256]);

[nlat, nlong] = size(lf_bw_sm);
% latitude -> theta is the angle from the z-axis (0,pi)
% longitude -> phi is the angle in the xy-plane from the x-axis (0,2*pi)
latstep = pi/nlat;
longstep = 2*pi/nlong;

lats = (latstep/2):latstep:(pi-(latstep/2));
longs = 0:longstep:(2*pi-longstep);
[phi, theta] = meshgrid(longs,lats);

% gridpts = [phi(:) theta(:)];
% wts = getVoronoiWeights(gridpts);
% Nord = 50;
% xform = leastSquaresSHT(Nord,lfgsm(:),gridpts,'complex',wts);

j = 8;
m = 4;
fprintf('\n\n\n****************\nTesting with j=%g, m=%g:\n',j,m)
SphHRep=FSTRep2SphHRep(FST_semi_fly_mex(Ylm(j,m,theta,phi)));
[nonzeroelements]=find(abs(SphHRep)>1e-14);
if numel(nonzeroelements)>0
    for idx=1:numel(nonzeroelements)
        [foundj,foundm]=idx2jm(nonzeroelements(idx));
        fprintf('Found: j=%g, m=%g. amplitude=%g\n',foundj,foundm, ...
            SphHRep(nonzeroelements(idx)))
    end
else
    fprintf('All coefficients are negligible!\n')
end
fprintf('----------------------\n');
