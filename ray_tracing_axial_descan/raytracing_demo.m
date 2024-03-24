
%% Demo for ray tracing of the axial-descan system, unit is mm except specifically defined 
% input parameters 
focallength1 = [5,200,300,10];      % focal lengths of the lenses before scan mirror
focallength2 = [10,300];            % focal lengths of the lenses after LFA mirror
d1 = [5,100,500,50,50,543,10];      % distances between the optical elements before LFA mirror
d_rfobj = -0.018;                   % offset of LFA mirror to the nominal focal plane of the remote objective
d2 = [10+d_rfobj,300];              % distances between the optical elements after LFA mirror, ignore the repeated distances in the unfold path, which are d1([4,5,6])
d0 = 0.045;                         % offset of the sample position to the nominal focal plane of the detection objective
nmed = 1.33;                        % refractive index of the sample medium
mirror_angle = [-1,0.3733];         % [slope of BS, slope of mirror relative to the BS coordinates], BS coordinate is defined by the x-axis along the splitting surface
aperturesize = 50;                  % diameter of optical elements
pupilsize = 8;                      % diameter of the pupil of the detection objective
fov = 0.1;                          % size of the object
ccdpixelsize = 6.5;                 % pixel size of the camera, unit: micron
labelflag = 1;                      % 0 or 1, show labels for beam sizes, beam size is calculated at each lens position.

% output parameters
% - dia_lens:                       diameter of the beam size at the each lens position
% - dia_image:                      diameter of the image size at each image plane
% - magnification:                  total magnification
% - pixelsize:                      pixel size at the sample space, unit: nm 


[dia_lens,dia_image,magnification,pixelsize] = raytracing_rf_descan(focallength1,focallength2,d1,d2,d0,...
                                                                       nmed,mirror_angle,pupilsize,aperturesize,...
                                                                       fov,ccdpixelsize,labelflag);
