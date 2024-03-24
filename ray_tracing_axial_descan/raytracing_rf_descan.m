function [dia_lens,dia_image,magnification,pixelsize] = raytracing_rf_descan(f1,f2,d1,d2,d0,nmed,mirror_angle,pupilsize,aperturesize,fov,ccdpixelsize,labelflag)
%raytracing_rf_descan  Ray tracing of the axial-descan system.
% [dia_lens,dia_image,magnification,pixelsize] = raytracing_rf_descan(f1,f2,d10,d20,d0,nmed,mirror_angle,pupilsize,aperturesize,fov,ccdpixelsize,labelflag)
%
% This function perform ray tracing of the axial-descan system.
%
%
% INPUTS:
%   f1:                 focal lengths of the lenses before scan mirror (mm)
%   f2:                 focal lengths of the lenses after LFA mirror (mm)
%   d1:                 distances between the optical elements before LFA mirror (mm)
%   d2:                 distances between the optical elements after LFA mirror, ignore the repeated distances in the unfold path (mm)
%   d0:                 offset of the sample position to the nominal focal plane of the detection objective (mm)
%   nmed:               refractive index of the sample medium
%   mirror_angle:       [slope of BS, slope of mirror relative to the BS coordinates], BS coordinate is defined by the x-axis along the splitting surface
%   pupilsize:          diameter of the pupil of the detection objective (mm)
%   aperturesize:       diameter of optical elements (mm)
%   fov:                size of the object (mm)    
%   ccdpixelsize:       pixel size of the camera, unit: micron (um) 
%   labelflag:          0 or 1, show labels for beam sizes, beam size is calculated at each lens position.
%
%
% OUTPUTS:
%   dia_lens:           diameter of the beam size at the each lens position
%   dia_image:          diameter of the image size at each image plane
%   magnification:      total magnification
%   pixelsize:          pixel size at the sample space, unit: nm 
%
% Author:
%    Sheng Liu, 2024

d_m = d1(4:5);                              % distances to BS and mirrors
d1 = d1([1,2,3,6:end]);                     % distances to lenses before LFA
d2 = [d2(1),d_m(1),d2(2)];                  % distances to lenses after LFA
[f1] = caldis(f1);                          % append zero 
f1(1) = f1(1)*nmed;                         % focal length in the sample medium
d1(1) = f1(1)+d0;                           % distance from object to the collection objective
n = [nmed,ones(1,numel(d1))];               % refractive indices between each optical elements
pupilpoints = [pupilsize/2,0,-pupilsize/2]; % define three points at the pupil plane
fieldpoints = [fov/2,0,-fov/2];             % define three points at the object plane
cm = {'r','g','b'};                         % define a color for each field points
% draw rays
h = figure;
h.Position = [10,10,1700,900];
ha = axes;
hold(ha,'on');
axis equal
% beam splitter position
d_bs = d_m(1);
x_mbs = sum(d1(1:3))+d_bs; % x
y_mbs = 0; % y
p_mbs = mirror_angle(1); % angle
plot([x_mbs,x_mbs]+[-1,1]*20,[y_mbs,y_mbs]+[-1,1]*20*p_mbs,'k')
% get conversion matrix to BS coordinate
[A,Xm] = get_conversionM(x_mbs,y_mbs,p_mbs);

% mirror position relative to BS coordinates
x_ms_r = [1,1].*d_m(2)/sqrt(2);
y_ms_r = [-1,1].*d_m(2)/sqrt(2);
p_ms_r = [-1,1]*mirror_angle(2);
% convert to global coordinates
for ss=1:2
    [x_ms(ss),y_ms(ss),p_ms(ss)]  = getpos_conversion(A,Xm,x_ms_r(ss),y_ms_r(ss),p_ms_r(ss),'bs2global');
end

% find image plane at LFA
x_m = [x_mbs,x_ms(1)];
y_m = [y_mbs,y_ms(1)];
p_m = [p_mbs,p_ms(1)];
V0_image = calrayparam(fieldpoints(1),pupilpoints(1:2),d1(1));
[~,dx] = findimgplane(f1,d1,n,V0_image,aperturesize,ha,[],[],x_m,y_m,p_m,A,Xm);
d1(end) =  dx;

% mirror position, in the unfolded path
x_ms_r1 = [-1,-1].*d_m(2)/sqrt(2)+sum(d1(4:end))+f2(1)+d1(4);
y_ms_r1 = [-1,1].*d_m(2)/sqrt(2);
p_ms_r1 = -p_ms_r;
% convert to global coordinates
for ss=1:2
    [x_ms1(ss),y_ms1(ss),p_ms1(ss)]  = getpos_conversion(A,Xm,x_ms_r1(ss),y_ms_r1(ss),p_ms_r1(ss),'bs2global');
end

% BS position, in the unfolded path
x_mbs_r1 = sum(d1(4:end))+f2(1)+d1(4);
y_mbs_r1 = 0;
p_mbs_r1 = 0;
% convert to global coordinates
[x_mbs1,y_mbs1,p_mbs1]  = getpos_conversion(A,Xm,x_mbs_r1,y_mbs_r1,p_mbs_r1,'bs2global');
plot([x_mbs1,x_mbs1]+[-1,1]*20,[y_mbs1,y_mbs1]+[-1,1]*20*p_mbs1,'k')

% set up distances and focal lengths for ray tracing
[f2] = caldis(f2);
d=[d1(1:end-1),d1(end)+d2(1),d2(2:end)];
f=[f1(1:end-1),f2];
NL = numel(d);
n = [nmed,ones(1,NL)];
y_pupil = -1000.*ones(1,NL);
y_lens = -1000.*ones(1,NL);
dia_lens = [];

% calculate image size from path reflection path
x_m = [x_mbs,x_ms(1),x_ms1(1),x_mbs1];
y_m = [y_mbs,y_ms(1),y_ms1(1),y_mbs1];
p_m = [p_mbs,p_ms(1),p_ms1(1),p_mbs1];
V0_image = calrayparam(fieldpoints(1),pupilpoints(1:2),d1(1));
[dia_image1,df] = findimgplane(f,d,n,V0_image,aperturesize,ha,[],[],x_m,y_m,p_m,A,Xm);
V0_image = calrayparam(fieldpoints(3),pupilpoints(1:2),d1(1));
[dia_image2,~,pos_image] = findimgplane(f,d,n,V0_image,aperturesize,ha,'image plane','m--',x_m,y_m,p_m,A,Xm);
dia_image = (dia_image2-dia_image1)/2;
d(end) =  df;
for nn = 1:size(dia_image,1)
    text(ha,pos_image(nn,1),pos_image(nn,2),['image size: ',num2str(dia_image(nn),2),' mm'])
end

% start ray tracing
for ii = 1:3 % loop through field points
    for jj = 1:numel(pupilpoints) % loop through pupil points
        y0 = fieldpoints(ii);
        u0 = pupilpoints(jj);
        p0 = u0/d(1);
        x0 = 0;
        for kk = 1:3 % from object to a position before BS
            T = [1,0;d(kk)/n(kk),1]; % translation matrix
            O = [1,-n(kk)/f(kk);0,1]; % thin lens matrix
            r = O*T*[n(kk)*p0;y0];
            p1 = r(1)/n(kk+1);
            y1 = r(2);
            x1 = x0+d(kk);
            y_pupil(kk) = max([y_pupil(kk),p1*f(kk)+y1]);
            y_lens(kk) = max([y_lens(kk),y1]);
            plot([x0,x1],[y0,y1],cm{ii})
            if jj == numel(pupilpoints) && ii == numel(fieldpoints) && kk<NL
                plot([x1,x1],[-aperturesize/2,aperturesize/2],'k-')
                if labelflag == 1
                    text(x1,aperturesize,['f=',num2str(f(kk)),' mm'])
                    text(x1,aperturesize+20,['beam size: ',num2str(y_lens(kk)*2,3),' mm'])
                end
                dia_lens = cat(1,dia_lens,y_lens(kk)*2);
            end
            x0 = x1;
            y0 = y1;
            p0 = p1;
        end
        % trace to BS
        [x1,y1,p1]=getray_mirror(x_mbs,y_mbs,p_mbs,x0,y0,p0);
        plot([x0,x1],[y0,y1],cm{ii})
        % ray parameters at the BS, [reflection, transmission]
        x0_bs = [x1,x1];
        y0_bs = [y1,y1];
        p0_bs = [p1,p0];
        % line style for reflection and transmission
        linestyle = {'--','-'};
        for ss = 1:2 % loop over transmission and reflection paths
            x_m = x_ms(ss);
            y_m = y_ms(ss);
            p_m = p_ms(ss);
            x0 = x0_bs(ss);
            y0 = y0_bs(ss);
            p0 = p0_bs(ss);
            plot([x_m,x_m]+[-1,1]*20,[y_m,y_m]+[-1,1]*20*p_m,'k')
            % trace to mirror after the BS
            [x1,y1,p1]=getray_mirror(x_m,y_m,p_m,x0,y0,p0);
            plot([x0,x1],[y0,y1],cm{ii},'linestyle',linestyle{ss})
            for kk=4:5 % through remote focusing lens (two times)
                D = sum(d(4:kk));
                [x2,y2,p2,x2m,y2m] = getray_conversion(A,Xm,D,x1,y1,p1,f(kk));
                plot([x1,x2],[y1,y2],cm{ii},'linestyle',linestyle{ss})
                y_lens(kk) = max([y_lens(kk),y2m]);
                if jj == numel(pupilpoints) && ii == numel(fieldpoints) && kk<NL && ss==2 % draw lens and calculate beam size at the lens
                    X1_aperture = A\[x2m;aperturesize]+Xm;
                    X2_aperture = A\[x2m;-aperturesize]+Xm;
                    plot([X1_aperture(1),X2_aperture(1)],[X1_aperture(2),X2_aperture(2)],'k-','parent',ha)
                    if labelflag == 1
                        text(X1_aperture(1)+10,X1_aperture(2)-(kk-4.5)*30,['f=',num2str(f(kk)),' mm'])
                        text(X1_aperture(1)+10,X1_aperture(2)+20-(kk-4.5)*30,['beam size: ',num2str(y_lens(kk)*2,3),' mm'])
                    end
                    dia_lens = cat(1,dia_lens,y_lens(kk)*2);
                end
                x1 = x2;
                y1 = y2;
                p1 = p2;
            end
            x_m = x_ms1(ss);
            y_m = y_ms1(ss);
            p_m = p_ms1(ss);
            x0 = x1;
            y0 = y1;
            p0 = p1;
            plot([x_m,x_m]+[-1,1]*20,[y_m,y_m]+[-1,1]*20*p_m,'k')
            % trace to mirror in the unfold path
            [x1,y1,p1]=getray_mirror(x_m,y_m,p_m,x0,y0,p0);
            plot([x0,x1],[y0,y1],cm{ii},'linestyle',linestyle{ss})
            % trace to BS in the unfold path
            [x2,y2,p2]=getray_mirror(x_mbs1,y_mbs1,p_mbs1,x1,y1,p1);
            plot([x1,x2],[y1,y2],cm{ii},'linestyle',linestyle{ss})
            x1 = x2;
            y1 = y2;
            if ss==2
                p1 = p2;
            end
            for kk=6:NL % from BS to camera
                D = sum(d(6:kk));
                A1 = [1,0;0,1];
                Xm1 = [x_mbs1;y_mbs1];
                % convert to BS coordinate in the unfold path
                [x2,y2,p2,x2m,y2m] = getray_conversion(A1,Xm1,D,x1,y1,p1,f(kk));
                plot([x1,x2],[y1,y2],cm{ii},'linestyle',linestyle{ss})
                y_lens(kk) = max([y_lens(kk),y2m]);
                if jj == numel(pupilpoints) && ii == numel(fieldpoints) && kk<NL&& ss==2
                    X1_aperture = A1\[x2m;aperturesize]+Xm1;
                    X2_aperture = A1\[x2m;-aperturesize]+Xm1;
                    plot([X1_aperture(1),X2_aperture(1)],[X1_aperture(2),X2_aperture(2)],'k-','parent',ha)
                    if labelflag == 1
                        text(X1_aperture(1)+10,X1_aperture(2),['f=',num2str(f(kk)),' mm'])
                        text(X1_aperture(1)+10,X1_aperture(2)+20,['beam size: ',num2str(y_lens(kk)*2,3),' mm'])
                    end
                    dia_lens = cat(1,dia_lens,y_lens(kk)*2);
                end
                x1 = x2;
                y1 = y2;
                p1 = p2;
            end
        end
    end
end
magnification = abs(dia_image(end))/fov;
pixelsize = ccdpixelsize/magnification*1e3;
ha.Title.String = ['Ray tracing diagram, total magnification: ',num2str(magnification,3),', pixel size: ',num2str(pixelsize,4),' nm'];
ha.XLabel.String = 'Z (mm)';
ha.YLabel.String = 'Y (mm)';
end

%% utility functions
% function for finding the image plane, by tracing two rays from two different pupil points
function [dia_image,dx,pos_image] = findimgplane(f,d,n,V0,aperturesize,ha,planename,linestyle,x_m,y_m,p_m,A,Xm)
NL = numel(f);
dia_image = [];
pos_image = [];
x0 = 0;
for kk = 1:3 % from object to a position before BS
    T = [1,0;d(kk)/n(kk),1]; % translation matrix
    O = [1,-n(kk)/f(kk);0,1]; % thin lens matrix
    r = O*T*[n(kk).*V0(1,:);V0(2,:)];
    p1 = r(1,:)./n(kk+1);
    y1 = r(2,:);
    x1 = x0+d(kk); % at the lens surface
    if abs(diff(p1))>1e-6
        dx = -diff(y1)/diff(p1); % distance from image plane to the lens
        if (dx>0) && (dx<d(kk+1)+0.1)
            x_image = x1+dx;
            y_image = p1(1)*dx+y1(1);
            if ~isempty(linestyle)
                text(ha,x_image,-aperturesize,planename)
                plot([x_image,x_image],[-aperturesize/2,aperturesize/2],linestyle,'parent',ha)
            end
            dia_image = cat(1,dia_image,y_image*2);
            pos_image = cat(1,pos_image,[x_image,-aperturesize-20]);
        end
    end
    x0 = x1;
    V0 = [p1;y1];

end
x0 = [x1,x1];
y0 = y1;
p0 = p1;
% trace from BS to the mirror
for ii=1:2 
    [x1(1),y1(1),p1(1)]=getray_mirror(x_m(ii),y_m(ii),p_m(ii),x0(1),y0(1),p0(1));
    [x1(2),y1(2),p1(2)]=getray_mirror(x_m(ii),y_m(ii),p_m(ii),x0(2),y0(2),p0(2));

    x0 = x1;
    y0 = y1;
    p0 = p1;
end
% trace through the remote focusing lens
for kk=4:min([NL-1,5]) 
    D = sum(d(4:kk));
    for ss = 1:2
        [x1(ss),y1(ss),p1(ss),x1m(ss),y1m(ss),p1m(ss)] = getray_conversion(A,Xm,D,x0(ss),y0(ss),p0(ss),f(kk));

    end
    if abs(diff(p1m))>1e-6
        dx = -diff(y1m)/diff(p1m);
        if (dx>0) && (dx<d(kk+1)+0.1)
            x_image = x1m(1)+dx;
            y_image = p1m(1)*dx+y1m(1);
            X1_aperture = A\[x_image;aperturesize]+Xm;
            X2_aperture = A\[x_image;-aperturesize]+Xm;
            if ~isempty(linestyle)
                text(ha,X2_aperture(1)-aperturesize-50,X1_aperture(2)-120,planename)
                plot([X1_aperture(1),X2_aperture(1)],[X1_aperture(2),X2_aperture(2)],linestyle,'parent',ha)
            end
            dia_image = cat(1,dia_image,y_image*2);
            pos_image = cat(1,pos_image,[X2_aperture(1)-aperturesize-50,X1_aperture(2)-100]);
        end
    end
    x0=x1;
    y0 = y1;
    p0 = p1;
end

if NL>5
    % trace from the mirror to the BS in the unfolded path
    for ii=3:4 
        [x1(1),y1(1),p1(1)]=getray_mirror(x_m(ii),y_m(ii),p_m(ii),x0(1),y0(1),p0(1));
        [x1(2),y1(2),p1(2)]=getray_mirror(x_m(ii),y_m(ii),p_m(ii),x0(2),y0(2),p0(2));

        x0 = x1;
        y0 = y1;
        if ii==3
            p0 = p1;
        end
    end
    %trace from BS to the camera
    for kk=6:NL-1 
        A1 = [1,0;0,1];
        Xm1 = [x_m(4);y_m(4)];
        D = sum(d(6:kk));
        for ss = 1:2
            [x1(ss),y1(ss),p1(ss),x1m(ss),y1m(ss),p1m(ss)] = getray_conversion(A1,Xm1,D,x0(ss),y0(ss),p0(ss),f(kk));
        end
        if abs(diff(p1m))>1e-6
            dx = -diff(y1m)/diff(p1m);
            if (dx>0)
                x_image = x1m(1)+dx;
                y_image = p1m(1)*dx+y1m(1);
                X1_aperture = A1\[x_image;aperturesize]+Xm1;
                X2_aperture = A1\[x_image;-aperturesize]+Xm1;
                if ~isempty(linestyle)
                    text(ha,X1_aperture(1)+10,X1_aperture(2)-20,planename)
                    plot([X1_aperture(1),X2_aperture(1)],[X1_aperture(2),X2_aperture(2)],linestyle,'parent',ha)
                end
                dia_image = cat(1,dia_image,y_image*2);
                pos_image = cat(1,pos_image,[X1_aperture(1)+10,X1_aperture(2)]);
            end
        end
        x0=x1;
        y0 = y1;
        p0 = p1;
    end
end
end

function [f,d] = caldis(f)
f = [f,0];
NL = numel(f);

d = zeros(1,NL);
d(1) = f(1);
for nn = 1:NL-1
    d(nn+1) = f(nn)+f(nn+1);
end
end

function V0 = calrayparam(fieldpoints,pupilpoints,d0)
V0 = [];
for ii = 1:numel(fieldpoints)
    for jj = 1:numel(pupilpoints)
        y0 = fieldpoints(ii);
        u0 = pupilpoints(jj);
        p0 = u0/d0;
        V0 = cat(2,V0,[p0;y0]);
    end
end
end

function [x1,y1,p1]=getray_mirror(x_m,y_m,p_m,x0,y0,p0)

Xm = [x_m;y_m];
X0 = [x0;y0];
Pm_v = [-p_m,1]./norm([-p_m,1]);
P0_v = [-p0,1]./norm([-p0,1]);
Pm = [1,p_m]./norm([1,p_m]);
P0 = [1,p0]./norm([1,p0]);

A = [P0_v;Pm_v];
X1 = A\[P0_v*X0;Pm_v*Xm];

a = P0*Pm_v';
P1 = P0-2*a.*Pm_v;

x1 = X1(1);
y1 = X1(2);
p1 = P1(2)/P1(1);


end

function [x2,y2,p2,x2m,y2m,p2m] = getray_conversion(A,Xm,D,x1,y1,p1,f)
% convert to BS coordinates
[x1m,y1m,p1m]  = getpos_conversion(A,Xm,x1,y1,p1,'global2bs');

% translate by distance D and through lens f
T = [1,0;D-x1m,1];
if nargin>6
    O = [1,-1/f;0,1];
    T = O*T;
end
r = T*[p1m;y1m];
p2m = r(1);
y2m = r(2);
x2m = x1m+T(2,1);

% convert global coordinates
[x2,y2,p2]  = getpos_conversion(A,Xm,x2m,y2m,p2m,'bs2global');
end

function [A,X] = get_conversionM(x0,y0,p0)
X = [x0; y0];
vm = [1,p0];
vm1 = [1,-1/p0];

A = [vm./norm(vm);vm1./norm(vm1)];
end

function [x1m,y1m,p1m]  = getpos_conversion(A,Xm,x1,y1,p1,direction)
switch direction
    case 'global2bs'
        X1m = A*([x1;y1]-Xm);
        P1m = A*[1;p1];
        x1m = X1m(1);
        y1m = X1m(2);
        p1m = P1m(2)/P1m(1);
    case 'bs2global'
        X2 = A\[x1;y1]+Xm;
        P2 = A\[1;p1];
        x1m = X2(1);
        y1m = X2(2);
        p1m = P2(2)/P2(1);
end
end






