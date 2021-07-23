function [ output_args ] = volshow(varargin)
a=0;
if  length(varargin)==1
    target_volume = varargin{1};
    axes_aspect_ratio = [1.0, 1.0, 1];
    alpha_val = 0.5
%     color_code = 'r'
    color_code = rand(1,3);
elseif length(varargin)==2
    target_volume = varargin{1};
    axes_aspect_ratio = varargin{2};
    alpha_val = 0.5
    color_code = 'r'
    
elseif length(varargin)==3
    target_volume = varargin{1};
    axes_aspect_ratio = varargin{2};
    alpha_val = varargin{3};
    color_code = 'r'
    
elseif length(varargin)==4
    target_volume = varargin{1};
    axes_aspect_ratio = varargin{2};
    alpha_val = varargin{3};
    color_code = varargin{4};
end



% volume_, axes_ratio_, alpha_
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

p = patch( isosurface(target_volume ,0) );                 %# create isosurface patch
isonormals(target_volume , p)                              %# compute and set normals
set(p, 'FaceColor',color_code , 'EdgeColor','none')   %# set surface props
daspect(axes_aspect_ratio)                              %# axes aspect ratio
view(3), axis vis3d tight, box on, grid on    %# set axes props
camproj perspective                           %# use perspective projection
camlight, lighting phong, alpha(alpha_val)           %# enable light, set transparency


end

