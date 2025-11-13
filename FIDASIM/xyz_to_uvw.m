function xyz = xyz_to_uvw(alpha, beta, gamma, xyz, origin)
% Express non-rotated coordinate 'uvw' in rotated 'xyz' coordinates
% Arguments:
%     alpha: Rotation angle about z [radians]
%     beta: Rotation angle about y' [radians]
%     gamma: Rotation angle about x" [radians]
%     uvw: Point in rotated coordinate system
% Keyword Arguments:
%     origin: Origin of rotated coordinate system in non-rotated (uvw) coordinates.
%Function converted from D3D FIDASIM idl routines
if nargin < 5
    origin = [0.0, 0.0, 0.0];
end

s = size(xyz);
if numel(s) ~= 2
    s = [s, 1];
end


R = tb_zyx(alpha, beta, gamma);

xyz = R * xyz.';

xyz = xyz.'+repmat(origin, s(1), 1);
end