function xyz = uvw_to_xyz(alpha, beta, gamma, uvw, origin)
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

s = size(uvw);
if numel(s) ~= 2
    s = [s, 1];
end

uvw_shifted = uvw - repmat(origin, s(1), 1);

R = tb_zyx(alpha, beta, gamma).';

xyz = R * uvw_shifted.';

xyz = xyz.';
end