clear all;
close all;

addpath('mfile')

[v, f] = readOBJ('davidhead.obj');

plot_mesh(v,f); view([-130 0])

map = spherical_conformal_map(v,f);

writeOBJ('davidhead_embedding.obj', map, f)

map_bad = distort_sphere_torsion(f, v, map, pi);  

plot_mesh(map,f); 

writeOBJ('davidhead_embedding_bad.obj', map_bad, f)

plot_mesh(map_bad,f); 



function [V, F] = readOBJ(filename)
    fid = fopen(filename);
    V = [];
    F = [];

    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, 'v ')
            V = [V; sscanf(line(3:end), '%f %f %f')'];
        elseif startsWith(line, 'f ')
            f = sscanf(line(3:end), '%d %d %d');
            F = [F; f'];
        end
    end

    fclose(fid);
end

function writeOBJ(filename, V, F)
    fid = fopen(filename, 'w');
    if fid == -1
        error('Failed to open file: %s', filename);
    end

    % Write vertices
    for i = 1:size(V, 1)
        fprintf(fid, 'v %.8f %.8f %.8f\n', V(i, 1), V(i, 2), V(i, 3));
    end

    % Write faces
    for i = 1:size(F, 1)
        fprintf(fid, 'f %d %d %d\n', F(i, 1), F(i, 2), F(i, 3));
    end

    fclose(fid);
end


function map_bad = distort_radial_sphere_map(f, v, map, epsilon, k)
    
    x = map(:,1);
    y = map(:,2);
    z = map(:,3);

    theta = acos(z);
    phi = atan2(y, x); % azimuthal angle (optional)
    
    r = 1 + epsilon * sin(k * theta);

    map_bad = map .* r;

    map_bad = map_bad ./ vecnorm(map_bad, 2, 2);
end

function map_bad = distort_sphere_torsion(f, v, map, twist_strength)
    % Extract coordinates
    x = map(:,1);
    y = map(:,2);
    z = map(:,3);  % z ∈ [-1,1] on unit sphere

    % Compute twist angle per vertex: smoothly vary with z
    twist_angle = twist_strength * z;

    % Apply rotation around Z-axis
    cos_a = cos(twist_angle);
    sin_a = sin(twist_angle);
    
    x_new = cos_a .* x - sin_a .* y;
    y_new = sin_a .* x + cos_a .* y;
    z_new = z;

    map_bad = [x_new, y_new, z_new];

    % Normalize to remain on unit sphere
    map_bad = map_bad ./ vecnorm(map_bad, 2, 2);
end

