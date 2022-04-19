%spec_struct

clear; clc;
% Initialize fusome/ring system

x_voxel = 0.0699;
z_voxel = 0.2098;
voxelSize = [x_voxel x_voxel z_voxel];

%all samples should be 16 cells - get file for analysis
sample = '2_1';
num_cells = '2';
smooth = 5; %smoothing factor for fusome/ring splitting (can adjust)
uiopen(strcat('C:\Users\rockyd\Desktop\Manuscripts in Progress\Fusome_Evolution_CurrBiol\Spectrin_Pav_Samples\Structure\Trained\PavCherry_Spec546_',num_cells,'cell_',sample,'_Probabilities.h5'),1);

% preparing probabilities:
image_germ = double(I.img);
image_germ = permute(image_germ, [1 2 4 3]);
image_germ = image_germ/max(image_germ(:));
ring = image_germ(:,:,:,2);
fusome = image_germ(:,:,:,1);

% Clean up rings and segment fusome
bw_fusome = fusome > 0.5 | ring > 0.5;
frame = true(size(bw_fusome));
frame(2:end-1, 2:end-1, :) = false;
labels = bwlabeln(bw_fusome);
to_remove = integer_unique(labels(frame));
labels(ismember(labels, to_remove)) = 0;
labels(~isolate_lcc(labels)) = 0;
newfus = isolate_lcc(labels > 0 & fusome > 0.5);
newring = bwlabeln(labels > 0 & ring > 0.5);

%remove small objects that are not rings
s = bwskel(newring > 0);
change_occur = true;
while change_occur
    new_s = s;
    new_s(bwmorph3(s,'endpoints')) = false;
    change_occur = ~isequal(new_s, s);
    s = new_s;
end

% finish removing things that are not rings...
s = bwmorph3(s, 'clean');
rings_to_keep = integer_unique(nonzeros(newring(s)));
newring(~ismember(newring, rings_to_keep)) = 0;
newring = bwlabeln(newring > 0);
num_rings = length(unique(newring))-1;
fus_parts = length(unique(newfus))-1;
[x,y,z] = ind2sub(size(s), find(s));
x = x * voxelSize(1);
y = y * voxelSize(2);
z = z * voxelSize(3);

%expand out size to make isotropic directions
in_x = (0:size(newfus,1)-1)*voxelSize(1);
in_y = (0:size(newfus,2)-1)*voxelSize(2);
in_z = (0:size(newfus,3)-1)*voxelSize(3);

out_x = (0:size(newfus,1)-1)*voxelSize(1);
out_y = (0:size(newfus,2)-1)*voxelSize(1);
out_z = (0:size(newfus,3)-1)*voxelSize(1);

interp_newfus = medfilt3(interp3(single(newfus), 1:size(newfus,2), (1:size(newfus,1))', 1 : voxelSize(1)/voxelSize(3) : size(newfus,3), 'linear') >= 0.5);
interp_newring = medfilt3(bwlabeln(interp3(single(newring > 0), 1:size(newring,2), (1:size(newring,1))', 1 : voxelSize(1)/voxelSize(3) : size(newfus,3), 'linear') >= 0.5));
%% visualization before segmenting

for i = 1:num_rings
    p(i) = isosurface(interp_newring == i);
    p(i).vertices(:,1:2) = p(i).vertices(:,1:2)*voxelSize(1);
    p(i).vertices(:,3) = p(i).vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
end

figure;
for j = 1:length(p)
    patch(p(j),'FaceColor','r','EdgeColor','none'); axis image; alpha(1);
end

p1 = isosurface(interp_newfus == 1);
p1.vertices(:,1:2) = p1.vertices(:,1:2)*voxelSize(1);
p1.vertices(:,3) = p1.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p1,'FaceColor',[0 0.6 0.6],'EdgeColor','none'); axis image; alpha(1)
box off; grid off; axis off;
view([8 52])
%% split fusome and rings into component parts
clc;
newfus_labeled = split_fusome(0, interp_newfus, interp_newring, smooth);
%%
%Find adjacencies between objects
bw = imdilate(interp_newring, create_ball(1));
newimage = newfus_labeled.*(bw > 0);

if exist('list','var')
    clear list
end

for i = 1:(num_rings)
    check_obj = bwconncomp(bw == i,26);
    props = regionprops(check_obj,'Area','BoundingBox','Centroid','PixelIdxList','PixelList');
    list{i} = props.PixelIdxList;
end

matchup = zeros(length(list),2);

for j = 1:length(list)
    output = newimage(list{j});
    
    number = zeros(1,num_rings);
    for k = 1:num_rings+1
        number(k) = sum(output == k);
    end
    
    if nnz(number) < 2
        matchup(j,1) = 0;
        matchup(j,2) = 0;
    elseif nnz(find(number == max(number))) == 2
        [indices] = find(number == max(number));
        matchup(j,1) = indices(1);
        matchup(j,2) = indices(2);
    elseif max(number) == max(number(number < max(number)))
        matchup(j,1) = find(number == max(number));
        matchup(j,2) = find(number == max(number));
    else
        matchup(j,1) = find(number == max(number));
        matchup(j,2) = find(number == max(number(number<max(number))));
    end
end

% Remove zero rows and rows with one element
rule = matchup(:,2) == 0;
matchup(rule,:) = [];

% Remove zero columns
matchup(:,all(~matchup,1)) = [];

matchupswitch = [matchup(:,2) matchup(:,1)];
matchup = unique([matchup; matchupswitch], 'rows');

%Remove duplicate pairs
pair = matchup(:,1) < matchup(:,2);
matchup = matchup.*pair;
rule = matchup(:,2) == 0;
matchup(rule,:) = [];

total_seg = max(matchup(:));

%get fusome volumes for each portion
fus_vol_out = zeros(total_seg,1);
for i = 1:total_seg
    fus_vol_out(i) = nnz(newfus_labeled(:) == i).*voxelSize(1).*voxelSize(2).*voxelSize(3).*voxelSize(1)./voxelSize(3);
end
oocyte = find(fus_vol_out == max(fus_vol_out)); %assume cell 1 has largest fusome piece

matchupmatrix = zeros(total_seg);
for i = 1:size(matchup,1)
    matchupmatrix(matchup(i,1),matchup(i,2)) = 1;
    matchupmatrix(matchup(i,2),matchup(i,1)) = 1;
end

%compare adjacency matrix with known 2,4,8,16-cell cyst adjancencies
if total_seg == 2
    truematchup = [1 2];
    truematchupmatrix = zeros(2);
elseif total_seg == 4
    truematchup = [1 2;1 3;2 4];
    truematchupmatrix = zeros(4);
elseif total_seg == 8
    truematchup = [1 2; 1 3; 1 5; 2 4; 2 6; 3 7;4 8];
    truematchupmatrix = zeros(8);
else
    truematchup = [1 2;1 3;1 5;1 9;2 4; 2 6;2 10;3 7;3 11;4 12;4 8;5 13;6 14;7 15;8 16];
    truematchupmatrix = zeros(16);
end

for i = 1:size(truematchup,1)
    truematchupmatrix(truematchup(i,1),truematchup(i,2)) = 1;
    truematchupmatrix(truematchup(i,2),truematchup(i,1)) = 1;
end

matchupmatrix = zeros(total_seg);
for i = 1:size(matchup,1)
    matchupmatrix(matchup(i,1),matchup(i,2)) = 1;
    matchupmatrix(matchup(i,2),matchup(i,1)) = 1;
end

matchupmatrix(oocyte,:) = 2*matchupmatrix(oocyte,:);
matchupmatrix(:,oocyte) = 2*matchupmatrix(:,oocyte);

truematchupmatrix(1,:) = 2*truematchupmatrix(1,:);
truematchupmatrix(:,1) = 2*truematchupmatrix(:,1);

%compare and get permutation to reassign identities for each object
[P1, ~] = eig(truematchupmatrix);
[P2, ~] = eig(matchupmatrix);

P1 = abs(P1);
P2 = abs(P2);

outputmatrix = P2*P1';
[permutation, ~] = munkres(outputmatrix);
permutation = ((total_seg+1) - permutation)';

if length(permutation) == 2 && permutation(1) == 2 && permutation(2) == 1
    permutation = [1;2];
end

% relabel fusome pieces
output_fus = zeros(size(newfus_labeled));
fusome_vol = zeros(size(fus_vol_out));
for m = 1:(num_rings+1)
    output_fus(newfus_labeled == m) = permutation(m);
    fusome_vol(permutation(m)) = fus_vol_out(m);
end
%%
% Create 3D reconstruction

for i = 1:num_rings
    p(i) = isosurface(interp_newring == i);
    p(i).vertices(:,1:2) = p(i).vertices(:,1:2)*voxelSize(1);
    p(i).vertices(:,3) = p(i).vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
end

figure;
for j = 1:length(p)
    patch(p(j),'FaceColor','k','EdgeColor','none'); axis image; alpha(1)
end

if max(unique(output_fus)) >= 2
    p1 = isosurface(output_fus == 1);
    p1.vertices(:,1:2) = p1.vertices(:,1:2)*voxelSize(1);
    p1.vertices(:,3) = p1.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p1,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none'); axis image; alpha(1)
    
    p2 = isosurface(output_fus == 2);
    p2.vertices(:,1:2) = p2.vertices(:,1:2)*voxelSize(1);
    p2.vertices(:,3) = p2.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p2,'FaceColor',[118/255 173/255 190/255],'EdgeColor','none'); axis image; alpha(1)
end

if max(unique(output_fus)) >= 4
    p3 = isosurface(output_fus == 3);
    p3.vertices(:,1:2) = p3.vertices(:,1:2)*voxelSize(1);
    p3.vertices(:,3) = p3.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p3,'FaceColor',[118/255 173/255 190/255],'EdgeColor','none'); axis image; alpha(1)
    
    p4 = isosurface(output_fus == 4);
    p4.vertices(:,1:2) = p4.vertices(:,1:2)*voxelSize(1);
    p4.vertices(:,3) = p4.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p4,'FaceColor',[192/255 97/255 90/255],'EdgeColor','none'); axis image; alpha(1)
end

if max(unique(output_fus)) >= 8
    p5 = isosurface(output_fus == 5);
    p5.vertices(:,1:2) = p5.vertices(:,1:2)*voxelSize(1);
    p5.vertices(:,3) = p5.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p5,'FaceColor',[118/255 173/255 190/255],'EdgeColor','none'); axis image; alpha(1)
    
    p6 = isosurface(output_fus == 6);
    p6.vertices(:,1:2) = p6.vertices(:,1:2)*voxelSize(1);
    p6.vertices(:,3) = p6.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p6,'FaceColor',[192/255 97/255 90/255],'EdgeColor','none'); axis image; alpha(1)
    
    p7 = isosurface(output_fus == 7);
    p7.vertices(:,1:2) = p7.vertices(:,1:2)*voxelSize(1);
    p7.vertices(:,3) = p7.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p7,'FaceColor',[192/255 97/255 90/255],'EdgeColor','none'); axis image; alpha(1)
    
    p8 = isosurface(output_fus == 8);
    p8.vertices(:,1:2) = p8.vertices(:,1:2)*voxelSize(1);
    p8.vertices(:,3) = p8.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p8,'FaceColor',[79/255 124/255 98/255],'EdgeColor','none'); axis image; alpha(1)
end

if max(unique(output_fus)) >= 16
    p9 = isosurface(output_fus == 9);
    p9.vertices(:,1:2) = p9.vertices(:,1:2)*voxelSize(1);
    p9.vertices(:,3) = p9.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p9,'FaceColor',[118/255 173/255 190/255],'EdgeColor','none'); axis image; alpha(1)
    
    p10 = isosurface(output_fus == 10);
    p10.vertices(:,1:2) = p10.vertices(:,1:2)*voxelSize(1);
    p10.vertices(:,3) = p10.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p10,'FaceColor',[192/255 97/255 90/255],'EdgeColor','none'); axis image; alpha(1)
    
    p11 = isosurface(output_fus == 11);
    p11.vertices(:,1:2) = p11.vertices(:,1:2)*voxelSize(1);
    p11.vertices(:,3) = p11.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p11,'FaceColor',[192/255 97/255 90/255],'EdgeColor','none'); axis image; alpha(1)
    
    p12 = isosurface(output_fus == 12);
    p12.vertices(:,1:2) = p12.vertices(:,1:2)*voxelSize(1);
    p12.vertices(:,3) = p12.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p12,'FaceColor',[79/255 124/255 98/255],'EdgeColor','none'); axis image; alpha(1)
    
    p13 = isosurface(output_fus == 13);
    p13.vertices(:,1:2) = p13.vertices(:,1:2)*voxelSize(1);
    p13.vertices(:,3) = p13.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p13,'FaceColor',[192/255 97/255 90/255],'EdgeColor','none'); axis image; alpha(1)
    
    p14 = isosurface(output_fus == 14);
    p14.vertices(:,1:2) = p14.vertices(:,1:2)*voxelSize(1);
    p14.vertices(:,3) = p14.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p14,'FaceColor',[79/255 124/255 98/255],'EdgeColor','none'); axis image; alpha(1)
    
    p15 = isosurface(output_fus == 15);
    p15.vertices(:,1:2) = p15.vertices(:,1:2)*voxelSize(1);
    p15.vertices(:,3) = p15.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p15,'FaceColor',[79/255 124/255 98/255],'EdgeColor','none'); axis image; alpha(1)
    
    p16 = isosurface(output_fus == 16);
    p16.vertices(:,1:2) = p16.vertices(:,1:2)*voxelSize(1);
    p16.vertices(:,3) = p16.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
    patch(p16,'FaceColor',[252/255 202/255 84/255],'EdgeColor','none'); axis image; alpha(1)
end
box on; grid on; axis off; axis image
view([8 52])