function [output] = load_aspect_output(output_directory,options)
% Read the .xdmf files
%
% Max Rudolph, August 2018

if( nargin < 2 )
    options.nskip=1;
    options.regular_grid_ny = 100;
    options.load_particles = false;
end
if( ismac || isunix )
    path_sep = '/';
else
    path_sep = '\';
end

% first load the 'particles' xdmf file
if( options.load_particles )
    file = dir([output_directory path_sep 'particles.xdmf']);
    
    % folder = [output_directory];
    folder = output_directory;
    copyfile('Xdmf.dtd',file.folder);
    
    % read the xdmf file and parse into a matlab structure:
    
    data_description = xml2struct([file.folder path_sep file.name]);
    
    xdmf = data_description.Xdmf;
    xdmf = xdmf{2};
    domain = xdmf.Domain;
    grid = domain.Grid.Grid;
    
    % parse the output times, output files, and corresponding grid files for
    % the particles
    nfiles = length(grid);
    grid_files = {};
    model_time = [];
    
    for ifile=1:nfiles
        output_time(ifile) = str2num(grid{ifile}.Time.Attributes.Value);
        tmp = strtrim(grid{ifile}.Topology.DataItem.Text);
        tmp1 = tmp(1:(find(tmp==':')-1));
        tmp2 = tmp(find(tmp==':')+1:end);
        output_grid_file{ifile} = [folder path_sep tmp1];
        output_grid_fieldname{ifile} = tmp2;
        output_toplogy{ifile} = strtrim(grid{ifile}.Topology.DataItem.Text);
        tmp = strtrim(grid{ifile}.Attribute{1}.DataItem.Text);
        tmp = tmp(1:find(tmp==':')-1);
        output_data_file{ifile} = [folder path_sep tmp];
    end
    
    % load first and last data files
    initial_id = h5read( output_data_file{1},'/id')';
    initial_xy = h5read( output_data_file{1},'/nodes')';
    initial_composition = h5read( output_data_file{1},'/initial basal_layer')';
    
    % sort particles by tracer id
    [~,i] = sort(initial_id);
    initial_id = initial_id(i);
    initial_xy = initial_xy(i,:);
    initial_composition = initial_composition(i);
    ny = length(unique(initial_xy(:,1)));
    nx = length(unique(initial_xy(:,2)));
    initial_x = reshape(initial_xy(:,1),[nx ny]);
    initial_y = reshape(initial_xy(:,2),[nx ny]);
    
    % sort the tracers by increasing y-coordinate and then by increasing
    % x-coordinates. This is necessary so that tracers are in the same
    % arrangement as would be produced from Paul Tackley's Convection2D_tracers
    % matlab program.
    [~,i2] = sortrows( initial_xy,[1 2]);
    % i2 tells us how to re-sort tracers from sequentially-increasing id to the
    % same order as the original xy grid with fastest increasing coordinate in
    % the y-direction, then the x-direction
    initial_id = initial_id(i2);
    initial_xy = initial_xy(i2,:);
    initial_composition = initial_composition(i2);
    initial_composition = reshape(initial_composition,[nx ny]);
    nt = length(1:options.nskip:length(output_data_file));
    Xn = zeros(nx,ny,nt);
    Yn = zeros(nx,ny,nt);
    C = zeros(nx,ny);
    Xn(:,:,1) = initial_x;
    Yn(:,:,1) = initial_y;
    C(:,:) = initial_composition;
    
    isave=2;
    for ifile=(1+options.nskip):options.nskip:length(output_data_file)
        disp(['loading particles from ' output_data_file{ifile}]);
        final_id = h5read( output_data_file{ifile},'/id')';
        final_xy = h5read( output_data_file{ifile},'/nodes')';
        % first sort based on the tracer id
        [~,i] = sort(final_id);
        final_id = final_id(i);
        final_xy = final_xy(i,:);
        % then re-order to match initial ordering
        final_xy = final_xy(i2,:);
        final_id = final_id(i2);
        % this assertion ensures that the initial and final tracers are all
        % accounted for. Now the tracers should be ordered by id
        assert(all(initial_id==final_id))
        
        final_x = reshape(final_xy(:,1),[nx,ny]);
        final_y = reshape(final_xy(:,2),[nx,ny]);
        
        Xn(:,:,isave) = final_x;
        Yn(:,:,isave) = final_y;
        
        isave=isave+1;
    end
    t=output_time(1:options.nskip:length(output_data_file));
end

% Now parse the xdmf for the 'solution'
file = dir([output_directory path_sep 'solution.xdmf']);
folder = [output_directory];
copyfile('Xdmf.dtd',file.folder);
% read the xdmf file and parse into a matlab structure:
data_description = xml2struct([file.folder path_sep file.name]);

xdmf = data_description.Xdmf;
xdmf = xdmf{2};
domain = xdmf.Domain;
grid = domain.Grid.Grid;

% parse the output times, output files, and corresponding grid files for
% the particles
nfiles = length(grid);
grid_files = {};
model_time = [];
clear output_grid_file output_grid_fieldname output_topology output_data_file
for ifile=1:nfiles
    output_time(ifile) = str2num(grid{ifile}.Time.Attributes.Value);
    output_time(ifile) = str2num(grid{ifile}.Time.Attributes.Value);
    tmp = strtrim(grid{ifile}.Topology.DataItem.Text);
    tmp1 = tmp(1:(find(tmp==':')-1));
    tmp2 = tmp(find(tmp==':')+1:end);
    output_grid_file{ifile} = [folder path_sep tmp1];
    output_grid_fieldname{ifile} = tmp2;
    output_toplogy{ifile} = strtrim(grid{ifile}.Topology.DataItem.Text);
    tmp = strtrim(grid{ifile}.Attribute{1}.DataItem.Text);
    tmp = tmp(1:find(tmp==':')-1);
    output_data_file{ifile} = [folder path_sep tmp];
end
t=output_time(1:options.nskip:length(output_data_file));
nt = length(t);
% read the nodal fields from the h5 files. the mesh changes between steps
% so we need to save both the nodal coordinates and the values.
nodes = h5read(output_grid_file{1},'/nodes')';
[~,i] = sortrows( nodes,[1 2]);
X = nodes(i,1);
% Xu = unique(X);
% Y = nodes(i,2);
% Yu = unique(Y);
% remove duplicates from list of nodes:
[~,inodeu] = unique(nodes,'rows');

% generate a regular grid for output
xmax = max(nodes(:,1));
ymax = max(nodes(:,2));
ny = options.regular_grid_ny;
nx = ny/ymax*xmax;
xrg = linspace(0,xmax,nx+1); xrg = (xrg(1:end-1)+xrg(2:end))/2;
yrg = linspace(0,ymax,ny+1); yrg = (yrg(1:end-1)+yrg(2:end))/2;
[xrg,yrg] = meshgrid(xrg,yrg);
comp_rg = zeros(ny,nx,nt);

isave=1;
for ifile=1:options.nskip:length(output_data_file)
    disp(['reading grid output from ' output_data_file{ifile}]);
    velocity = h5read(output_data_file{ifile},'/velocity');
    temperature = h5read(output_data_file{ifile},'/T');
    nodes = h5read(output_grid_file{ifile},'/nodes')';
    basal_layer = h5read(output_data_file{ifile},'/basal_layer')';
    % interpolate the nodal composition values from the structured,
    % irregular ASPECT mesh onto the regular grid:
    F = scatteredInterpolant(nodes(:,1),nodes(:,2),basal_layer);
    comp_rg(:,:,isave) = F(xrg,yrg);
    
    isave=isave+1;
end

output.composition = comp_rg;
output.x = xrg;
output.y = yrg;
output.time = t;
output.xmax = xmax;
output.ymax = ymax;

return