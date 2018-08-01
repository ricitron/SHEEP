clear
close all
% this is a script to calculate the volume of the chemical pile structures
% from the ASPECT calculations. This probably only works on mac or linux
% computers due to the way that paths are specified. It may require a
% recent version of matlab, since older versions of the 'dir' command do
% not return a 'folder' attribute.
%% settings and options:
aspect_output_folder = '~/Downloads/h300_B08_HPEx5_layered'; % NO trailing /
options.regular_grid_ny = 290; % 
options.nskip = 20; % number of output files to skip (saves time and memory if not needed)
options.load_particles = false; %we're not doing anything with the particles at the moment
options.output_fields = {'basal_layer','T','density'}; % list of field names to extract from output files
composition_threshold = 0.5; % composition greater than this number is treated as 'pile' material.
%% load the output and convert to a regular grid
aspect_output = load_aspect_output(aspect_output_folder,options);

%% postprocess output
% at each output step, calculate the fraction of the domain occupied by
% pile material
nt = length(aspect_output.time);
pile_fraction = zeros(nt,1);
for i=1:nt
    N = size(aspect_output.basal_layer,1)*size(aspect_output.basal_layer,2);
    pile_fraction(i) = sum(sum( aspect_output.basal_layer(:,:,i) >= 0.5 ))/N;
end

%% plot the fraction in the piles
figure;
plot(aspect_output.time,pile_fraction);
xlabel('Time (years)');
ylabel('Basal Layer Fraction (-)');

