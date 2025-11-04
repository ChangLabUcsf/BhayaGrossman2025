function [electrodes, a]=ctmr_gauss_plot(cortex,electrodes,weights,hemi)
% function [electrodes]=ctmr_gauss_plot(cortex,electrodes,weights)
% projects electrode locations onto their cortical spots in the 
% left hemisphere and plots about them using a gaussian kernel
% for only cortex use: 
% ctmr_gauss_plot(cortex,[0 0 0],0)
% rel_dir=which('loc_plot');
% rel_dir((length(rel_dir)-10):length(rel_dir))=[];
% addpath(rel_dir)

%     Copyright (C) 2009  K.J. Miller & D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%   Version 1.1.0, released 26-11-2009

load('loc_colormap_thresh');
brain=cortex.vert;
v='l';

if length(weights)~=length(electrodes(:,1))
    error('you sent a different number of weights than electrodes (perhaps a whole matrix instead of vector)')
end
gsp=10; %zg edited from 50

c=zeros(length(cortex(:,1)),1);
for i=1:length(electrodes(:,1))
    b_z=abs(brain(:,3)-electrodes(i,3));
    b_y=abs(brain(:,2)-electrodes(i,2));
    b_x=abs(brain(:,1)-electrodes(i,1));
    d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian 
    c=c+d';
end

% c=(c/max(c));
a=tripatch(cortex, 'nofigure', c');
shading interp;
a=get(gca);
%%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
%%RELEVANCE OF CHANNEL FOR COMPARISON's ACROSS CORTICES
d=a.CLim;
set(gca,'CLim',[-max(abs(d)) max(abs(d))])
l=light;
colormap(cm);
a.FaceColor = [0.9 0.9 0.9];
lighting phong; %play with lighting...
material shiny;
% material dull;
alpha 0.9;
material([.3 .8 .1 10 1]);
a.HandleVisility = 'off';
axis off
if strcmp(hemi, 'lh')
    view(270, 0);
    set(l,'Position',[-1 0 0],'Color',[0.8 0.8 0.8]);
elseif strcmp(hemi, 'rh')
    view(90, 0);
    set(l,'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
end
