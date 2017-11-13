clear
load('final1__Rmatrix'); % load the response matrix
% 
% n = [0 1  1  2  2 2  3 3  3 3 4 4  4 4  4];
% m = [0 1 -1  0 -2 2 -1 1 -3 3 0 2 -2 4 -4];
n = [0  1  1  2  2 2  3 3  3 3 4 4  4 4  4  5 5  5 5  5 5 6  6 6  6 6  6 6];
m = [0  1 -1  0 -2 2 -1 1 -3 3 0 2 -2 4 -4 -1 1 -3 3 -5 5 0 -2 2 -4 4 -6 6];
radius = 44;
% normalize the axes
num_points = 2*(radius+5)+1;%size(DOP,1);
[x, y]=meshgrid(linspace((-radius-5)/radius,(radius+5)/radius, num_points),linspace((-radius-5)/radius,(radius+5)/radius, num_points));
[qi,ri] = cart2pol(x,y);
IOI = ri<=1;
Z = zernfun(n,m,ri(IOI),qi(IOI));
wavefront_0 = zeros(num_points,num_points,10);
DOP=zeros(size(num_points,1));

for i=2:28
    %28 or 15
    load(strcat('final1__wavefront_Z',num2str(i),'.mat'));
    wavefront_0(:,:,i-1)=squeeze(wavefront(5,:,:));
end  
base = (mean(wavefront_0,3));%fliplr
xxxx=base(IOI);
RMS = std(xxxx)   
diff_I = DOP-base;  % get the difference between reference and the changed wavefront
diff_I=squeeze(diff_I);
diff_I=diff_I-mean(diff_I(IOI));
% diff_I=imgaussfilt(diff_I,0.5);

figure
imagesc(-diff_I)
title(strcat('RMS = ', num2str(RMS),' nm'))
set(gca,'YDir','normal')
colormap Jet
caxis([-100 100])
colorbar

% decompose the wavefront
a = Z\diff_I(IOI);

% get the changed Zernike term
Change_Zernike = (Rmatrix(:,1)).*(a-(Rmatrix(:,2)));

correction = zeros(num_points,num_points);
for i = 1:28
    % 11 or 15
    correction(IOI)=correction(IOI)+ Change_Zernike(i).*Z(:,i);
end
a = Change_Zernike'
figure
imagesc(correction)
title('correction')
set(gca,'YDir','normal')
colormap Jet
caxis([-100 100])
colorbar