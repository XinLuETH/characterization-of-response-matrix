clear
%set geometry
X = 775;        %6X6 779   9X9 779    11X11: 777
Y = 667;        %6X6 669 9X9 668    11X11: 669
center = [X,Y];
radius = 44;   %6X6 22 9X9 44     11X11: 54

max_zernike = 28;  % 11 or 15  
% image size
num_points = 2*(radius+5)+1;
[x, y]=meshgrid(linspace((-radius-5)/radius,(radius+5)/radius, num_points),linspace((-radius-5)/radius,(radius+5)/radius, num_points));
[qi,ri] = cart2pol(x,y);
IOI = ri<=1;
% spatial frequency filter size
filter_size = 8;
[F_x, F_y]=meshgrid(linspace((-radius-5)/filter_size,(radius+5)/filter_size, num_points),linspace((-radius-5)/filter_size,(radius+5)/filter_size, num_points));
[F_qi,F_ri] = cart2pol(F_x,F_y);
F_IOI = F_ri<=1;


%%load gain map of the camera
% load('gain'); % load the response matrix
% Gain = gain(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5);
%set x axis for plot 
x_plot = [-20 -15 -10 -5 0 5 10 15 20];

% get the first 15 zernike terms
% n = [0 1  1  2  2 2  3 3  3 3 4];
% m = [0 1 -1  0 -2 2 -1 1 -3 3 0];
% n = [0  1  1  2  2 2  3 3  3 3 4 4  4 4  4];
% % m = [0  1 -1  0 -2 2 -1 1 -3 3 0 2 -2 4 -4];
% n = [0  1  1  2  2 2  3 3  3 3 4 4  4 4  4  5 5  5 5  5 5 6  6 6  6 6  6 6  7 7  7 7  7 7  7 7 8  8 8  8 8  8 8  8 8  9 9 9  9  9 9  9 9  9 9 10 10 10 10 10 10 10 10 10  10 10 11 11 11 11 11 11 11 11 11 11  11 11];
% m = [0  1 -1  0 -2 2 -1 1 -3 3 0 2 -2 4 -4 -1 1 -3 3 -5 5 0 -2 2 -4 4 -6 6 -1 1 -3 3 -5 5 -7 7 0 -2 2 -4 4 -6 6 -8 8 -1 1 -3 3 -5 5 -7 7 -9 9  0 -2  2 -4  4 -6  6 -8  8 -10 10 -1  1 -3  3 -5  5 -7  7 -9  9 -11 11];
n = [0  1  1  2  2 2  3 3  3 3 4 4  4 4  4  5 5  5 5  5 5 6  6 6  6 6  6 6];
m = [0  1 -1  0 -2 2 -1 1 -3 3 0 2 -2 4 -4 -1 1 -3 3 -5 5 0 -2 2 -4 4 -6 6];
% n = [0  1  1  2  2  2  3  3  3  3  4  4 4 4 4];
% m = [0 -1  1 -2  0  2 -3 -1  1  3 -4 -2 0 2 4];

% read all the fits files
% DC = fitsread('Bdark_50ms.fits');   % read the dark current
% DC = DC(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5, :);
% set a tring to read fits
header1 = 'final1_';
header = {strcat(header1,'N20_Z'), strcat(header1,'N15_Z'),strcat(header1,'N10_Z'),strcat(header1,'N5_Z'),strcat(header1,'0_Z'),strcat(header1,'5_Z'),strcat(header1,'10_Z'),strcat(header1,'15_Z'),strcat(header1,'20_Z')};

% set the wavefront matrix
wavefront = zeros(9,2*(radius+5)+1,2*(radius+5)+1);
% set the response matrix
Rmatrix = zeros(28,2);

for j = 2:max_zernike
    
    for i = 1:9
        data0 = fitsread(strcat(header{i}, num2str(j),'_0Pi','.fits'));   % read data of 0*pi/2
        data1 = fitsread(strcat(header{i}, num2str(j),'_1Pi','.fits'));  % read data of 1*pi/2
        data2 = fitsread(strcat(header{i}, num2str(j),'_2Pi','.fits'));   % read data of 2*pi/2
        data3 = fitsread(strcat(header{i}, num2str(j),'_3Pi','.fits'));   % read data of 3*pi/2

        % reduce all the data size, exchange X and Y 
        DATA0 = data0(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5);
        DATA1 = data1(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5);
        DATA2 = data2(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5);
        DATA3 = data3(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5);
    
        % flip DATA along the x axis
        DATA0 = fliplr(DATA0);
        DATA1 = fliplr(DATA1);
        DATA2 = fliplr(DATA2);
        DATA3 = fliplr(DATA3);
%         DATA0 = fliplr(DATA0)./Gain;
%         DATA1 = fliplr(DATA1)./Gain;
%         DATA2 = fliplr(DATA2)./Gain;
%         DATA3 = fliplr(DATA3)./Gain;
        %% remove the patters 
        % set data outside the circle to zero
%         for k=1:size(DATA3,3)
%             f0 = fftshift(fft2(DATA0(:,:,k)));
%             f1 = fftshift(fft2(DATA1(:,:,k)));
%             f2 = fftshift(fft2(DATA2(:,:,k)));
%             f3 = fftshift(fft2(DATA3(:,:,k)));
%         % used to plot the image
%             fLog0 = log(1 + abs(f0));
%             fLog1 = log(1 + abs(f1));
%             fLog2 = log(1 + abs(f2));
%             fLog3 = log(1 + abs(f3));
%         % filter by a range based on fLog
%             filter0 = fLog0 < .7*max(fLog0(:));
%             filter1 = fLog1 < .7*max(fLog1(:));
%             filter2 = fLog2 < .7*max(fLog2(:));
%             filter3 = fLog3 < .7*max(fLog3(:));
%             filter0(54:78,54:78) = 1;
%             filter1(54:78,54:78) = 1;
%             filter2(54:78,54:78) = 1;
%             filter3(54:78,54:78) = 1;
%             DATA0(:,:,k) = abs(ifft2(f0.*filter0));
%             DATA1(:,:,k) = abs(ifft2(f1.*filter1));
%             DATA2(:,:,k) = abs(ifft2(f2.*filter2));
%             DATA3(:,:,k) = abs(ifft2(f3.*filter3));
%         end
        %% spatial filter, set data outside the circle to zero
        for k=1:size(DATA0,3)
            f0 = fftshift(fft2(DATA0(:,:,k)));
            f1 = fftshift(fft2(DATA1(:,:,k)));
            f2 = fftshift(fft2(DATA2(:,:,k)));
            f3 = fftshift(fft2(DATA3(:,:,k)));
        % filter by a range based on fLog
            f0(~F_IOI) = 0;
            f1(~F_IOI) = 0;
            f2(~F_IOI) = 0;
            f3(~F_IOI) = 0;
            
            DATA0(:,:,k) = real(ifft2(ifftshift(f0)));
            DATA1(:,:,k) = real(ifft2(ifftshift(f1)));
            DATA2(:,:,k) = real(ifft2(ifftshift(f2)));
            DATA3(:,:,k) = real(ifft2(ifftshift(f3)));
        end
%         
%           F = zeros(size(DATA0));
%           F(F_IOI)=1;
%           f0 = fftshift(fft2(DATA0));
%           f1 = fftshift(fft2(DATA1));
%           f2 = fftshift(fft2(DATA2));
%           f3 = fftshift(fft2(DATA3));
%           f0=f0.*F;  
%           f1=f1.*F; 
%           f2=f2.*F; 
%           f3=f3.*F; 
%           figure
%           imagesc(10*log(abs(f1.^2)))
                
% figure
% imagesc(DATA1(:,:,1))%
% colormap Jet
% set(gca,'YDir','normal')
% hold on



        % get the average of DATA0 in the ROI
%         [x, y]=meshgrid(1:size(DATA0,1),1:size(DATA0,2));
%         % the cernter of the circle changes after resize the data
%         XX = radius + 5 + 1;
%         YY = radius + 5 + 1;
%         % set data outside the circle to zero
%         Circle = (x-XX).^2+(y-YY).^2 <= radius^2;
%         DATA0(~Circle) = 0;
%         DATA1(~Circle) = 0;
%         DATA2(~Circle) = 0;
%         DATA3(~Circle) = 0;

        % Zernike sensing phase-shifting
        I0 = (DATA0+DATA1+DATA2+DATA3);
        phase_1 = (DATA1-DATA3)./I0;
%         phase_1 = (DATA0-DATA2)./I0;  amplitude
        phase_1(~IOI) = 0;
% %         
        % Zernike sensor
%         phase_1 = zeros(size(DATA0));
%         ave_DATA0 = mean(DATA0(IOI));
%         DATA1 = DATA1/ave_DATA0;
%         phase_1(IOI) = (-1+sqrt(3-2*0.5-(1-DATA1(IOI))/0.5));
%         phase_1(IOI) = DATA1(IOI)-0.5;
%         phase_1(~IOI) = 0;
% % 
% 
        DOP1 = phase_1*633/2/pi;
%         h = ones(5,5) / 25;
%         DOP = imfilter(DOP1,h);  change the the four images
%         DOP(~IOI) = 0;
%         DOP=imgaussfilt(DOP1,0.3);
%         DOP=imrotate(DOP,-3,'crop');
        wavefront(i,:,:) = DOP1;   
    end
%     figure
%     imagesc(squeeze(wavefront(5,:,:)))
%     set(gca,'YDir','normal')
%     colormap Gray
%     colorbar
    save(strcat(header1,'_wavefront_Z',num2str(j),'.mat'),'wavefront');
    % normalize the axes
    num_points = size(DOP1,1);
    [x, y]=meshgrid(linspace((-radius-5)/radius,(radius+5)/radius, num_points),linspace((-radius-5)/radius,(radius+5)/radius, num_points));
    [qi,ri] = cart2pol(x,y);
    IOI = ri<=1;
%     get the Zernike terms
    Z = zernfun(n,m,ri(IOI),qi(IOI));
    decompose_zernike = zeros(max_zernike,9);
    for l = 1:9
        diff_wavefront = wavefront(l,:,:)-wavefront(5,:,:);
        decompose_zernike(:,l) = Z\diff_wavefront(IOI);
    end
    save(strcat(header1,'_Decomposed_coef_of_Z',num2str(j),'.mat'),'decompose_zernike');
%     get the linear fitting
    P = polyfit(x_plot,decompose_zernike(j,:),1);
    yfit = P(1)*x_plot+P(2);
%     get thr r^2 of the linear fitting
    y_ave=mean(decompose_zernike(j,:));
    SSE = 0;
    SST = 0;
    for r=1:size(x_plot)
        SSE=SSE+(decompose_zernike(j,r)-yfit(r))^2;
        SST=SST+(decompose_zernike(j,r)-y_ave)^2;
    end
    R_square = 1-SSE/SST;
    
    figure
    scatter(x_plot, decompose_zernike(j,:),10)
    box on
    set(gca,'FontSize',16)
    colormap Jet
    title(strcat(num2str(j),'th Zernike term'),'FontSize',16)
    xlabel('Applied [nm]','FontSize',16)
    ylabel('Measured [nm]','FontSize',16)
    hold on
    %plot linear fitting
    plot(x_plot,yfit,'LineWidth',2);
    % add the linear relation to figure
    txt1 = strcat('r^2=', num2str(R_square),'     y=', num2str(P(1),4),'*x+',num2str(P(2),4));
    text(0,0,txt1,'FontSize',16);
    hold off
    Rmatrix(j,:)=[P(1)^(-1) P(2)];
end 
save(strcat(header1,'_Rmatrix','.mat'),'Rmatrix');
