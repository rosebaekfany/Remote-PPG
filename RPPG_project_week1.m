D2_DATE = meanRefLowosRef;
X = (D2_DATE(:,1))' ;
Y = (D2_DATE(:,3))' ;
n= length(X) ;
T= 0:1:n-1 ;
plot(T,X);
hold on;
plot(T,Y);
xlabel('fram');
ylabel('y');
title('2D plot');
grid on;
%% 
figure
mat = meanRefLowosRef ;
t = 0:numel(mat)-1;                                                                                        
L = numel(t); 
mat2=mat(:,2:3);
imagesc(mat2);
colormap jet;
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('2D Signal');
%% 
fs=2000 ;
f= fs/2*linspace(-1,1,length(X));
% plot(f, abs(fftshift(fft(X)/L)));
mat2=mat(:,2:3);
FTmat = fft(mat2)/L;
figure
plot(f, abs(fftshift(FTmat)));
xlim([-100,100]);
grid
title('Two-Sided Fourier Transform');
xlabel('Frequency');
ylabel('Amplitude');
figure;
frequencyDomain = fftshift(fft2(mat2));
imagesc(f, f, abs(frequencyDomain));
colorbar;
xlabel('Frequency 1 (Hz)');
ylabel('Frequency 2 (Hz)');
title('Frequency Domain of the 2D Signal');
%% 
cwtout = cwtft2(mat2);
sca = 1;
imagesc(abs(cwtout.cfs(:,:,1,1,sca)));
colormap jet;
colorbar;
xlabel('Time (X-axis)');
ylabel('Scale (Y-axis)');
title('CWT of 2D Signal');
figure 
cwtmexh = cwtft2(mat2,'wavelet','mexh','scales',1,...
    'angles',[0 pi/2]);
surf(real(cwtmexh.cfs(:,:,1,1,1)));
shading interp; 
colormap jet;
colorbar;
figure
cwt(X);

%% 
%for 1D-x
[pks,locs] = findpeaks(X,T,'MinPeakDistance',1000);
figure;
findpeaks(X,T,'MinPeakDistance',1000);
text(locs+.02,pks,num2str((1:numel(pks))'))
time_diff = diff(locs);
matrices = {X(locs(1):locs(2)), X(locs(2):locs(3)), X(locs(3):locs(4)), X(locs(4):locs(5)), X(locs(5):locs(6)), X(locs(6):locs(7)), X(locs(7):locs(8)), X(locs(8):locs(9))};
max_size = size(matrices{1});
for i = 2:numel(matrices)
    current_size = size(matrices{i});
    max_size = max(max_size, current_size);
end
padded_matrices = cell(size(matrices));
for i = 1:numel(matrices)
    current_size = size(matrices{i});
    padded_matrices{i} = padarray(matrices{i}, max_size - current_size, NaN, 'post');
end
numeric_array = cat(3, padded_matrices{:});
average_matrix_X = nanmean(numeric_array, 3);

%for 1D-y

[pks,locs] = findpeaks(Y,T,'MinPeakDistance',1000);
figure;
findpeaks(Y,T,'MinPeakDistance',1000);
text(locs+.02,pks,num2str((1:numel(pks))'))
time_diff = diff(locs);
matrices = {X(locs(1):locs(2)), X(locs(2):locs(3)), X(locs(3):locs(4)), X(locs(4):locs(5)), X(locs(5):locs(6)), X(locs(6):locs(7)), X(locs(7):locs(8)), X(locs(8):locs(9))};
max_size = size(matrices{1});
for i = 2:numel(matrices)
    current_size = size(matrices{i});
    max_size = max(max_size, current_size);
end
padded_matrices = cell(size(matrices));
for i = 1:numel(matrices)
    current_size = size(matrices{i});
    padded_matrices{i} = padarray(matrices{i}, max_size - current_size, NaN, 'post');
end
numeric_array = cat(3, padded_matrices{:});
average_matrix_Y = nanmean(numeric_array, 3);

%ploting
ny= length(average_matrix_Y) ;
Ty= 0:1:ny-1 ;
nx= length(average_matrix_X) ;
Tx= 0:1:nx-1 ;
figure;
ax1 = nexttile;
plot(ax1,Ty,average_matrix_Y);
xlabel('Time');
ylabel('Amplitude Y');
ax2 = nexttile;
plot(ax2,Tx,average_matrix_X);
linkaxes([ax1,ax2],'x');
xlabel('Time');
ylabel('Amplitude X');

%% 
figure;
S=Y.*X;
volume = cumsum(mat2);
Tv= 0:1:length(volume)-1 ;
plot(Tv,volume(:,1));
hold on;
plot(Tv,volume(:,2));
xlabel('T');
title('integerate 1D');
figure;
volume2 = cumsum(S);
plot(Tv,volume2);
xlabel('T');
title('integerate 2D -X*Y');

