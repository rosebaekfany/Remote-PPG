D2_DATE = meanRefLowosRef;
X = (D2_DATE(:,1))' ;
Y = (D2_DATE(:,3))' ;
n= length(X) ;
T= 0:1:n-1 ;
%% 
%for 1D-x
[pks2,locs2] = findpeaks(-X,T,'MinPeakDistance',500);
figure;
findpeaks(-X,T,'MinPeakDistance',500);
text(locs2+.02,pks2,num2str((1:numel(pks2))'));
pks = [];
locs = [];
for i = 1:length(pks2)
    if i < length(pks2)
        if abs(pks2(i)-pks2(i+1))<0.02 
            pks = [pks , pks2(i)];
            locs = [locs , locs2(i)];
        end
    end
end
time_diff = diff(locs);
matrices = {X(locs(1):locs(2)), X(locs(2):locs(3)), X(locs(3):locs(4)), X(locs(4):locs(5)), X(locs(5):locs(6)), X(locs(6):locs(7))};
max_size = size(matrices{1});
for i = 2:numel(matrices)
    current_size = size(matrices{i});
    max_size = max(max_size, current_size);
end
reffrence_pulse=matrices{4};
aligned_matrices_X = cell(size(matrices));
for i = 1:numel(matrices)
    [aligned_matrices_X{i},not] = alignsignals(matrices{i}, reffrence_pulse);
end
padded_matrices_X = cell(size(aligned_matrices_X));
for i = 1:numel(aligned_matrices_X)
    current_size = size(aligned_matrices_X{i});
    padded_matrices_X{i} = padarray(aligned_matrices_X{i}, max_size - current_size, NaN, 'post');
end
numeric_array = cat(3, padded_matrices_X{:});
average_matrix_X = nanmean(numeric_array, 3);
% 
% for 1D-y

[pks,locs] = findpeaks(-Y,T,'MinPeakDistance',1000);
figure;
findpeaks(-Y,T,'MinPeakDistance',1000);
text(locs+.02,pks,num2str((1:numel(pks))'))
time_diff = diff(locs);
matrices = {Y(locs(2):locs(3)), Y(locs(3):locs(4)), Y(locs(4):locs(5)), Y(locs(5):locs(6)), Y(locs(6):locs(7)), Y(locs(7):locs(8))};
max_size = size(matrices{1});
for i = 2:numel(matrices)
    current_size = size(matrices{i});
    max_size = max(max_size, current_size);
end
reffrence_pulse=matrices{4};
aligned_matrices_Y = cell(size(matrices));
for i = 1:numel(matrices)
    [aligned_matrices_Y{i},not] = alignsignals(matrices{i}, reffrence_pulse);
end
padded_matrices_Y = cell(size(aligned_matrices_Y));
for i = 1:numel(aligned_matrices_Y)
    current_size = size(aligned_matrices_Y{i});
    padded_matrices_Y{i} = padarray(aligned_matrices_Y{i}, max_size - current_size, NaN, 'post');
end
numeric_array = cat(3, padded_matrices_Y{:});
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
[cx,lagsx] = xcorr(X,average_matrix_X);
plot(lagsx,cx);
figure;
[cy,lagsy] = xcorr(Y,average_matrix_Y);
plot(lagsy,cy);
%% 
coefficient_Y = [];
coefficient_X = [];
conv_self_Y = max(conv(average_matrix_Y,average_matrix_Y));
conv_self_X = max(conv(average_matrix_X,average_matrix_X));
for i = 1:numel(padded_matrices_Y)
    conv_Y=conv(average_matrix_Y,padded_matrices_Y{i});
    coefficient_Y = [coefficient_Y , 100*abs(1-max(conv_Y)/conv_self_Y)] ;
end

for i = 1:numel(padded_matrices_X)
    conv_X=conv(average_matrix_X,padded_matrices_X{i});
    coefficient_X = [coefficient_X , 100*abs(1-max(conv_X)/conv_self_X)] ;
end
figure;
stem(coefficient_X);
ylabel('Difference coefficient_X %');
disp(mean(coefficient_X));
figure;
stem(coefficient_Y);
ylabel('Difference coefficient_Y %');
disp(mean(coefficient_Y));

%% 
D4_data = reginalRefSignalAll;
coefficient_all = zeros(6,7);
for i=1:6
    for j=1:7
        coefficient_R = 0 ;
        R= squeeze(D4_data(i,j,1,:));
        [pks,locs] = findpeaks(-R,T,'MinPeakDistance',2000);
        time_diff = diff(locs);
        matrices = {R(locs(2):locs(3)),R(locs(3):locs(4)),R(locs(4):locs(5)),R(locs(5):locs(6))};

        for h = 1:numel(matrices)
        conv_Y=conv(average_matrix_Y,matrices{h});
        coefficient_R = coefficient_R + 100*abs(1-max(conv_Y)/conv_self_Y) ;
        end
        coefficient_R = coefficient_R/numel(matrices) ;
        coefficient_all(i,j)= coefficient_R ;
    end
end
disp(coefficient_all);
figure;
imagesc(coefficient_all);
colorbar;
%% 
conv_Y=conv(average_matrix_Y,padded_matrices_Y{1});
ny= length(conv_Y) ;
Ty= 0:1:ny-1 ;

%% Week 3
%%
coefficient_Cc_Y = [];
coefficient_Cc_X = [];
for i = 1:numel(padded_matrices_Y)
    Cc_Y=corrcoef(average_matrix_Y,padded_matrices_Y{i},'rows','complete');
    coefficient_Cc_Y = [coefficient_Cc_Y , Cc_Y(1,2)] ;
end

for i = 1:numel(padded_matrices_X)
    Cc_X=corrcoef(average_matrix_X,padded_matrices_X{i},'rows','complete');
    coefficient_Cc_X = [coefficient_Cc_X , Cc_X(1,2)] ;
end

disp(mean(coefficient_Cc_X));
disp(mean(coefficient_Cc_Y));
%% 
coefficient_all_Cc = zeros(6,7);
for i=1:6
    for j=1:7
        coefficient_R = 0 ;
        R= squeeze(D4_data(i,j,2,:));
        [pks,locs] = findpeaks(-R,T,'MinPeakDistance',2000);
        time_diff = diff(locs);
        matrices = {R(locs(2):locs(3)),R(locs(3):locs(4)),R(locs(4):locs(5)),R(locs(5):locs(6))};
        for h = 1:numel(matrices)
            current_size = size(matrices{h}');
            max_size = max(size(average_matrix_Y), current_size);
            average_matrix_Y_Cc = padarray(average_matrix_Y, max_size - size(average_matrix_Y), NaN, 'post');
            matrices_i = padarray(matrices{h}', max_size - size(matrices{h}'), NaN, 'post');
            Cc=corrcoef(average_matrix_Y_Cc,matrices_i,'rows','complete');
            coefficient_R = coefficient_R + Cc(1,2) ;
        end
        coefficient_R = coefficient_R/numel(matrices) ;
        coefficient_all_Cc(i,j)= coefficient_R ;
    end
end

disp(coefficient_all_Cc);
figure;
imagesc(coefficient_all_Cc);
colorbar;
%% 
