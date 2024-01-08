%% 
%clc;
%clear;
person = RawSignalOri;

%% 
regionLeny = 10;
regionLenx = 10; 
m = 1:regionLeny:size(RawSignalOri,1)+1; % height or row = y
n = 1:regionLenx:size(RawSignalOri,2)+1; % width or column= x
vertical = 1;
filtered_person = zeros(length(m),length(n),3,19999);
for i = 1:length(m)-1 % by row
    horizontal = 1;
    for j = 1:length(n)-1  % by column
        reginalSignal = BUTTER3D(double(squeeze(mean(RawSignalOri(m(i):m(i+1)-1,n(j):n(j+1)-1,:,:),[1 2]))'),fps);
        filtered_person(horizontal,vertical,:,:)=reginalSignal'  ;
        horizontal = horizontal +1;
    end 
    vertical = vertical +1;
end 
%% 
n= 19999 ;
T= 0:1:n-1 ;
figure;
plot(T,squeeze(RawSignalOri(15,15,2,:)));
title('Raw Signal(15,15,2,:)');
xlabel('Fram');


n= 19999 ;
T= 0:1:n-1 ;
figure;
plot(T,squeeze(filtered_person(15,15,2,:)));
title('filtered Signal(15,15,2,:)');
xlabel('Fram');

%%
max_size=10000;
average_all = zeros(32,51,max_size);
average_signal = zeros(max_size,1) ;
for h = 1:32
    for j = 1 :51
        X=squeeze(filtered_person(h,j,2,:));
        [pks,locs] = findpeaks(-X,T,'MinPeakDistance',1000);
        time_diff = diff(locs);
        if length(locs)>=6
            matrices = {X(locs(1):locs(2)),X(locs(2):locs(3)), X(locs(3):locs(4)), X(locs(4):locs(5))};
            reffrence_pulse=matrices{2};
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
            average_matrix_Y = nanmean(numeric_array,3);
            average_all(h,j,:) = average_matrix_Y(:,1);
            average_signal = (average_matrix_Y(:,1)+average_signal)/2;
        end  
    end
end
%% 
figure; 
ny1= length(average_all(15,15,:)) ;
Ty1= 0:1:ny1-1 ;
plot(Ty1,squeeze(average_all(15,15,:)));
title('Local mask(15,15,2,:)');
xlabel('Fram');

figure;
ny1= length(X(locs(1):locs(2))) ;
Ty1= 0:1:ny1-1 ;
plot(Ty1,squeeze(X(locs(1):locs(2))));
title('pulse one');
xlabel('Fram');

figure;
ny1= length(X(locs(3):locs(4))) ;
Ty1= 0:1:ny1-1 ;
plot(Ty1,squeeze(X(locs(3):locs(4))));
title('pulse three');
xlabel('Fram');

figure;
ny1= length(padded_matrices_X{1}) ;
Ty1= 0:1:ny1-1 ;
plot(Ty1,squeeze(padded_matrices_X{1}));
title('aligned and resized pulse one');
xlabel('Fram');

figure;
ny1= length(padded_matrices_X{3}) ;
Ty1= 0:1:ny1-1 ;
plot(Ty1,squeeze(padded_matrices_X{3}));
title('aligned and resized pulse three');
xlabel('Fram');


%% 

figure;
ny= length(average_signal) ;
Ty= 0:1:ny-1 ;
plot(Ty,average_signal);
title('Total mask');
xlabel('Fram');

figure;
ny= length(average_all(5,5,:)) ;
Ty= 0:1:ny-1 ;
plot(Ty,squeeze(average_all(5,5,:)));
title('Local mask(5,5,2,:)');
xlabel('Fram');

figure;
ny1= length(average_all(10,10,:)) ;
Ty1= 0:1:ny1-1 ;
plot(Ty1,squeeze(average_all(10,10,:)));
title('Local mask(10,10,2,:)');
xlabel('Fram');

figure;
ny2= length(average_all(20,20,:)) ;
Ty2= 0:1:ny2-1 ;
plot(Ty2,squeeze(average_all(20,20,:)));
title('Local mask(20,20,2,:)');
xlabel('Fram');


%% meth1

coefficient_all_Cc = zeros(32,51);
 for i=1:32
     for j=1:51
         coefficient_R = 0 ;
         R= squeeze(average_all(i,j,:));   
         Cc=corrcoef(average_signal,R,'rows','complete');
         coefficient_R = coefficient_R + Cc(1,2) ;
         coefficient_all_Cc(i,j)= coefficient_R ;
         
     end
end

figure;
imagesc(coefficient_all_Cc);
colorbar;
title('correlation coefficient between total mask and local masks');

% %% meth 2
%  coefficient_all_Cc = zeros(6,7);
%  for i=1:52
%      for j=1:51
%          coefficient_R = 0 ;
%          R= squeeze(filtered_person(i,j,2,:));
%          [pks,locs] = findpeaks(-R,T,'MinPeakDistance',1000);
%          time_diff = diff(locs);
%          if length(locs)>=6
%              matrices = {R(locs(2):locs(3)),R(locs(3):locs(4)),R(locs(4):locs(5)),R(locs(5):locs(6))};
%              for h = 1:numel(matrices)
%                  current_size = size(matrices{h}');
%                  max_size = max(size(average_signal), current_size);
%                  average_matrix_Y_Cc = padarray(average_signal, max_size - size(average_signal), NaN, 'post');
%                  matrices_i = padarray(matrices{h}', max_size - size(matrices{h}'), NaN, 'post');
%                  Cc=corrcoef(average_matrix_Y_Cc,matrices_i,'rows','complete');
%                  coefficient_R = coefficient_R + Cc(1,2) ;
%              end
%              coefficient_R = coefficient_R/numel(matrices) ;
%              coefficient_all_Cc(i,j)= coefficient_R ;
%          end
%      end
% end
% 
% imagesc(coefficient_all_Cc);
% colorbar;
% 
% %% meth3
% 
% coefficient_all_Cc = zeros(6,7);
%  for i=1:52
%      for j=1:51
%          coefficient_R = 0 ;
%          R= squeeze(filtered_person(i,j,2,:));
%          [pks,locs] = findpeaks(-R,T,'MinPeakDistance',1000);
%          time_diff = diff(locs);
%          if length(locs)>=6
%              matrices = {R(locs(2):locs(3)),R(locs(3):locs(4)),R(locs(4):locs(5)),R(locs(5):locs(6))};
%              for h = 1:numel(matrices)
%                  current_size = size(matrices{h}');
%                  max_size = max(size(average_all(i,j,:)), current_size);
%                  average_matrix_Y_Cc = padarray(average_all(i,j,:), max_size - size(average_all(i,j,:)), NaN, 'post');
%                  matrices_i = padarray(matrices{h}', max_size - size(matrices{h}'), NaN, 'post');
%                  Cc=corrcoef(average_matrix_Y_Cc,matrices_i,'rows','complete');
%                  coefficient_R = coefficient_R + Cc(1,2) ;
%              end
%              coefficient_R = coefficient_R/numel(matrices) ;
%              coefficient_all_Cc(i,j)= coefficient_R ;
%          end
%      end
% end
% 
% imagesc(coefficient_all_Cc);
% colorbar;



%% 

function signal1DButter = BUTTER1D(raw1Dsignal,fps)
    freqBounder = [0.5 3.5]; 
     % *** get ride of trend again, can delete!!
    [A,B,C,D] = butter(2,[freqBounder(1) freqBounder(2)]/(fps/2)); % IMPORTANT****butterworth not fit for narrow bandpass, or use low order for narrow band.
    [filter_SOS,g] = ss2sos(A,B,C,D);
    signal1DButter = filtfilt(filter_SOS,g,raw1Dsignal); % remove lag
end

function signal3DBut = BUTTER3D(raw3DSignal,fps)    % without movemean
    freqBounder = [0.5 3.5];
    for i = 1:3
%         raw1DSignal = raw3DSignal(:,i) - movmean(raw3DSignal(:,i),fps*0.8);
        [A,B,C,D] = butter(2,[freqBounder(1) freqBounder(2)]/(fps/2)); % IMPORTANT****butterworth not fit for narrow bandpass, or use low order for narrow band.
        [filter_SOS,g] = ss2sos(A,B,C,D);
        raw1DSignal = filtfilt(filter_SOS,g,raw3DSignal(:,i));
        signal3DBut(:,i) = raw1DSignal;
    end
end 

