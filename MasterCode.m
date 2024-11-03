% Specify the filename
filename = 'Group0_Case1.xls';

% Check if the file exists
if isfile(filename)
    % Load the data from the Excel file with VariableNamingRule set to preserve
    try
        opts = detectImportOptions(filename, 'VariableNamingRule', 'preserve');
        data = readtable(filename, opts);
        fprintf('Successfully read %s\n', filename);
    catch ME
        fprintf('Error reading file %s: %s\n', filename, ME.message);
        return;
    end
    
    % Display the variable names to help with debugging
    disp(data.Properties.VariableNames);
    
    % Extract the time and acceleration data using the exact column names
    try
        time = data{:, 'Time (s)'}; % Adjust the column name as per your data
        accX = data{:, 'Linear Acceleration x (m/s^2)'}; % Adjust the column name as per your data
        accY = data{:, 'Linear Acceleration y (m/s^2)'};
        accZ = data{:, 'Linear Acceleration z (m/s^2)'};
        absAcc = data{:, 'Absolute acceleration (m/s^2)'};
    catch ME
        fprintf('Error extracting data from file %s: %s. Check column names.\n', filename, ME.message);
        return;
    end
    
    % Check if time vector is valid
    if length(time) < 2
        fprintf('Invalid time data in file %s. Skipping.\n', filename);
        return;
    end
    
    % Design a low-pass filter to remove high-frequency noise
    Fs = 1 / mean(diff(time)); % Sampling frequency
    
    % Check if Fs is valid
    if isnan(Fs) || isinf(Fs) || Fs <= 0
        fprintf('Invalid sampling frequency in file %s. Skipping.\n', filename);
        return;
    end
    
    % Define initial cutoff frequency
    initialCutoffFrequency = 50; % Initial cutoff frequency in Hz
    nyquistFrequency = Fs / 2; % Nyquist frequency
    
    % Adjust cutoff frequency if it's too high
    if initialCutoffFrequency >= nyquistFrequency
        cutoffFrequency = 0.9 * nyquistFrequency;
        fprintf('Cutoff frequency adjusted to %.2f Hz\n', cutoffFrequency);
    else
        cutoffFrequency = initialCutoffFrequency;
    end
    
    [b, a] = butter(4, cutoffFrequency / nyquistFrequency, 'low'); % 4th order Butterworth filter
    
    % Apply the filter to the data
    accX_filtered = filtfilt(b, a, accX);
    accY_filtered = filtfilt(b, a, accY);
    accZ_filtered = filtfilt(b, a, accZ);
    absAcc_filtered = filtfilt(b, a, absAcc);
    
    % Perform FFT on Linear Acceleration in X (Filtered)
    N = length(time);
    f = Fs * (0:(N/2)) / N; % Frequency range
    fftAccX_filtered = fft(accX_filtered);
    P2 = abs(fftAccX_filtered / N);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Create a table with the processed data
    processedTable = table(time, accX, accY, accZ, absAcc, accX_filtered, accY_filtered, accZ_filtered, absAcc_filtered, ...
                           'VariableNames', {'Time', 'accX', 'accY', 'accZ', 'absAcc', 'AccX_Filtered', 'AccY_Filtered', 'AccZ_Filtered', 'AbsAcc_Filtered'});
    
    % Add group number and case number columns
    processedTable.GroupNumber = zeros(size(processedTable, 1), 1);
    processedTable.CaseNumber = ones(size(processedTable, 1), 1);
    
    % Prompt user to save the output file
    [file, path] = uiputfile('Group0_Case1_fft.csv', 'Save Processed Data As');
    if isequal(file, 0) || isequal(path, 0)
        fprintf('User canceled the file save operation.\n');
    else
        newFilename = fullfile(path, file);
        try
            writetable(processedTable, newFilename);
            fprintf('Processed data saved to %s\n', newFilename);
        catch ME
            fprintf('Error writing to file %s: %s\n', newFilename, ME.message);
        end
    end
else
    fprintf('File %s does not exist.\n', filename);
end
