function [ ] = play( path, fr )
%PLAY Summary of this function goes here
%   Detailed explanation goes here
    
    disp(['Play ', path]);
    
    % Open output .bin file
    bin_pathname = [path, '.bin'];
    binFid = fopen(bin_pathname, 'r', 'l');
    % Read header
    binHeader = fread(binFid, 4, 'int16');
    width = binHeader(1);
    height = binHeader(2);
    % nb_images = header(3);
    % nb_bits = header(4);
    
    % Open output .vec file
    vec_pathname = [path, '.vec'];
    vecFid = fopen(vec_pathname, 'r', 'l');
    vecHeader = fscanf(vecFid, '%g %g %g %g %g\n', 5);
    nb_frames = vecHeader(2);
    
    disp(['  number of frames  : ', num2str(nb_frames)]);
    disp(['  sampling frequency: ', num2str(fr), ' Hz']);
    disp(['  duration          : ', num2str(nb_frames / fr), ' sec']);
    
    try
        figure();
        colormap(gray(256));
        % Display first frame.
        vecFrame = fscanf(vecFid, '%g %g %g %g %g\n', 5);
        imageId = vecFrame(2) + 1;
        headerSize = 8;
        imageSize = width * height;
        offset = headerSize + (imageId - 1) * imageSize;
        fseek(binFid, offset, 'bof');
        frame = fread(binFid, imageSize, 'uint8');
        frame = reshape(frame, width, height);
        frame = frame.';
        hi = image(frame);
        ht = title(['Image ', num2str(imageId + 1), ', Frame 1/', num2str(nb_frames)]);
        for frameId = 2:nb_frames
            % Diplay next frame
            vecFrame = fscanf(vecFid, '%g %g %g %g %g\n', 5);
            imageId = vecFrame(2) + 1;
            headerSize = 8;
            imageSize = width * height;
            offset = headerSize + (imageId - 1) * imageSize;
            fseek(binFid, offset, 'bof');
            frame = fread(binFid, imageSize, 'uint8');
            frame = reshape(frame, width, height);
            frame = frame.';
            set(hi, 'CData', frame);
            set(ht, 'String', ['Image ', num2str(imageId + 1), ', Frame ', num2str(frameId), '/', num2str(nb_frames)]);
            pause(1 / fr);
            drawnow;
            % if mod(frameId, fr) == 15
            %     id = floor(frameId / fr);
            %     disp(all(frame(:) == 128));
            %     frame_ = frame / 255;
            %     imwrite(frame_, ['~/data/dmd_stimuli/snapshots/', num2str(id),'.bmp'], 'bmp');
            % end
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:class:InvalidHandle')
            ;
        else
            disp(ME);
        end
    end
    
    return
    
end

