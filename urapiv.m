%%%Last modified EB 12/14

function [varargout] = urapiv(dirname,dirout,name,method,itt,spc,Niteracion,s2nm,s2nl,sclt,outl,outG,crop_vector,isub,iord,subpR,subpr,npix,iref,vecind,ndiv)


% URAPIV - Program reads all TIFF files in provided directory (in pairs),
% and calculates the velocity field, according to the following steps:
%
% 1. Read images in TIFF formats (original version reads in BMP)
% 2. Calcuate Cross-correlation for each itterogation block, according to
%    the grid spacing
% 3. Calculates Signal-To-Noise ratio and removes problematic results
% 4. Removes outlayers with the velocity values bigger than the average of the
%    matrix times Outlayer-limit, that is provided by the user
% 5. Removes local outlayers by means of adaptive local median filter
% 6. Fills in removed vectors by interpolation of the neighboor vectors
% 7. Presents all the results of position (x,y)  velocity (vx,vy) and S2N
%
% Usage:
%
% 				UraPIV(DIR,ITT,SPC,S2NM,S2NL,SCLT,OUTL)
%
% Inputs:
%
%        DIR - name of the directory with session images, (string)
%              like '/u/liberzon/roi/tifs', or shorter './tifs'
%        ITT - intterogation size in pixels, same on in X and Y directions (have to be 2^N)
%        SPC - spacing is the grid of the image in pixels
%        S2NM - Signal-To-Noise-Ratio calculation method:
%                1 - Peak-by-peak search with 3x3 pixels kernel
%                  or
%                2 - Peak-to-mean ratio in one interrogation area
%        S2NL - limit of signal-to-noise-ratio
%        SCLT - scaling*time in units [meter/pixel/seconds] (e.g., 1e5/26200)
%
%        OUTL - global filter, taking apart values that are bigger
%               that OUTL times the average value of the whole matrix.
%        CROPVEC - Vector of crop values in [left top right bottom] order,
%                  run URAPIV without this argument to enter it manually after
%                  seeing your images.
%
%
%   Usage:
%       vel = urapiv('./images',32,32,2,1,1,10,[0 0 0 0]);
%
% Authors: Alex Liberzon and Roi Gurka
% Co-Author: Uri Shavit
%
% Started at: 30-Jul-98


warning off
more off

% All arguments should be supplied, no defaults are defined:
if nargin < 7
    error('Usage: vel = urapiv(''./images'',32,32,2,1,1,10,[0 0 0 0]);');
end

% Argument check:
if s2nm ~=1 & s2nm ~= 2
    s2nm = 2;
end

[filenames,amount,ext] = ReadImDir(dirname,name,'TIF');
% [filenames,amount,ext] = ReadImDir(dirname,name,'tif');% % % for janice >>> USE THIS LINE AND COMMENT THE ONE ABOVE (use .tif and not .TIF)
% filenames, amount, ext,
filenames(1),
[a,b] = read_pair_of_images(fullfile(dirname,filenames{iref(1)}),fullfile(dirname,filenames{2}),[0 0 0 0],itt,spc,ext);

if isempty(a) | isempty(b)
    error('Something wrong with your images')
end

if ~exist('crop_vector','var')
    ah = figure;imshow(a),
    axis on,ax=axis; grid on; set(gca,'xtick',[0:itt:ax(2)],'xticklabel',[],...
        'ytick',[0:itt:ax(4)],'yticklabel',[]);
    bh = figure;imshow(b);
    axis on,ax=axis; grid on; set(gca,'xtick',[0:itt:ax(2)],'xticklabel',[],...
        'ytick',[0:itt:ax(4)],'yticklabel',[]);
    disp('')
    crop_vector = input(sprintf('%s \n %s \t ',...
        'Enter the number of interrogation lines to crop',...
        '[Left,Top,Right,Bottom], Enter for none '));
    if isempty(crop_vector),crop_vector = zeros(4,1); end

end

% matrix "veloci" includes the results of all the pairs of images
% veloci = zeros(;
% veloci = []; % too expensive in memory/speed
% learn pre-allocation:
% Prepare the results storage;

reslenx = (sx-itt)/spc+1;
resleny = (sy-itt)/spc+1;
%veloci = zeros(reslenx*resleny*amount,5);
velociInd = 0;

nfiles = length(iref);
for i=1:nfiles
    imf(i).A = fullfile(dirname,filenames{iref(i)});
    ima(i).AA = read_image(imf(i).A,crop_vector,itt,spc,ext);
end

if vecind==0
    vecind = 1:amount;
end
for fileind = vecind	% main loop, for whole file list

    velociInd = velociInd + 1;

    for i=1:nfiles

        fullfile(dirname,filenames{fileind}),
        imf(i).B = fullfile(dirname,filenames{fileind});

        [ima(i).B]=read_image(imf(i).B,crop_vector,itt,spc,ext);
	if     ndiv == 4;
            [ima(i).A1,ima(i).A2,ima(i).A3,ima(i).A4,ima(i).B] = im_interp_correc_quad ...
              (ima(i).AA,ima(i).B,s2nm,isub,iord,subpR,subpr,npix);
	elseif ndiv == 1;
            [ima(i).A,ima(i).B] = im_interp_correc(ima(i).AA,ima(i).B,s2nm,isub,iord,subpR,subpr,npix);
	end;
    end

    %    [a1,a2,a3,b1,b2,b3] = read_3_images(imageA1,imageA2,imageA3,imageB1,imageB2,imageB3,crop_vector,itt,spc,ext);

    [sx,sy]= size(ima(1).AA);

    %%%%%% Start the loop for each interrogation block %%%%%%%

    iteracion = 0;
    while (iteracion<Niteracion)
        iteracion = iteracion + 1;

        itt0 = (Niteracion - iteracion + 1)*itt;
        spc0 = (Niteracion - iteracion + 1)*spc;
        Nfft0 = 2*itt0;

        reslenx = floor((sx-itt0)/spc0+1);
        resleny = floor((sy-itt0)/spc0+1);

        res = zeros(reslenx*resleny,5);

        for i=1:nfiles
            aw(i).A = zeros(itt0);
            aw(i).B = zeros(itt0);
            cw(i).c = zeros(Nfft0);
        end

        % a1_2 = zeros(itt0);
        % a2_2 = zeros(itt0);
        % a3_2 = zeros(itt0);
        % b1_2 = zeros(itt0);
        % b2_2 = zeros(itt0);
        % b3_2 = zeros(itt0);

        %       c1 = zeros(Nfft0);
        %       c2 = zeros(Nfft0);
        %       c3 = zeros(Nfft0);

        resind = 0;
        xind = [1:spc0:sx-itt0+1];
        yind = [1:spc0:sy-itt0+1];
        for k=xind
            disp(sprintf('\n Working on %d pixels row',k))
            for m=yind
                % Remove following line if you like 'silent' run
                %               fprintf(1,'.');

                % interpolates from previous iteration

                if (iteracion==1);
                    shiftx = 0;
                    shifty = 0;
                else
                    shift = interp2(xind0,yind0,vector,k,m,'cubic',0);
                    shiftx = real(shift);
                    shifty = imag(shift);
                end

                shifti = round(shiftx);
                shiftj = round(shifty);

                vecxa = [k:k+itt0-1];
                vecya = [m:m+itt0-1];

                vecxb = [k:k+itt0-1]-shifti;
                vecyb = [m:m+itt0-1]-shiftj;

                if vecxb(end)>sx
                    vecxb = [sx-itt0+1:sx];
                end
                if vecxb(1)<1
                    vecxb = [1:itt0];
                end

                if vecyb(end)>sy
                    vecyb = [sy-itt0+1:sy];
                end
                if vecyb(1)<1
                    vecyb = [1:itt0];
                end

                if method==1; c  = zeros(Nfft0); end
                if method==2; c  = zeros(itt0); end
                for i=1:nfiles
                    [H W] = size(ima(1).AA); 
		    if ndiv == 4;
                        % quadrants
                        if k<=H/2 & m<=W/2	%1
                            aw(i).A = ima(i).A1(vecxa,vecya);
                        end
                        if k>H/2 & m<=W/2	%2
                            aw(i).A = ima(i).A2(vecxa,vecya);
                        end
                        if k<=H/2 & m>W/2	%3
                            aw(i).A = ima(i).A3(vecxa,vecya);
                        end
                        if k>H/2 & m>W/2	%4
                            aw(i).A = ima(i).A4(vecxa,vecya);
                        end
		    elseif ndiv == 1;
                        aw(i).A = ima(i).A(vecxa,vecya);
		    end;

                    aw(i).B = ima(i).B(vecxa,vecya);

                    if (method == 1)
                        cw(i).c = cross_correlate(aw(i).A,aw(i).B,Nfft0);
                    end
                    if (method == 2)
                        cw(i).c = mqd(aw(i).A,aw(i).B,itt0,itt0);
                    end

                    c = c + cw(i).c;

                end
                c = c/nfiles;

                %               a1_2 = a1(vecxa,vecya);
                %               a2_2 = a2(vecxa,vecya);
                %               a3_2 = a3(vecxa,vecya);
                %
                %               b1_2 = b1(vecxb,vecyb);
                %               b2_2 = b2(vecxb,vecyb);
                %               b3_2 = b3(vecxb,vecyb);

                %	       if method==1
                %                  c1 = cross_correlate(a1_2,b1_2,Nfft0);
                %                  c2 = cross_correlate(a2_2,b2_2,Nfft0);
                %                  c3 = cross_correlate(a3_2,b3_2,Nfft0);
                %	       end
                %
                %               if method==2
                %                  c1 = mqd(a1_2,b1_2,itt0);
                %                  c2 = mqd(a2_2,b2_2,itt0);
                %                  c3 = mqd(a3_2,b3_2,itt0);
                %	       end
                %
                %               c = (c1+c2+c3)/3;

                %	       for i=1
                %%                 jvec=1:Nfft0;
                %                jvec=floor([Nfft0/2-10:Nfft0/2+10]);
                %                 subplot(2,1,1)
                %                 imagesc(c(jvec,jvec));
                %                 axis equal;  axis tight
                %                 cax = caxis;
                %                 colorbar
                %                 subplot(2,1,2)
                %                 imagesc(cw(i).c(jvec,jvec));
                %                 axis equal;  axis tight
                %                 caxis(cax);
                %                 colorbar
                %                axis equal; axis tight;
                %		pause
                %%                  subplot(2,2,3)
                %%                 imagesc(c2(jvec,jvec));
                %%%                 caxis(cax);
                %%                 colorbar
                %%                 axis equal; axis tight;
                %%                  subplot(2,2,4)
                %%                 imagesc(c3(jvec,jvec));
                %%%                 caxis(cax);
                %%                 colorbar
                %%                 axis equal; axis tight;
                %                 pause
                %	       end


                [peak1,peak2,pixi,pixj] = find_displacement(c,s2nm);

                [peakx,peaky,s2n] = sub_pixel_velocity(c,pixi,pixj,peak1,peak2,s2nl,sclt,itt0,itt0,isub,iord,subpR,subpr);

                % Scale the pixel displacement to the velocity
                u = (itt0+1-peaky+shifti)*sclt;
                v = (itt0+1-peakx+shiftj)*sclt;

                x = m+itt0/2-1;
                y = k+itt0/2-1;

                resind = resind + 1;
                res(resind,:) = [x y u v s2n];

            end
        end

        % NO_FILT_RES will be stored in '.._noflt.txt' file at the end of program
        no_filt_res = res;
        %% Unfiltered, uninterpolated: (comment with % sign if you don't need it)
        %fid = fopen([fullfile(dirname,filenames{fileind}(1:end-4)),'_noflt.txt'],'w');
        %fprintf(fid,'%3d %3d %7.4f %7.4f %7.4f\n',no_filt_res');
        %fclose(fid);

        % Reshape U and V matrices in two-dimensional grid and produce
        % velocity vector in U + i*V form (real and imaginary parts):

        u = reshape(res(:,3),resleny,reslenx);
        v = reshape(res(:,4), resleny,reslenx);
        vector = u + sqrt(-1)*v;

        %
        %      FILTERING PROCEDURES AND REMOVAL OF OUTLIERS ONLY IF OUTL>0
        %

        if outG>0			% filtering

            % Remove outliers - GLOBAL FILTERING
            vector(abs(vector)>mean(abs(vector(find(vector))))*outG) = 0;
            vector(vector~=vector)=0;
            u = real(vector);
            v = imag(vector);

        end
        if outl>0

            % Adaptive Local Median filtering

            w = 3; rad = 1;
            kernel = zeros(2*w+1);
            B = sum(exp(-([-w:w].^2)/(rad^2)))^2;
            for i=1:2*w+1;
                kernel(i,:) =  exp(-((i-w-1)^2+[-w:w].^2)/(rad^2))/B - 1./(2*w+1).^2;
            end

            tmpv = abs(conv2(v,kernel,'same'));
            tmpu = abs(conv2(u,kernel,'same'));

            % WE HAVE TO DECIDE WHICH LIMIT TO USE:
            % 1. Mean + 3*STD for each one separately OR
            % 2. For velocity vector length (and angle)
            % 3. OR OTHER.

            lmtv = mean(tmpv(find(tmpv))) + outl*std(tmpv(find(tmpv)));
            lmtu = mean(tmpu(find(tmpu))) + outl*std(tmpu(find(tmpu)));
            u_out = find(tmpu>lmtu);
            v_out = find(tmpv>lmtv);

            % Let's throw the outlayers out:
            u(u_out) = 0; u(v_out) = 0;
            v(v_out) = 0; v(u_out) = 0;
            vector = u + sqrt(-1)*v;

            res(:,3) = reshape(real(vector),resleny*reslenx,1);
            res(:,4) = reshape(imag(vector),resleny*reslenx,1);

            % Filtered results will be stored in '.._flt.txt' file
            filt_res = res;
            % Filtered, but not interpolated:
            fid = fopen([fullfile(dirname,filenames{fileind}(1:end-4)),'_flt.txt'],'w');
            fprintf(fid,'%3d %3d %7.4f %7.4f %7.4f\n',filt_res');
            fclose(fid);

            % Interpolation of the data:

            [indx,indy] = find(abs(vector)==0);

            while ~isempty(indx)

                for z=1:length(indx)
                    k = [max(3,indx(z))-2:min(resleny-2,indx(z))+2];
                    m = [max(3,indy(z))-2:min(reslenx-2,indy(z))+2];
                    tmpvec = vector(k,m);
                    tmpvec = tmpvec(find(tmpvec));
                    vector(indx(z),indy(z)) = mean(real(tmpvec))+ sqrt(-1)*mean(imag(tmpvec));
                end
                [indx,indy] = find(abs(vector)==0);
            end
            res(:,3) = reshape(real(vector),resleny*reslenx,1);
            res(:,4) = reshape(imag(vector),resleny*reslenx,1);

        end 					% end filtering

        xind0 = xind;
        yind0 = yind;
    end  % end of multigrid iteration

    % Save results as ASCII (text) files:
    % Final (filtered, interpolated) results
    fid = fopen([fullfile(dirout,filenames{fileind}(1:end-4)),'.txt'],'w'); fid,
    fprintf(fid,'%3d %3d %7.4f %7.4f %7.4f\n',res');
    fclose(fid);

    % Results visualization
    % Only for final, filtered and interpolated data
    %figure; imshow(a);
    %hold on
    %quiverm(res,'r');

    % Results for operating in other applications
    %veloci(reslenx*resleny*(velociInd-1)+1:velociInd*reslenx*resleny,:) = res;

end % of the loop

%if nargout == 1
%    varargout{1} = veloci;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    EXTERNAL FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [B] = read_image(image2,crop_vector,itt,spc,ext)
        % READ__IMAGE - reads 1 image (image2) as tif file
        % and crops it according to 'crop_vector'
        % Inputs:
        %         image2 - tif file names (string)
        %         crop_vector - 4 x 1 vector of follwoing values:
        %         [left,top,right,bottom] - each value is a number of lines
        %                                   of interrogation areas (ITTxITT pixels)
        %                                   which should be removed before the analysis.
        %         itt - interrogation area size in pixels
        %         spc - grid spacing (overlapping) size in pixels
        %	  A   - reference image
        %

        B = imread(image2,ext);
        B = double(B(:,:,1));

        [sxb,syb]=size(B);
        % A & B matrices HAVE to be of the same size, we take smallest:
        sx = min(sx,sxb); sy = min(sy,syb);

        % Crop the images to the desired size and
        % cut the last couple of pixels, so we'll get the
        % integer number of interrogation areas
        %
        %       ---- t ---
        %      |          |
        %      |          |
        %      l          r
        %      |          |
        %      |          |
        %       --- b ----
        %
        %
        l = crop_vector(1); % left side of the image
        t = crop_vector(2); % top side of the image
        r = crop_vector(3); % right side of the image
        b = crop_vector(4); % bottom of the image

        B = B(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);
    end

    function [A,B] = read_pair_of_images(image1,image2,crop_vector,itt,spc,ext)
        % READ_PAIR_OF_IMAGES - reads two images (image1,image2) as tif files
        % and crops them according to 'crop_vector'
        % Inputs:
        %         image1,image2 - tif file names (string)
        %         crop_vector - 4 x 1 vector of follwoing values:
        %         [left,top,right,bottom] - each value is a number of lines
        %                                   of interrogation areas (ITTxITT pixels)
        %                                   which should be removed before the analysis.
        %         itt - interrogation area size in pixels
        %         spc - grid spacing (overlapping) size in pixels
        %
        % Authors: Alex Liberzon & Roi Gurka
        % Date: 20-Jul-99
        % Last modified:
        % Copyright(c) 1999, Alex Liberzon

        A = imread(image1,ext);
        B = imread(image2,ext);
        A = double(A(:,:,1));
        B = double(B(:,:,1));

        %         A = rgb2gray(imread(image1,'tif'));
        %         B = rgb2gray(imread(image2,'tif'));
        [sx,sy]=size(A);[sxb,syb]=size(B);
        % A & B matrices HAVE to be of the same size, we take smallest:
        sx = min(sx,sxb); sy = min(sy,syb);

        % Crop the images to the desired size and
        % cut the last couple of pixels, so we'll get the
        % integer number of interrogation areas
        %
        %       ---- t ---
        %      |          |
        %      |          |
        %      l          r
        %      |          |
        %      |          |
        %       --- b ----
        %
        %
        l = crop_vector(1); % left side of the image
        t = crop_vector(2); % top side of the image
        r = crop_vector(3); % right side of the image
        b = crop_vector(4); % bottom of the image

        A = A(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);
        B = B(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);
        [sx,sy]=size(A);
    end


    function [im,im0] = im_interp_correc(im,im0,s2nm,isub,iord,subpR,subpr,npix)
        if (npix>0)
            se=strel('disk',npix+1);nhood=double(getnhood(se));

imm0 = im0; %figure, imagesc(im0); colorbar,

            im0(npix+1:end-npix,npix+1:end-npix) = im0(npix+1:end-npix,npix+1:end-npix)./filter2(nhood,im0,'valid');
            %
            im0(1:npix        ,:             )=0;
            im0(end-npix+1:end,:             )=0;
            im0(:             ,1:npix        )=0;
            im0(:             ,end-npix+1:end)=0;

imm = im; %figure, imagesc(im); colorbar,

            im(npix+1:end-npix,npix+1:end-npix) = im(npix+1:end-npix,npix+1:end-npix)./filter2(nhood,im,'valid');
            %
            im(1:npix        ,:             )=0;
            im(end-npix+1:end,:             )=0;
            im(:             ,1:npix        )=0;
            im(:             ,end-npix+1:end)=0;
            %
        end
        %
        [H W] = size(im);
        [x,y] = meshgrid(1:W,1:H);
        %
	% pause
        imfft  = fft2(im (npix+1:end-npix,npix+1:end-npix)-mean(mean(im (npix+1:end-npix,npix+1:end-npix))));
        im0fft = fft2(im0(npix+1:end-npix,npix+1:end-npix)-mean(mean(im0(npix+1:end-npix,npix+1:end-npix))));
        ccc=ifftshift(real(ifft2(imfft.*conj(im0fft))));
        %      ccc = mqd(im,im0,H,W);
        %
        [peak1,peak2,pixi,pixj] = find_displacement(ccc,s2nm); disp(sprintf('peak1=%d ; peak2=%d ; pixi=%d ; pixj=%d',peak1,peak2,pixi,pixj))
%save cccmatrix.mat imm0 imm im0 im H W imfft im0fft ccc, 
        %
        [disy,disx,s2n] = sub_pixel_velocity(ccc,pixi,pixj,peak1,peak2,0,1,H/2-npix+2,W/2-npix+2,isub,iord,subpR,subpr);
        %
        disx = disx-W/2-1+npix,,
        disy = disy-H/2-1+npix,,
        im = interp2(x,y,im,x+disx,y+disy,'spline',0);
        %
        im = real((im));
        im0 = real((im0));
	%
%        im = real(log(im));
%        im0 = real(log(im0));
        %
        A0=max(max(im0(npix+1:end-npix,npix+1:end-npix)));
        B0=min(min(im0(npix+1:end-npix,npix+1:end-npix)));
        %
        A=max(max(im(npix+1:end-npix,npix+1:end-npix)));
        B=min(min(im(npix+1:end-npix,npix+1:end-npix)));

        im0 = ((im0)-B0)/(A0 - B0);
        im  = ((im )-B )/(A  - B );

        j0 = find(im0~=im0);
        im0(j0) = 0;
        j0 = find(im~=im);
        im(j0) = 0;

        im0(1:npix      ,:             )=0;
        im0(end-npix+1:end,:             )=0;
        im0(:             ,1:npix      )=0;
        im0(:             ,end-npix+1:end)=0;

        %
        im(1:npix      ,:             )=0;
        im(end-npix+1:end,:             )=0;
        im(:             ,1:npix      )=0;
        im(:             ,end-npix+1:end)=0;

        %      imagesc(im)
        %      pause
        %      imagesc(im0)
        %      pause

    end

    function [im1,im2,im3,im4,im0] = ...
            im_interp_correc_quad(im,im0,s2nm,isub,iord,subpR,subpr,npix)

        [H W] = size(im);
        if npix>0

            se=strel('disk',npix+1);nhood=double(getnhood(se));

            im0(npix+1:end-npix,npix+1:end-npix) = im0(npix+1:end-npix,npix+1:end-npix)./filter2(nhood,im0,'valid');
            %
            im0(1:npix      ,:             )=0;
            im0(end-npix+1:end,:             )=0;
            im0(:             ,1:npix      )=0;
            im0(:             ,end-npix+1:end)=0;


            im(npix+1:end-npix,npix+1:end-npix) = im(npix+1:end-npix,npix+1:end-npix)./filter2(nhood,im,'valid');
            %
            im(1:npix      ,:             )=0;
            im(end-npix+1:end,:             )=0;
            im(:             ,1:npix      )=0;
            im(:             ,end-npix+1:end)=0;

            %

            [x,y] = meshgrid(1:W,1:H);

            %

            im0q1 = im0(1:H/2    ,1:W/2    );
            im0q2 = im0(H/2+1:end,1:W/2    );
            im0q3 = im0(1:H/2    ,W/2+1:end);
            im0q4 = im0(H/2+1:end,W/2+1:end);

            imq1  = im (1:H/2    ,1:W/2    );
            imq2  = im (H/2+1:end,1:W/2    );
            imq3  = im (1:H/2    ,W/2+1:end);
            imq4  = im (H/2+1:end,W/2+1:end);

            %

            imq1fft  = fft2(imq1 (npix+1:end-npix,npix+1:end-npix)-mean(mean(imq1 (npix+1:end-npix,npix+1:end-npix))));
            im0q1fft = fft2(im0q1(npix+1:end-npix,npix+1:end-npix)-mean(mean(im0q1(npix+1:end-npix,npix+1:end-npix))));
                                                                                                         
            imq2fft  = fft2(imq2 (npix+1:end-npix,npix+1:end-npix)-mean(mean(imq2 (npix+1:end-npix,npix+1:end-npix))));
            im0q2fft = fft2(im0q2(npix+1:end-npix,npix+1:end-npix)-mean(mean(im0q2(npix+1:end-npix,npix+1:end-npix))));
                                                                                                         
            imq3fft  = fft2(imq3 (npix+1:end-npix,npix+1:end-npix)-mean(mean(imq3 (npix+1:end-npix,npix+1:end-npix))));
            im0q3fft = fft2(im0q3(npix+1:end-npix,npix+1:end-npix)-mean(mean(im0q3(npix+1:end-npix,npix+1:end-npix))));
                                                                                                         
            imq4fft  = fft2(imq4 (npix+1:end-npix,npix+1:end-npix)-mean(mean(imq4 (npix+1:end-npix,npix+1:end-npix))));
            im0q4fft = fft2(im0q4(npix+1:end-npix,npix+1:end-npix)-mean(mean(im0q4(npix+1:end-npix,npix+1:end-npix))));
            %

            cccq1=ifftshift(real(ifft2(imq1fft.*conj(im0q1fft))));
            cccq2=ifftshift(real(ifft2(imq2fft.*conj(im0q2fft))));
            cccq3=ifftshift(real(ifft2(imq3fft.*conj(im0q3fft))));
            cccq4=ifftshift(real(ifft2(imq4fft.*conj(im0q4fft))));

            %
            [peak1q1,peak2q1,pixiq1,pixjq1] = find_displacement(cccq1,s2nm);
            [peak1q2,peak2q2,pixiq2,pixjq2] = find_displacement(cccq2,s2nm);
            [peak1q3,peak2q3,pixiq3,pixjq3] = find_displacement(cccq3,s2nm);
            [peak1q4,peak2q4,pixiq4,pixjq4] = find_displacement(cccq4,s2nm);
            %
            [disyq1,disxq1,s2n] = sub_pixel_velocity(cccq1,pixiq1,pixjq1,peak1q1,peak2q1,0,1,H/4-npix+2,W/4-npix+2,isub,iord,subpR,subpr);
            [disyq2,disxq2,s2n] = sub_pixel_velocity(cccq2,pixiq2,pixjq2,peak1q2,peak2q2,0,1,H/4-npix+2,W/4-npix+2,isub,iord,subpR,subpr);
            [disyq3,disxq3,s2n] = sub_pixel_velocity(cccq3,pixiq3,pixjq3,peak1q3,peak2q3,0,1,H/4-npix+2,W/4-npix+2,isub,iord,subpR,subpr);
            [disyq4,disxq4,s2n] = sub_pixel_velocity(cccq4,pixiq4,pixjq4,peak1q4,peak2q4,0,1,H/4-npix+2,W/4-npix+2,isub,iord,subpR,subpr);
            %
            disxq1 = disxq1-W/4-1+npix,,
            disyq1 = disyq1-H/4-1+npix,,
            disxq2 = disxq2-W/4-1+npix,,
            disyq2 = disyq2-H/4-1+npix,,
            disxq3 = disxq3-W/4-1+npix,,
            disyq3 = disyq3-H/4-1+npix,,
            disxq4 = disxq4-W/4-1+npix,,
            disyq4 = disyq4-H/4-1+npix,,
            %
            im1 = interp2(x,y,im,x+disxq1,y+disyq1,'spline',0);
            im2 = interp2(x,y,im,x+disxq2,y+disyq2,'spline',0);
            im3 = interp2(x,y,im,x+disxq3,y+disyq3,'spline',0);
            im4 = interp2(x,y,im,x+disxq4,y+disyq4,'spline',0);
            %
            %      im1 = real(log(im1));
            %      im2 = real(log(im2));
            %      im3 = real(log(im3));
            %      im4 = real(log(im4));
            %      im0 = real(log(im0));
            %
            A=max(max(im0(npix+1:end-npix,npix+1:end-npix)));
            B=min(min(im0(npix+1:end-npix,npix+1:end-npix)));
            %
            im0 = ((im0)-B)/(A - B);
            im1 = ((im1)-B)/(A - B);
            im2 = ((im2)-B)/(A - B);
            im3 = ((im3)-B)/(A - B);
            im4 = ((im4)-B)/(A - B);

            jw = find(im0~=im0);
            im0(jw) = 0;
            jw = find(im1~=im1);
            im1(jw) = 0;

            im0(1:npix      ,:             )=0;
            im0(end-npix+1:end,:             )=0;
            im0(:             ,1:npix      )=0;
            im0(:             ,end-npix+1:end)=0;

            im1(1:npix      ,:             )=0;
            im1(end-npix+1:end,:             )=0;
            im1(:             ,1:npix      )=0;
            im1(:             ,end-npix+1:end)=0;

            im2(1:npix      ,:             )=0;
            im2(end-npix+1:end,:             )=0;
            im2(:             ,1:npix      )=0;
            im2(:             ,end-npix+1:end)=0;

            im3(1:npix      ,:             )=0;
            im3(end-npix+1:end,:             )=0;
            im3(:             ,1:npix      )=0;
            im3(:             ,end-npix+1:end)=0;

            im4(1:npix      ,:             )=0;
            im4(end-npix+1:end,:             )=0;
            im4(:             ,1:npix      )=0;
            im4(:             ,end-npix+1:end)=0;

        else
            im1 = im;
            im2 = im;
            im3 = im;
            im4 = im;
        end

        %   imagesc(im)
        %   pause
        %   imagesc(im0)
        %   pause
        %
    end


    function [c] = cross_correlate(a2,b2,Nfft)
        % CROSS_CORRELATE - calculates the cross-correlation
        % matrix of two interrogation areas: 'a2' and 'b2' using
        % IFFT(FFT.*Conj(FFT)) method.
        % Modified version of 'xcorrf.m' function from ftp.mathworks.com
        % site.
        % Authors: Alex Liberzon & Roi Gurka
        %

        c = zeros(Nfft,Nfft);
        % Remove Mean Intensity from each image
        a2 = a2 - mean2(a2);
        b2 = b2 - mean2(b2);
        % FFT of both:
        ffta=fft2(a2,Nfft,Nfft);
        fftb=fft2(b2,Nfft,Nfft);
        % Real part of an Inverse FFT of a conjugate multiplication:
        c = ifftshift(real(ifft2(ffta.*conj(fftb))));
    end

    function [peak1,peak2,pixi,pixj] = find_displacement(c,s2nm)
        % FIND_DISPLACEMENT - Finds the highest peak in cross-correlation
        % matrix and the second peak (or mean value) for signal-to-noise
        % ratio calculation.
        % Inputs:
        %         c - cross-correlation matrix
        %         s2nm - method (1 or 2) of S2N ratio calculation
        % Outputs:
        %         peak1 = highest peak
        %         peak2 = second highest peak (or mean value)
        %         pixi,pixj = row,column indeces of the peak1
        %
        % Authors: Alex Liberzon & Roi Gurka
        % Date: 20-Jul-99
        % Last modified:
        %

        % Find your majour peak = mean pixel displacement between
        % two interrogation areas:

        [Nfft,junk] = size(c);

        peak1 = max(c(:));
        [pixi,pixj]=find(c==peak1);

        % Temproraly matrix without the maximum peak:
        tmp = c;
        tmp(pixi,pixj) = 0;
        % If the peak is found on the border, we should not accept it:
        if pixi==1 | pixj==1 | pixi==Nfft | pixj==Nfft
            peak2 = peak1; % we'll not accept this peak later, by means of S2N
        else
            % Look for the Signal-To-Noise ratio by
            % 1. Peak detectability method: First-to-second peak ratio
            % 2. Peak-to-mean ratio - Signal-to-noise estimation

            if s2nm == 1		% First-to-second peak ratio
                % Remove 3x3 pixels neighbourhood around the peak
                tmp(pixi-1:pixi+1,pixj-1:pixj+1) = NaN;
                % Look for the second highest peak
                peak2 = max(tmp(:));
                [x2,y2] = find(tmp==peak2);
                tmp(x2,y2) = NaN;
                % Only if second peak is within the borders
                if x2 > 1 & y2 > 1 & x2 < Nfft & y2 < Nfft

                    % Look for the clear (global) peak, not for a local maximum:
                    while peak2 < max(max(c(x2-1:x2+1,y2-1:y2+1)))
                        peak2 = max(tmp(:));
                        [x2,y2] = find(tmp==peak2);
                        if x2 == 1 | y2==1 | x2 == Nfft | y2 == Nfft
                            peak2 = peak1;	% will throw this one out later
                            break;
                        end
                        tmp(x2,y2) = NaN;
                    end		% end of while
                else			% second peak on the border means "second peak doesn't exist"
                    peak2 = peak1;
                end    % if x2 >1 ......end
                % PEAK-TO-MEAN VALUE RATIO:
            elseif s2nm == 2
                peak2 = mean2(abs(tmp));
            end		% end of second peak search, both methods.
        end				% end of if highest peak on the border
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    EXTERNAL FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [A1,A2,A3,B1,B2,B3] = read_3_images(imagea1,imagea2,imagea3,imageb1,imageb2,imageb3,crop_vector,itt,spc,ext)
        % READ_PAIR_OF_IMAGES - reads 4 images (image1,... image4) as tif files
        % and crops them according to 'crop_vector'
        % Inputs:
        %         image1,... image4 - tif file names (string)
        %         crop_vector - 4 x 1 vector of follwoing values:
        %         [left,top,right,bottom] - each value is a number of lines
        %                                   of interrogation areas (ITTxITT pixels)
        %                                   which should be removed before the analysis.
        %         itt - interrogation area size in pixels
        %         spc - grid spacing (overlapping) size in pixels
        %
        % juanc: adapted from read_pair_of_images, 06/14/07


        A1 = imread(imagea1,ext);
        A2 = imread(imagea2,ext);
        A3 = imread(imagea3,ext);

        B1 =  imread(imageb1,ext);
        B2 =  imread(imageb2,ext);
        B3 =  imread(imageb3,ext);

        A1 = double(A1(:,:,1));
        A2 = double(A2(:,:,1));
        A3 = double(A3(:,:,1));

        B1  = double(B1(:,:,1));
        B2  = double(B2(:,:,1));
        B3  = double(B3(:,:,1));

        %         A = rgb2gray(imread(image1,'tif'));
        %         B = rgb2gray(imread(image2,'tif'));
        [sxa1,sya1]=size(A1);
        [sxa2,sya2]=size(A2);
        [sxa3,sya3]=size(A3);

        [sxb1,syb1]=size(B1);
        [sxb2,syb2]=size(B2);
        [sxb3,syb3]=size(B3);

        % A & B matrices HAVE to be of the same size, we take smallest:
        sx = min(sxa1,sxa2);
        sx = min(sx  ,sxa3);
        sx = min(sx  ,sxb1);
        sx = min(sx  ,sxb2);
        sx = min(sx  ,sxb3);

        sy = min(sya1,sya2);
        sy = min(sy  ,sya3);
        sy = min(sy  ,syb1);
        sy = min(sy  ,syb2);
        sy = min(sy  ,syb3);

        % Crop the images to the desired size and
        % cut the last couple of pixels, so we'll get the
        % integer number of interrogation areas
        %
        %       ---- t ---
        %      |          |
        %      |          |
        %      l          r
        %      |          |
        %      |          |
        %       --- b ----
        %
        %
        l = crop_vector(1); % left side of the image
        t = crop_vector(2); % top side of the image
        r = crop_vector(3); % right side of the image
        b = crop_vector(4); % bottom of the image

        A1 = A1(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);
        A2 = A2(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);
        A3 = A3(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);

        B1 = B1(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);
        B2 = B2(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);
        B3 = B3(1+t*itt:spc*floor(sx/spc)-b*itt,1+l*itt:spc*floor(sy/spc)-r*itt);
    end

    function [peakx,peaky,s2n] = sub_pixel_velocity(c,pixi,pixj,peak1,peak2,s2nl,sclt,ittx,itty,isub,iord,subpR,subpr)
        % SUB_PIXEL_VELOCITY - Calculates Signal-To-Noise Ratio, fits Gaussian
        % bell, find sub-pixel displacement and scales it to the real velocity
        % according the the time interval and real-world-to-image-scale.
        %
        % Authors: Alex Liberzon & Roi Gurka
        % Date: Jul-20-99
        % Last Modified:

        % If peak2 equals to zero, it means that nothing was found,
        % and we'll divide by zero:

        % isub = 0 ==> Quick and Dirty Gaussian (original)
        % isub = 1 ==> Log poly order iord
        % isub = 2 ==> poly order iord

	if pixi ==2*ittx | pixj==2*itty; 
           peakx = ittx;
           peaky = itty;
           s2n = Inf;
	   return; 
	end
        if max(abs(c(:))) < 1e-3;
            peakx = ittx;
            peaky = itty;
            s2n = Inf;
            return
        end

        if ~peak2
            s2n = Inf;		% Just to protect from zero dividing.
        else
            s2n = peak1/peak2;
        end

        % If Signal-To-Noise ratio is lower than the limit, "mark" it:
        if s2n < s2nl
            peakx = ittx;
            peaky = itty;
	    return
        else            % otherwise, calculate the velocity

            if isub ==0;
                % Sub-pixel displacement definition by means of
                % Gaussian bell.

                f0 = log(c(pixi,pixj));
                f1 = log(c(pixi-1,pixj));
                f2 = log(c(pixi+1,pixj));
                peakx = pixi+ (f1-f2)/(2*f1-4*f0+2*f2);
                f0 = log(c(pixi,pixj));
                f1 = log(c(pixi,pixj-1));
                f2 = log(c(pixi,pixj+1));
                peaky = pixj+ (f1-f2)/(2*f1-4*f0+2*f2);

                if ~isreal(peakx) | ~isreal(peaky)
                    peakx = ittx;
                    peaky = itty;
                end
            else
                w=subpR;
                rad=subpr;
                x0 = [-w:w];
                y0 = [-w:w];
                vecx0 = pixi + x0;
                vecy0 = pixj + y0;
                if isub == 1; z = log(abs(c)); end
                if isub == 2; z = c; end
                kernel = zeros(2*w+1);
                for i=1:2*w+1;
                    kernel(i,:) =  exp(-((i-w-1)^2+[-w:w].^2)/(rad^2));
                end
                p = polyfitweighted2(x0,y0,z(vecx0,vecy0),iord,kernel);
                %
                options = optimset('display','none','TolX',1e-8,'TolFun',1e-8);
                v = fminsearch(@(x) wrapfun(x,p),[0 0],options);
                %               v = fminsearch(@(x) wrapint2(x,x0,y0,c),[0 0],options);
                peakx = v(2) + pixi;
                peaky = v(1) + pixj;
                %
                %x = [-w:w];
                %y = [-w:w];
                %z = polyval2(p,x,y);
                %subplot(2,1,1)
                %imagesc(c(vecx0,vecy0));
                %cax=caxis;
                %hold on
                %plot(v(1)+w+1,v(2)+w+1,'o');
                %colorbar
                %hold off
                %axis tight
                %subplot(2,1,2)
                %imagesc(z);
                %caxis(cax);
                %colorbar
                %hold on
                %plot(v(1)+w+1,v(2)+w+1,'o');
                %hold off
                %axis tight
                %pause(1)

            end

        end
    end


    function [filenames,amount,ext] = ReadImDir(directory,name,ext)

        % graphical extensions
        knownExtensions = {'bmp','jpg','tif','tiff','jpeg'};

        if nargin == 1 | isempty(ext)
            for i = 1:length(knownExtensions)
                direc = dir([directory,filesep,name,'*.',knownExtensions{i}]);
                if ~isempty(direc), ext = knownExtensions{i}; break, end
            end
        end
        direc = dir([directory,filesep,name,'*.',ext]);
        amount = length(direc);
        filenames={};
        [filenames{1:amount,1}] = deal(direc.name);
        dates={};
        [dates{1:amount,1}] = deal(direc.date);
        [dates,isordates] = sortrows(dates);
        filenames = filenames(isordates);
        %         filebase = filename(1:regexpi(filename,'\d','once')-1);

        %         if ~isempty(findstr(filenames{1},'_b')) % case 1 _b, _c
        % %             amount = max(str2num(filenames{1,end-8:end-6}));
        %             amount = str2num(filenames{1}(regexp(filenames{1},'\s*\d')));
        %             filebase = filenames{1}(1:regexp(filenames{1},'_','once'));
        %         else % sequential numbering, UTAH case
        %             amount = length(direc)-1; % max(str2num(filenames(:,end-9:end-4)));
        %             %             filebase = filenames(1,regexp(filenames(1,:),'\w\d*[.]')+1:regexp(filenames(1,:),'[.]')-1);
        %             filebase = filenames(1,1:regexp(filenames(1,:),'\w\d*[.]'));
        %
        %         end
    end

    function [] = quiverm(x,varargin)
        % QUIVERM - plots quiver plot of matrix,
        % assuming first column as X, second as Y
        % third as U, and forth as V.
        %
        % QUIVERM(A,'r') - plots quiver plot
        % of matrix A.
        % Used by QUIVERTXT function.
        %
        %
        %
        % Author: Alex Liberzon
        %

        if isstr(x)
            x = eval(x);
        end
        quiver(x(:,1),x(:,2),x(:,3),x(:,4),3,varargin{:});
    end

    function cmqd = mqd(a,b,nx,ny);

        xmargin = 5;
        ymargin = 5;

        cmqd = zeros(nx,ny);

        for iM=-xmargin:xmargin;
            vecxa = [xmargin+1:nx-xmargin];
            vecxb = vecxa - iM;
            for jM=-ymargin:ymargin;
                vecya = [ymargin+1:ny-ymargin];
                vecyb = vecya - jM;
                cmqd(iM+nx/2+1,jM+ny/2+1) = sum(sum((a(vecxa,vecya)-b(vecxb,vecyb)).^2));
            end;
        end;
        cmax = max(max(cmqd));
        jjj = find(cmqd);
        cmqd(jjj) = cmax - cmqd(jjj);

    end
end
