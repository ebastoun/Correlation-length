% readu_correlation_length_WT - Script to read txt files with cellular
% displacement vectors strored, plot them on top of the phase contrast image of cells (or image of their 
% Hoechst-stained nuclei) and calculate the correlation length of movement between neighbouring 
% cells same as in Angenlini et al. PRL 104, 168104 (2010)
%
% 1. Read the phase contrast images and plot on top of them the
% displacement vectors that cells undergo (frames taken every 10 min)
% 2. Construct kymographs of the average with respect to the wound
% displacement vectors, uy as a function of time, t (min) and vertical position, y (μm)  
% 3. Calculation of mean uy along a distance of 100 μm away from the wound 
% and calculated immediately after ablation
% 4. Storage of the results as a .mat file
% Last modified: F. Serrano-Alcalde and E. Bastounis: 2020-06-10

clear all; close all; clc
% First and last frames
kin=2; kfin=40;
% Time between frame intervals
delta=10;
% Calibration factor (um/pixel)
fcal=0.443;

all_velocities=[];
kymograph=zeros(26,kfin);% 26 because of R_max=25 and R_step=1
% Note: below commented are the parameters we used to run the urapiv.m code to
% compare subsequent frames. We used windows of 48x48 and 50% overlap (24
% pixels)
%  for k=2:40
%  urapiv('/Users/effiebastounis/Documents/MDCKs/My_paper/final_revision_2020_11_21/Cell_star_protocols_2021_01_26/scripts_correlation_length/Pos0/','/Users/effiebastounis/Documents/MDCKs/My_paper/final_revision_2020_11_21/Cell_star_protocols_2021_01_26/scripts_correlation_length/urapiv/Pos0/','nuclei',1,48,24,1,2,0,1,0,0,0*[1 1 1 1],0,2,3,3,8,[k-1],[k],1);
%  end  
overlap=24;

for k= kin:kfin
   
fil1   = '/Users/effiebastounis/Documents/MDCKs/My_paper/final_revision_2020_11_21/Cell_star_protocols_2021_01_26/scripts_correlation_length/urapiv/Pos0/nuclei'; %txt files from PIV
dirname='/Users/effiebastounis/Documents/MDCKs/My_paper/final_revision_2020_11_21/Cell_star_protocols_2021_01_26/scripts_correlation_length/Pos0/nuclei';% image of nuclei
img = double(imread([dirname sprintf('%3.3d.tif', k)]));

figure;imagesc(img);colormap gray;axis image;hold on;

filename = ([fil1 sprintf('%3.3d.txt', k ) ]);
outl = 10;
umax=100;
s6nl= 0; %normally 0
vec = load(filename); 
   n1 = min(find(diff(vec(:,1))<0));
   n2 = length(vec)/n1;

   vec=reshape(vec,n1,n2,5);
   
   x    = vec(:,:,1);
   y    = vec(:,:,2);
   u    = vec(:,:,3);
   v    = vec(:,:,4);
   s6n = vec(:,:,5);

   if s6nl>0

      % signal to noise check
      
      s6n_low = find(s6n<s6nl);
      
      w = 3; rad = 2;
      kernel = zeros(2*w+1);
      for i=1:2*w+1;
          kernel(i,:) =  exp(-((i-w-1)^2+[-w:w].^2)/rad^2);
      end
      kernel(w+1,w+1)=0;
      kernel = kernel/sum(sum(kernel));
   
      u(s6n_low) = 0;
      v(s6n_low) = 0;
      
      tmpv = (conv2(v,kernel,'same'));
      tmpu = (conv2(u,kernel,'same'));
      
      % Let's throw the outlayers out:

      u(s6n_low) = tmpu(s6n_low); 
      v(s6n_low) = tmpv(s6n_low); 
      u(s6n<s6nl)=NaN;
      v(s6n<s6nl)=NaN;

     uu=sqrt(u.^2+v.^2);
     u(uu>umax)=NaN;
     v(uu>umax)=NaN;
  end 
   
  if outl>0

      % Adaptive Local Median filtering
   
      w = 2; rad = 1;
      kernel = zeros(2*w+1);
      for i=1:2*w+1;
          kernel(i,:) =  exp(-((i-w-1)^2+[-w:w].^2)/rad^2);
      end
      kernel(w+1,w+1)=0;
      kernel = kernel/sum(sum(kernel));
   
      tmpv = (conv2(v,kernel,'same'));
      tmpu = (conv2(u,kernel,'same'));
   
      lmtv_p = mean(mean(tmpv(2:end-1,2:end-1))) + ...
             outl*std(reshape(tmpv(2:end-1,2:end-1),(n1-2)*(n2-2),1));
      lmtv_m = mean(mean(tmpv(2:end-1,2:end-1))) - ...
             outl*std(reshape(tmpv(2:end-1,2:end-1),(n1-2)*(n2-2),1));
      lmtu_p = mean(mean(u(2:end-1,2:end-1))) + ...
             outl*std(reshape(u(2:end-1,2:end-1),(n1-2)*(n2-2),1));
      lmtu_m = mean(mean(u(2:end-1,2:end-1))) - ...
             outl*std(reshape(u(2:end-1,2:end-1),(n1-2)*(n2-2),1));

      u_out_p = find(u>lmtu_p);
      u_out_m = find(u<lmtu_m);
      v_out_p = find(v>lmtv_p);
      v_out_m = find(v<lmtv_m);
   
      % Let's throw the outlayers out:

      u(u_out_m) = tmpu(u_out_m); 
      u(v_out_m) = tmpu(v_out_m); 
      v(u_out_m) = tmpv(u_out_m); 
      v(v_out_m) = tmpv(v_out_m); 

      u(u_out_p) = tmpu(u_out_p); 
      u(v_out_p) = tmpu(v_out_p); 
      v(u_out_p) = tmpv(u_out_p); 
      v(v_out_p) = tmpv(v_out_p); 

  end  
  
qq=1;
hold on
quiver(x(1:qq:end,1:qq:end),y(1:qq:end,1:qq:end),u(1:qq:end,1:qq:end)*10*fcal,v(1:qq:end,1:qq:end)*10*fcal,'AutoScale','off','color',[0 1 0]);%htt=text(110,860,'20 \mum','fontsize',24); 

u=u';v=v';

u(isnan(u))=0;
v(isnan(v))=0;


%% Correlacion Angelini et al. (v1)
%Subtract of mean
mean_u=mean(mean(u));
mean_v=mean(mean(v));
u=u-mean_u;
v=v-mean_v;

%Generate matrix with distances
Generar_distancias_i=1:1:length(u);
Generar_distancias_i=repmat(Generar_distancias_i,length(Generar_distancias_i),1);
Generar_distancias_j=(1:1:length(u))';
Generar_distancias_j=repmat(Generar_distancias_j,1,length(Generar_distancias_j));

%Initilization of matrices
P_escalar_superior=zeros(length(u));
P_escalar_inferior=zeros(length(u));
P_escalar=[];
numero_R=1;
R_ant=0;
bordes=0; % Changing this value we can remove the correlation of the edges. (pixels)
R_step=1; % Step of the radius discretization (pixels)
R_max=25; % Maximum radius (pixels). This corresponds to distance R_max*overal*fcal
for R=0:R_step:R_max
    for i=1+bordes:length(u)-bordes
        for j=1+bordes:length(u)-bordes
            Distancias=sqrt((i-Generar_distancias_i).^2+(j-Generar_distancias_j).^2);
            [Posiciones_y,Posiciones_x]=find(Distancias<=R & Distancias>=R_ant);
            vector_0=[u(i,j) v(i,j)];
            for escalar=1:length(Posiciones_x)
                vector_Posicion=[u(Posiciones_x(escalar),Posiciones_y(escalar)) v(Posiciones_x(escalar),Posiciones_y(escalar))]';
                P_escalar(i,j,escalar)=vector_0*vector_Posicion;
            end
            P_escalar_superior(i,j)=mean(P_escalar(i,j,find(P_escalar(i,j,:)~=0 & P_escalar(i,j,:)~=NaN)));
            P_escalar_inferior(i,j)=vector_0*vector_0';
        end
    end
    Coef(numero_R, k)=sum(sum(P_escalar_superior))/sum(sum(P_escalar_inferior)); %Correlation for each R iteration for each image(k)
    valor_Radio(numero_R,k)=R; % Radius for each iteration for each image(k)(the same for all images)
    numero_R=numero_R+1;
    
    R_ant=R;
    
    %%%%
    %Image of the correlation for a given radius (R)
    %Be careful (and comment) if R_step and/or k (number of images) is too large!!!
    Coef_2D=P_escalar_superior./P_escalar_inferior; %Matrix of the individual correlation of each pixel for the radius R
    % Plot the 2D map of the correlation coefficinet for each R considered
%     figure
%     imagesc([1:length(Coef_2D)]*24*fcal,[1:length(Coef_2D)]*24*fcal, Coef_2D); caxis([-1 1]);colorbar;axis image
%     set(gca,'FontSize',18,'DefaultAxesFontName', 'Arial');title(R);xlabel('um');ylabel('um');
    %%%%
end

valor_Radio(:,k)=valor_Radio(:,k)*overlap*fcal; %Change of units from pixels to um
time(k)=k*delta; %time in 'delta' units

%% Length from equation
%Fit equation
eqexpon=fittype('exp(-R/R0)', 'independent',{'R'}, 'coefficients',{'R0'});
myfit=fit(valor_Radio(:,k),Coef(:,k),eqexpon,'startpoint', [50], 'maxiter',[20000]);
R0_equation(k)=myfit.R0 %Value of R0 from equation

fig=figure;
plot(myfit,valor_Radio(:,k),Coef(:,k))
xlabel('Radius (um)','FontSize',10,'FontName','Arial','FontWeight','bold')
ylabel('Correlation','FontSize',10,'FontName','Arial','FontWeight','bold')
name=strcat('Time_WT=',int2str(k));
title(name)
print(fig,'-dtiff',name);
close(fig)
kymograph(:,k)=Coef(:,k);

%% Length from minimum value
[min_cor, position]=min(Coef(:,k));
R0_minimum(k)=valor_Radio(position,k);
close all;
end 

% Plot of evolution of mean correlation length in each frame
fig2=figure
hold on
plot(time,R0_equation)
plot(time,R0_minimum)
xlabel('Time (min)','FontSize',10,'FontName','Arial','FontWeight','bold')
ylabel('Correlation length (um)','FontSize',10,'FontName','Arial','FontWeight','bold')
title('Correlation')
legend('equation','minimum value')
hold off
print(fig2,'-dtiff','Correlation_WT');

% Kymograph where you can see how the correlation coefficient changes as a
% function of time and space
fig3=figure
imagesc(time/60,[1:length(kymograph)]*24*fcal, kymograph); colormap redblue;colorbar;
caxis([-0.5 1])
xlabel('Time (h)','FontSize',10,'FontName','Arial','FontWeight','bold')
ylabel('Correlation coefficient','FontSize',10,'FontName','Arial','FontWeight','bold')
title('Correlation kymograph')
set(gca,'FontSize',18,'DefaultAxesFontName', 'Arial');
print(fig3,'-dtiff','Kymograph_WT');
