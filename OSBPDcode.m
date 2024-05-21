%   This code is written for ordinary state-based peridynamic theory that solves the problem given in Section 5.1.
%
%
%   Dorduncu M, Ren H, Zhuang X, Silling S, Madenci E, Rabczuk T. "A review of peridynamic theory and nonlocal operators 
%   along with their computer implementations." Computers & Structures 2024;299:107395. https://doi.org/10.1016/j.compstruc.2024.107395.
%

clc;
clear;
close all;
echo off;
global icount
global ndivx ndivy totnode
global nbdx nbdy
global coord dx area vol delta
global ntype numfam family
global fncstd fncstm

icount=0;
%length: Length of the plate
length = 100.0;
%width: Width of the plate
width = 100.0;
% Specify the number of intervals in x direction
ndivx = 80;
% Specify the number of intervals in y direction
ndivy = 80;
% Specify the number of ghost points in x direction
nbdx = 0;
% Specify the number of ghost points in y direction
nbdy = 3;
% Compute total number of PD points
totnode= (ndivx+ 2 * nbdx) * (ndivy+ 2 * nbdy);
% Compute interval size
dx = length/ndivx;
% Cross-sectional area
area = dx * dx;
% Thickness
thick = 1;
% Volume of a material point
vol = area * thick;
% Material point radius
radij = dx / 2.0;
% Compute horizon size in each dimension
delta = (nbdy + 0.015d0) * dx;
% Initialize coordinate array
coord=zeros(totnode,2);
% Node type
ntype = zeros(totnode,1);
% Initialize family array of a material point.
% Number of family nodes
numfam = zeros(totnode,1);
% Force normalization constants
fncstd = zeros(totnode,2);
fncstm = zeros(totnode,2);
% Peridynamic force
pforce = zeros(totnode,2);
% Peridynamic force (old)
pforceold = zeros(totnode,2);
% Dilatation
dilm = zeros(totnode,1);
% Body force
bforce = zeros(totnode,2);
% Displacement
disp = zeros(totnode,2);
% Velocity
vel = zeros(totnode,2);
% Velocity half time increment
velhalf = zeros(totnode,2);
% Velocity after half time increment (old)
velhalfold = zeros(totnode,2);
% Acceleration
acc = zeros(totnode,2);
% Mass vector
massvec = zeros(totnode,2);
%dmg: Damage parameter
dmg = zeros(totnode,1);
% Generate peridynamic grid and store the coordinates to coord array.
GeneratePoints
% Generate family of each material point
GenerateFamily
% Apply slit for a pre-existing crack (please go to the subroutine)
ApplySlit
% Specify the density
rho = 1.0;
% Specify the elastic modulus
emod = 65.8e3;
% Specify the Poisson's ratio
pratio = 0.38;
% Shear modulus
smod = emod / (2.0 * (1.0 + pratio));
% Bulk modulus
bmod = emod / (2.0 * (1.0 - pratio));
% PD parameters
alpha= 0.5*(bmod-2.0*smod);
% Bond constant
bcd = 2.0/(pi*thick*delta^3);
% Bond constant
bcs = 6.0*smod/(pi*thick*delta^4);
% Critical energy release rate
gc = 4/1000;
%scr0: Critical stretch
% scr0 = sqrt(4.0d0 * pi * gc / (9.0d0 * emod * delta));
scr0 = sqrt(gc / (6 * smod / pi + 16 * (bmod - 2*smod) / (9*pi^2))/delta);
%Dilatation for surface correction factor calculation
dilm1 = 0.001;
dilm2 = 0.001;
%sedload1: Strain energy density of a material point under uniaxial tension in x direction
sedloadm1 = 1.0 / (2.0*(1.0-pratio^2)) * emod * (dilm1)^2 - alpha * (dilm1)^2;
sedloadm2 = 1.0 / (2.0*(1.0-pratio^2)) * emod * (dilm2)^2 - alpha * (dilm2)^2;
%dt: Time interval
dt = 1.0;
%dtload: Applied displacement constraint
u0 = 3e-6;
%totime: Total time
nt = 2000;
totime = nt * dt;

% Reaction force vector
fplane1 = zeros(nt,1);
fplane1t = zeros(nt,1);
fplane2 = zeros(nt,1);
fplane2t = zeros(nt,1);

% Compute the surface correction factors in x and y directions for each material
SCorrection(sedloadm1, sedloadm2, bcd, bcs, radij)

% Stable mass vector computation
 massvec(:,:) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * 4 * ...
                bcs * delta / dx  ;
           

% Initialize the displacements and velocities
vel(:,:) = 0.0;
disp(:,:) = 0.0;

% Start time integration
for tt = 1:nt
    ctime = tt * dt;
    fprintf('Itr = %8.0d  Time(s) = %10.8f \n',tt,ctime);
    
    icount = 0;
    
    for i = 1:totnode
        
        % Apply the initial conditions at the top and bottom edges of the plate
        
        if (coord(i,2) <= 0.0)
            vel(i,2) = -u0;
            disp(i,2) = (-u0 + 2.0 * u0 * coord(i,2)/length)* tt * dt; 
        elseif (coord(i,2) >= length)
            vel(i,2) = u0;
            disp(i,2) = (-u0 + 2.0 * u0 * coord(i,2)/length) * tt * dt; 
        end
        
        if (ntype(i) == 1)
            icount= icount + 1;
            results(icount,:) = [coord(i,1) coord(i,2) disp(i,1) ...
                                 disp(i,2) dmg(i,1)];
        end
    end
    
    
    for i = 1:totnode
        
        dilm(i) = 0.0d0;
        
        for jmem = 1:numfam(i,1)
            
            j = family(i).node(jmem);
            jflag = family(i).nodeFlag(jmem);
            ix1 = (coord(j,1) - coord(i,1));
            ix2 = (coord(j,2) - coord(i,2));
            nx1 = (coord(j,1)+disp(j,1)-coord(i,1)-disp(i,1));
            nx2 = (coord(j,2)+disp(j,2)-coord(i,2)-disp(i,2));
            
            idist   = sqrt(ix1^2 + ix2^2);
            nlength = sqrt(nx1^2+nx2^2);
            str = (nlength - idist) / idist;
            pplen = nlength * idist;
            
            ppx = nx1 * ix1;
            ppy = nx2 * ix2;
            pptot = 1; %(ppx + ppy) / pplen;
            
            if (idist <= delta-radij)
                fac = 1.0;
            elseif (idist <= delta+radij)
                fac = (delta+radij-idist)/(2.0*radij);
            else
                fac = 0.0;
            end
            
            
            if (abs(coord(j,2) - coord(i,2)) <= 1.0d-10)
                theta = 0.0;
            elseif (abs(ix1) <= 1.0d-10)
                theta = 90.0 * pi / 180.0;
            else
                theta = atan(abs(ix2) / abs(ix1));
            end
            
            scmx = (fncstd(i,1) + fncstd(j,1)) / 2.0;
            scmy = (fncstd(i,2) + fncstd(j,2)) / 2.0;
            scmr = 1.0d0 / ( ( cos(theta)^2.0 / scmx^2.0)...
                + (sin(theta)^2.0 / scmy^2.0d0) );
            scmr = sqrt(scmr);
            
            dilm(i) = dilm(i) + bcd * delta * str * pptot * vol ...
                                * scmr * fac * jflag;
                            
        end
    end
    
    for i = 1:totnode
        dmgpar1 = 0.0;
        dmgpar2 = 0.0;
        pforce(i,1) = 0.0;
        pforce(i,2) = 0.0;
        for jmem = 1:numfam(i,1)
            j = family(i).node(jmem);
            jflag = family(i).nodeFlag(jmem);   
            
            ix1 = (coord(j,1) - coord(i,1));
            ix2 = (coord(j,2) - coord(i,2));
            nx1 = (coord(j,1)+disp(j,1)-coord(i,1)-disp(i,1));
            nx2 = (coord(j,2)+disp(j,2)-coord(i,2)-disp(i,2));
            
            idist   = sqrt(ix1^2 + ix2^2);
            nlength = sqrt(nx1^2+nx2^2);
            str = (nlength - idist) / idist;
            pplen = nlength * idist;
            
            ppx = nx1 * ix1;
            ppy = nx2 * ix2;
            pptot = (ppx + ppy) / pplen;
            
            %Volume correction
            if (idist <= delta-radij)
                fac = 1.0;
            elseif (idist <= delta+radij)
                fac = (delta+radij-idist)/(2.0*radij);
            else
                fac = 0.0;
            end
            
            
            if (abs(coord(j,2) - coord(i,2)) <= 1.0d-10)
                theta = 0.0;
            elseif (abs(ix1) <= 1.0d-10)
                theta = 90.0 * pi / 180.0;
            else
                theta = atan(abs(ix2) / abs(ix1));
            end
            
            scmx = (fncstm(i,1) + fncstm(j,1)) / 2.0;
            scmy = (fncstm(i,2) + fncstm(j,2)) / 2.0;
            scmr = 1.0d0 / ( ( cos(theta)^2.0 / scmx^2.0)...
                + (sin(theta)^2.0 / scmy^2.0d0) );
            scmr = sqrt(scmr);
            
            
            scmx = (fncstd(i,1) + fncstd(j,1)) / 2.0;
            scmy = (fncstd(i,2) + fncstd(j,2)) / 2.0;
            scdr = 1.0d0 / ( ( cos(theta)^2.0 / scmx^2.0)...
                + (sin(theta)^2.0 / scmy^2.0d0) );
            scdr = sqrt(scdr);
            
            
            if (jflag==1)
                %Calculate the peridynamic forces in x and y
                %directions acting on material point i due to material point j
                                             
                dforce1 = (2*bcd* delta *alpha/idist*(dilm(i)+dilm(j))*scdr ...
                    + 4.0 * bcs * delta * str * scmr) * vol * fac * nx1 / nlength;
                dforce2 = (2*bcd* delta *alpha/idist*(dilm(i)+dilm(j))*scdr ...
                    + 4.0 * bcs * delta * str * scmr) * vol * fac * nx2 / nlength;     
                

            else
                dforce1 = 0.0;
                dforce2 = 0.0;
            end
            pforce(i,1) = pforce(i,1) + dforce1;
            pforce(i,2) = pforce(i,2) + dforce2;
            
            if ((coord(j,2) > length/4) && (coord(i,2) < length/4)) 
               fplane2(tt) = fplane2(tt) + dforce2 * vol;  
               fplane2t(tt) = fplane2t(tt) + dforce1 * vol; 
            end     
            
    %       Break the PD bond if its stretch exceeds the critical stretch of the bond.
            if (abs(str) > scr0)
                family(i).nodeFlag(jmem) = 0;
            end
            
            dmgpar1 = dmgpar1 + jflag * vol * fac;
            dmgpar2 = dmgpar2 + vol * fac;
        end
        dmg(i,1) = 1.0 - dmgpar1 / dmgpar2;
    end
    
    %Compute the adaptive dynamic relaxation technique parameters
    
    cn = 0.0;
    cn1 = 0.0;
    cn2 = 0.0;
    for i = 1:totnode
        if (velhalfold(i,1)~=0.0)
            cn1 = cn1 - disp(i,1) * disp(i,1) * (pforce(i,1) / massvec(i,1) ...
                - pforceold(i,1) / massvec(i,1)) / (dt * velhalfold(i,1));
        end
        if (velhalfold(i,2)~=0.0)
            cn1 = cn1 - disp(i,2) * disp(i,2) * (pforce(i,2) / massvec(i,2) ...
                - pforceold(i,2) / massvec(i,2)) / (dt * velhalfold(i,2));
        end
        cn2 = cn2 + disp(i,1) * disp(i,1);
        cn2 = cn2 + disp(i,2) * disp(i,2);
    end
    
    if (cn2~=0.0)
        if ((cn1 / cn2) > 0.0)
            cn = 2.0 * sqrt(cn1 / cn2);
        else
            cn = 0.0;
        end
    else
        cn = 0.0;
    end
    
    if (cn > 2.0)
        cn = 1.9;
    end

    
    % Apply boundary conditions and perform time integration to obtain displacements
    % and velocities.
    
    for i = 1:totnode
        if (tt == 1)
            velhalf(i,1) = 1.0 * dt / massvec(i,1) * (pforce(i,1) + bforce(i,1)) / 2.0;
            velhalf(i,2) = 1.0 * dt / massvec(i,2) * (pforce(i,2) + bforce(i,2)) / 2.0;
        else
            velhalf(i,1) = ((2.0 - cn * dt) * velhalfold(i,1) + ...
                2.0 * dt / massvec(i,1) * (pforce(i,1) + bforce(i,1))) / (2.0 + cn * dt);
            velhalf(i,2) = ((2.0 - cn * dt) * velhalfold(i,2) + ...
                2.0 * dt / massvec(i,2) * (pforce(i,2) + bforce(i,2))) / (2.0 + cn * dt);
        end
        
        vel(i,1) = 0.5 * (velhalfold(i,1) + velhalf(i,1));
        vel(i,2) = 0.5 * (velhalfold(i,2) + velhalf(i,2));
        disp(i,1) = disp(i,1) + velhalf(i,1) * dt;
        disp(i,2) = disp(i,2) + velhalf(i,2) * dt;
        
        velhalfold(i,1) = velhalf(i,1);
        velhalfold(i,2) = velhalf(i,2);
        pforceold(i,1) = pforce(i,1);
        pforceold(i,2) = pforce(i,2);
    end
    
%     drecord(tt,:) = [tt disp(2038,1) disp(2038,2)];
    
    if (tt==1000)
        save -ascii results.dat results
        scatter(coord(:,1),coord(:,2),100,dmg(:,1),'filled')
    end
    
end


% Function to generate the grid of the multidimensional domain
function GeneratePoints( )
global icount;
global ndivx ndivy nbdx nbdy
global coord ntype dx

for i = (1-nbdy):ndivy+(nbdy)
    for j = (1-nbdx):ndivx+(nbdx)
        icount=icount+1;
        coord(icount,1) = (dx / 2.0) + (j-1) * dx;
        coord(icount,2) = (dx / 2.0) + (i-1) * dx;
        
        if (i>0 && j >0)
            if (i <= ndivx && j <= ndivy)
                ntype(icount) = 1;
            end
        end
        
    end
end

end


% Function to generate the grid of the multidimensional domain
function GenerateFamily( )
global coord delta totnode
global numfam family

% Increase the size of this array if family members exceed 10000.
nodefam = zeros(10000,1);

% Determine the material points inside the horizon of each material point

for i = 1:totnode
    for j = 1:totnode
        if (i~=j)
            xdist = coord(j,1) - coord(i,1);
            ydist = coord(j,2) - coord(i,2);
            
            idist = sqrt( xdist^2 + ydist^2 );
            
            if(idist <= delta)
                
                numfam(i) = numfam(i) + 1;
                
                nodefam(numfam(i)) = j;
            end
        end
        
    end
    
    for j = 1 : numfam(i)
        family(i).node(j) = nodefam(j);
        family(i).nodeFlag(j) = 1;
    end
    
    
end
end





function SCorrection(sedloadm1, sedloadm2, bcd, bcs, radij)
global coord delta vol totnode
global numfam family
global fncstd fncstm

tdisp = zeros(totnode,2);
dilm1 = 0.001;
dilm2 = 0.001;


tdisp(:,1) = dilm1 * coord(:,1);
tdisp(:,2) = 0.0;


for i = 1:totnode
    
    dload1 = 0.0;
    stedens1 = 0.0;
    
    for jmem = 1:numfam(i,1)
        j = family(i).node(jmem);
        ix1 = (coord(j,1) - coord(i,1));
        ix2 = (coord(j,2) - coord(i,2));
        nx1 = (coord(j,1)+tdisp(j,1)-coord(i,1)-tdisp(i,1));
        nx2 = (coord(j,2)+tdisp(j,2)-coord(i,2)-tdisp(i,2));
        
        idist   = sqrt(ix1^2 + ix2^2);
        nlength = sqrt(nx1^2+nx2^2);
        str = (nlength - idist) / idist;
        pplen = nlength * idist;
        
        ppx = nx1 * ix1;
        ppy = nx2 * ix2;
        pptot = 1; %(ppx + ppy) / pplen;
        
        if (idist <= delta-radij)
            fac = 1.0;
        elseif (idist <= delta+radij)
            fac = (delta+radij-idist)/(2.0*radij);
        else
            fac = 0.0;
        end
        
        dload1 = dload1 + bcd * delta * str * pptot * vol * fac;
        stedens1 = stedens1 + bcs * delta * (str^2) * (idist) * vol * fac;
        
    end
    fncstd(i,1) = dilm1 / dload1;
    fncstm(i,1) = sedloadm1 / stedens1;
end



tdisp(:,1) = 0.0;
tdisp(:,2) = dilm2 * coord(:,2);


for i = 1:totnode
    
    dload2 = 0.0;
    stedens2 = 0.0;
    
    for jmem = 1:numfam(i,1)
        j = family(i).node(jmem);
        ix1 = (coord(j,1) - coord(i,1));
        ix2 = (coord(j,2) - coord(i,2));
        nx1 = (coord(j,1)+tdisp(j,1)-coord(i,1)-tdisp(i,1));
        nx2 = (coord(j,2)+tdisp(j,2)-coord(i,2)-tdisp(i,2));
        
        idist   = sqrt(ix1^2 + ix2^2);
        nlength = sqrt(nx1^2+nx2^2);
        str = (nlength - idist) / idist;
        pplen = nlength * idist;
        
        ppx = nx1 * ix1;
        ppy = nx2 * ix2;
        pptot = 1; %(ppx + ppy) / pplen;
        
        if (idist <= delta-radij)
            fac = 1.0;
        elseif (idist <= delta+radij)
            fac = (delta+radij-idist)/(2.0*radij);
        else
            fac = 0.0;
        end
        
        dload2 = dload2 + bcd * delta * str * pptot * vol * fac;
        stedens2 = stedens2 + bcs*delta * (str^2) * (idist) * vol * fac ;
        
    end
    fncstd(i,2) = dilm2 / dload2;
    fncstm(i,2) = sedloadm2 / stedens2;
    
end

end

function ApplySlit()
global coord totnode
global numfam family

% % %crack information
x1=40;
x2=60;
y1=50;
y2=50;
for k = 1:totnode
    for jmem =  1 : numfam(k,1)
        
        j = family(k).node(jmem);
        
        a0x= coord(j,1)-coord(k,1);
        a0y= coord(j,2)-coord(k,2);
        b0x= x1 - coord(k,1);
        b0y= y1 - coord(k,2);
        c0x= x2 - coord(k,1);
        c0y= y2 - coord(k,2);
        ax= x2 - x1;
        ay= y2 - y1;
        bx=coord(k,1) - x1;
        by=coord(k,2) - y1;
        cx=coord(j,1) - x1;
        cy=coord(j,2) - y1;
        delta1=(a0x*b0y-a0y*b0x)*(a0x*c0y-a0y*c0x);
        delta2=(ax*by-ay*bx)*(ax*cy-ay*cx);
        
        if(delta1<=0 && delta2<=0)
            family(k).nodeFlag(jmem)=0;
        end
        
    end
end

end
