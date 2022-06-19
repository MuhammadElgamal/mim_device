classdef mim_device<handle               

    properties
        
        %% Device Structures
        % can be diode, single-gate transistor, double-gate transistor
        geometry;
        % w1 specifies the position of gate electrode, W2 specifies the
        % widths of different dielectrics of the channel region
        widths; % channelwidths, 3-element gate oxide width
        % h1 the hights of several oxide layers below the gate
        % h2 the hight of the channel region
        heights; % 1-element channel height and layer thickness of oxide
        % per1 must be same as h1, per 2 must be same as w2
        permitivity;
        channel_affinity;
        channel_lengths;
        
        %% Device Potential models
        pdemodel;
        potential;
        point_count;
        mesh_longest_edge;
        
        % IV solution charcateristics
        VdsStart ;
        VdsEnd ;
        VdsNumPoints;
        VdsSweepValues;
        VoltageLimit;
        VoltagePoints;
        i_tlm;
        i_tmm;
        i_simmon;
        i_wkb;
        sweepVoltage;
        
        %% Constants
        q=1.6021892e-19 ;
        h=6.62617e-34 ;
        kB  = 1.380662e-23;                 % J/K
        hbar= 1.0545e-34;                   % J.sec
        m0  = 9.10953e-31;                  % Kg
        eps0= 8.854e-12;                    % Permitivity (F/m)
        Temp=300;                         % Temperature (K.elv)
        
        %% Device Current Calculation models
        %Device Left/Right Metals Specs
        WFL;
        WFR;
        Phi_B;
 
        % Tunneling Window Energy Limits
        Ef=0;
        EL=10;
        ER;
        
        EF_D;
        EF_S;
        E_ref_S;
        E_ref_D;
        
        Emin;
        Emax;
        dE=26e-5;
        NE;
        E;
        
        Ns;
        Np;
        % Tunneling Calculation Properties
        
        alfa_S;
        alfa_D;
        
        % the value of this property depends on temprature
        kBT;       


        % Extracted Features
        volt;
        curr;
        resistivity;
        responsivity;
        nonlinearity;
        assymetry;           % assymetry
        
    end
    
    methods
        % Constructor function
        function device=mim_device(geometry,widths, heights, permitivity, affinity, MeshPoints)
            device.point_count = MeshPoints;
            if ~(length(widths)==length(heights) && length(heights)==length(permitivity) && length (widths)==2)
                error('Widths must be 2 cell arrays for channel and oxide regions');
            elseif length(heights{1})~=1
                error('Channel region has a single height');
            elseif length(widths{2})~=3
                error('Oxide region is divided into three only for placing the gate electrode');
            elseif ( abs(sum(widths{2})-sum(widths{1}))  > 1e-6 )
                error('Device cannot have more than a single width');
            elseif any(permitivity{1}<0) || any(permitivity{2}<0)
                error('Dielectric Constants cannot be negative at all');
            elseif length(widths{1})~=length(permitivity{1}) || length(heights{2})~=length(permitivity{2})
                error('There should be a matching of number of permitivities and numver of layers');
            elseif length(affinity)~=length(permitivity{1})
                error('The channel affinity vector must be same length as permitivity of channel')
            else
                factor=50;
                device.geometry=geometry;
                if length(widths{1})<20
                    fac=ceil(20/length(widths{1}));
                    widths{1}=repelem(widths{1}./fac,1,fac);
                    permitivity{1}=repelem(permitivity{1}, 1,fac);
                    affinity=repelem(affinity,1,fac);
                end
                device.widths=widths;
                device.permitivity=permitivity;
                device.channel_affinity=affinity;
                device.channel_lengths = widths{1};
                
                if strcmp(geometry, "diode")
                    %heights{2}=ones(size(permitivity{2}));
                    heights{2}=0.5*sum(widths{1})/factor;
                end
                device.heights=heights;
                % original meshing factor 50
                device.mesh_longest_edge=sum(widths{1})/factor;
            end
           
        end
        
        function VdsWindowLimit(device,VdStartVal,VdEndVal,NumPoints)
            
            %% setConstatnts
            device.kBT=(1.380662e-23/device.q)*device.Temp;
          
            %% Set Start and End Potentials
            device.VdsStart=VdStartVal;
            device.VdsEnd=VdEndVal;
            device.VdsNumPoints = NumPoints;
            device.VdsSweepValues = linspace(VdStartVal,VdEndVal,NumPoints);            
        end
        
        function exit=find_potential(device, potentials)
            if length(potentials)~=3
                error('Enter [vg vd vs] replace vg with 0 for diode');
            end
            %trial_time=1;
            trial_time=0.3;
            exit = 0;
            numOfTrials = 5;
            for kkk=1:numOfTrials
                try
                    switch (device.geometry)
                        case "diode"
                            % S is senstivity number of points e.g. 1000;
                            % W is the width vector of several materials
                            % h is the hieght of the MnIm
                            % e is the vector of permitivity
                            S=device.point_count;
                            w=device.widths{1};
                            h=device.heights{1};
                            e=device.permitivity{1};
                            meshparam=device.mesh_longest_edge;
                            vd=potentials(2);
                            vs=potentials(3);
                            
                            model = createpde();
                            n=length(w);
                            minw=0;
                            for i=1:n-1
                                minw=[minw sum(w(1:i))];
                            end
                            
                            maxw=minw+w;
                            gd=[];
                            ns=[];
                            sf=[];
                            for i=1:n
                                Ri=[3;4;minw(i);maxw(i);maxw(i);minw(i);0;0;h;h];
                                gd=[gd,Ri];
                                Rstr=strcat('R',num2str(i,'%05d'));
                                ns=[ns, Rstr'];
                                if i<n
                                    sf=[sf Rstr '+'];
                                else sf=[sf Rstr];
                                end
                                
                            end
                            
                            [dl,~] = decsg(gd,sf,ns);
                            geometryFromEdges(model,dl);
                            % pdegplot(model);
                            if n~=1
                                facevector=[n-1, (1:n-2), n];
                            end
                            
                            if n==1
                                applyBoundaryCondition(model,'neumann','Edge',[1,3],'g',0,'q',0);
                                applyBoundaryCondition(model,'dirichlet','Edge',4,'h',1,'r',vd);
                                applyBoundaryCondition(model,'dirichlet','Edge',2,'h',1,'r',vs);
                                specifyCoefficients(model,'m',0,'d',0,'c',e,'a',0,'f',0,'Face',1);
                            else
                                applyBoundaryCondition(model,'neumann','Edge',(3:(3*n+1)),'g',0,'q',0);
                                applyBoundaryCondition(model,'dirichlet','Edge',1,'h',1,'r',vd);
                                applyBoundaryCondition(model,'dirichlet','Edge',2,'h',1,'r',vs);
                                for i=1:n
                                    specifyCoefficients(model,'m',0,'d',0,'c',e(i),'a',0,'f',0,'Face',facevector(i)); %rho=0 because of laplace equation not poisson
                                end
                            end
                            
                            generateMesh(model,'Hmax',meshparam);
                            results = solvepde(model);
                            xq=linspace(0,sum(w),S);
                            yq=linspace(0,h,S);
                            [X,Y] =meshgrid(xq,yq);
                            querypoints=[X(:),Y(:)]';
                            V=interpolateSolution(results,querypoints);
                            V = reshape(V,size(X));
                        case "single-gate transistor"
                            device.Ef=device.EL;
                            % Onoff controls the appearance of structure figure
                            meshparam=device.mesh_longest_edge;
                            vg=potentials(1);
                            vd=potentials(2);
                            vs=potentials(3);
                            w1=device.widths{1};
                            w2=device.widths{2};
                            h1=device.heights{1};
                            h2=device.heights{2};
                            e1=device.permitivity{1};
                            e2=device.permitivity{2};
                            S=device.point_count;
                            
                            model = createpde();
                            n=length(w1);  % How many horizontal sections
                            m=length(h2);  % How many Vertical Sections
                            minw1=0;
                            for i=1:n-1
                                minw1=[minw1 sum(w1(1:i))];
                            end
                            maxw1=minw1+w1; % specifying the start and end of each  horizontal region
                            
                            minw2=0;
                            for i=1:2
                                minw2=[minw2 sum(w2(1:i))];
                            end
                            maxw2=minw2+w2;  % specifying the start and end of each  horizontal region
                            
                            minh2=0;
                            for i=1:m-1
                                minh2=[minh2 sum(h2(1:i))];
                            end
                            maxh2=minh2+h2;  % specifying the start and end of each  vertical region
                            
                            gd=[];
                            ns=[];
                            sf=[];
                            for i=1:n
                                Ri=[3;4;minw1(i);maxw1(i);maxw1(i);minw1(i);0;0;-h1;-h1];
                                gd=[gd,Ri];
                                ll=randi([97,122],1,1); % In order to add a small charcter that designates every region
                                mm=randi([97,122],1,1);
                                Rstr=['R',char(ll) char(mm)];
                                ns=[ns, Rstr'];
                                if i<n
                                    sf=[sf Rstr '+'];
                                else sf=[sf Rstr];
                                end
                                
                            end
                            for i=1:m
                                Rj1=[3;4;minw2(1);maxw2(1);maxw2(1);minw2(1);minh2(i);minh2(i);maxh2(i);maxh2(i)];
                                Rj2=[3;4;minw2(2);maxw2(2);maxw2(2);minw2(2);minh2(i);minh2(i);maxh2(i);maxh2(i)];
                                Rj3=[3;4;minw2(3);maxw2(3);maxw2(3);minw2(3);minh2(i);minh2(i);maxh2(i);maxh2(i)];
                                
                                gd=[gd,Rj1,Rj2,Rj3];
                                ll=randi([97,122],1,1); % In order to add a capital charcter that designates every region
                                mm=randi([97,122],1,1);
                                Rstr1=['X',char(ll) char(mm)];
                                Rstr2=['Y',char(ll) char(mm)];
                                Rstr3=['Z',char(ll) char(mm)];
                                ns=[ns, Rstr1',Rstr2',Rstr3'];
                                sf=[sf '+' Rstr1 '+' Rstr2 '+' Rstr3];
                                
                            end
                            
                            
                            [dl,~] = decsg(gd,sf,ns);
                            geometryFromEdges(model,dl);
                            figure(9)
                            pdegplot(model,'EdgeLabels','on', 'FaceLabels', 'on');
                            j=SegName(model,dl,w2(1),sum(h2),w2(1)+w2(2),sum(h2));
                            applyBoundaryCondition(model,'dirichlet','Edge',j,'h',1,'r',vg);
                            
                            j=SegName(model,dl,0,0,0,-h1);
                            applyBoundaryCondition(model,'dirichlet','Edge',j,'h',1,'r',vd);
                            j=SegName(model,dl,sum(w1),0,sum(w1),-h1);
                            applyBoundaryCondition(model,'dirichlet','Edge',j,'h',1,'r',vs);
                            facevector=[n-2,n-1,(1:n-3),n];
                            % Assigning the permitivity for the channel region
                            for i=1:n
                                specifyCoefficients(model,'m',0,'d',0,'c',e1(i),'a',0,'f',0,'Face',facevector(i)); %rho=0 because of laplace equation not poisson
                            end
                            a=(n+4:3*m+n);
                            facematrix=[n+(1:3);reshape(a,length(a)/3,3)];
                            for i=1:m
                                specifyCoefficients(model,'m',0,'d',0,'c',e2(m-i+1),'a',0,'f',0,'Face',facematrix(i,:)); %rho=0 because of laplace equation not poisson
                            end
                            % the Hgrad option specifies the growth rate outside
                            % the tiny-meshed regions
                            generateMesh(model,'Hmax',meshparam, 'Hmin', 0.5*meshparam, 'Hgrad', 1.7);
                            results = solvepde(model);
                            xq=linspace(0,sum(w1),S);
                            yq=linspace(-h1,sum(h2),S); %%plots the voltage accross all device
                            %yq=linspace(-h1,0,S); %only voltage accross the
                            [X,Y] =meshgrid(xq,yq);
                            querypoints=[X(:),Y(:)]';
                            V=interpolateSolution(results,querypoints);
                            V = reshape(V,size(X));
                        case "double-gate transistor"
                    end
                    device.pdemodel=model;
                    device.potential.v=V;
                    device.potential.x=X;
                    device.potential.y=Y;
                    break;
                catch
                    pause(trial_time);
                    fprintf('Fail Trial %d\n', kkk);
                    if(kkk == numOfTrials)
                        exit=-1;
                    end
                end
            end
        end
        
        function plot_device(device, struct)
            % if struc is 1 it draws the structure if 2 it draws the
            % dielectric constant of the device
            switch (nargin)
                case 1
                    %device.potential;
                    mesh(device.potential.x, device.potential.y, device.potential.v);
                    xlabel('Width of Device');
                    ylabel('Thickness of Device');
                    view(0,90);
                    title('Channel Potential (V)')
                    colorbar;
                    
                case 2
                    if struct==1
                        pdegplot(device.pdemodel,'EdgeLabels','on', 'FaceLabels', 'on');
                        xlabel('Width of Device');
                        ylabel('Thickness of Device');
                        view(0,90);
                    else
                        X=device.potential.x;
                        Y=device.potential.y;
                        p=zeros(size(X));
                        switch (device.geometry)
                            case 'diode'
                                w1=device.widths{1};
                                start=0;
                                e1=device.permitivity{1};
                                for i=1:length(e1)
                                    p(start<=X & X<start+w1(i))=e1(i);
                                    start=start+w1(i);
                                end
                            case "single-gate transistor"
                                w1=device.widths{1};
                                start=0;
                                e1=device.permitivity{1};
                                for i=1:length(e1)
                                    p(start<=X & X<start+w1(i) & Y<0)=e1(i);
                                    start=start+w1(i);
                                end
                                
                                h2=device.heights{2};
                                start=0;
                                e2=device.permitivity{2};
                                for i=1:length(e2)
                                    p(start<=Y & Y<start+h2(i))=e2(i);
                                    start=start+h2(i);
                                end
                                h2=device.heights{2};
                                we=device.widths{2};
                                p (we(1)<X & X<we(1)+we(2) & sum(h2)-0.1*h2(end)<Y)=max(p(:))*0.01;
                            case "double-gate transistor"
                        end
                        
                        
                        
                        mesh(X,Y,p)
                        xlabel('Width of Device');
                        ylabel('Thickness of Device');
                        view(0,90);
                        title('Permitivity')
                        colorbar;
                        %axis([0 sum(device.widths{1}) -device.heights{1}, sum(device.heights{2})])
                        
                    end
            end
            switch (device.geometry)
                case "diode"
                    axis([0 sum(device.widths{1}) 0 device.heights{1}]);
                case "single-gate transistor"
                    axis([0 sum(device.widths{1}) -device.heights{1}, sum(device.heights{2})]);
            end
        end
        
        function [dis,v] = takeslice(device, xory, loc)
            % the function takes a slice at a specific x or specific y
            X=device.potential.x;
            Y=device.potential.y;
            xq=X(1,:);
            yq=Y(:,1)';
            pot=device.potential.v;
            switch (xory)
                case 'x'
                    dis=yq;
                    tol=0.5*(abs(xq(2)-xq(1)));
                    v=pot(abs(X-loc)<= tol)';
                    
                case 'y'
                    dis=xq;
                    tol=0.5*(abs(yq(2)-yq(1)));
                    v=pot(abs(Y-loc)<= tol)';
            end
        end
        

        % Current Calculation
        function setMetalsSpecs(obj,LeftFermi,LeftWorkFunction,RightWorkFunction,Vdrain)
            % Define Constatnts MIM Device Metals
            obj.WFL = LeftWorkFunction ;
            obj.WFR = RightWorkFunction ;
            obj.Phi_B= obj.WFL;
            obj.EL = LeftFermi;
            obj.ER = LeftFermi;
            obj.E_ref_S = 0; % The reference potential energy at the source side (Ec at source contact) is taken as reference
            obj.EF_S = obj.E_ref_S + obj.EL; %E_ref_S + 20*kBT; % The Fermi level inside the source contact is assumed to be 20kT above Ec (i.e. 20 KT above E_ref_S)

        end
        
        function Energy_states = TunnelingWindow(obj,UB,Vdrain,Vdrain_window,calculationType)
            
            obj.dE = 0.001;

            if( calculationType  == "TLM_TMM")
                
                
                %% Defiene Drain and Source Fermi Levels
                obj.E_ref_D =  Vdrain; % E_ref_D is the potentail energy in drain contact, which is Ec inside drain
                obj.EF_D = obj.EF_S - Vdrain; % The Fermi level at the drain in eV
            
                % obj.EF_S =obj.E_ref_S+ obj.EL; %E_ref_S + 20*kBT; % The Fermi level inside the source contact is assumed to be 20kT above Ec (i.e. 20 KT above E_ref_S)

                % Ec inside the drain (E_ref_D) is assumed to be  20kT below EF_D
                %Emin = 9; Emax = 13; % maximum energy at which carriers may be found

                %VV=[ obj.VdsStart obj.VdsEnd ];
                VV = Vdrain_window;
                % Calcualte max energy level and minimum levels
                obj.Emin = min(obj.EF_S, obj.EF_S - max(VV))-100* obj.kBT;
                obj.Emax = max(obj.EF_S, obj.EF_S - min(VV))+100*obj.kBT; % you have at least add (1.4 - 2 )over max potential val
                obj.E = (obj.Emin:obj.dE:obj.Emax)';

                
            elseif ( calculationType == "transistor" )
            

                % For Tunneling 
                obj.Emin =obj.Ef-Vdrain-(10*obj.kBT); 
                obj.Emax = UB(2) %obj.Ef+obj.Phi_B;  %the energy of the tunnel electrons only

                % For Therminoic
                %Emin=Ef+phib; 
                %Emax=Ef+phib+15*kbT;  %thermoionic current component only

                %the energy of both the tunnel and thermoionic current componenents  
                % Emin=Ef-Vd-10*kbT; 
                % Emax=Ef++phib+15*kbT;  

                obj.E = (obj.Emin:obj.dE:obj.Emax);
            end
            obj.NE = length(obj.E);
            Energy_states = obj.E;
%             obj.E_ref_S = 0; % The reference potential energy at the source side (Ec at source contact) is taken as reference
%             obj.EF_S =obj.E_ref_S+ obj.EL; %E_ref_S + 20*kBT; % The Fermi level inside the source contact is assumed to be 20kT above Ec (i.e. 20 KT above E_ref_S)
%             
%             obj.EF_D = obj.EF_S - Vdrain; % The Fermi level at the drain in eV
%             obj.E_ref_D =  Vdrain; % E_ref_D is the potentail energy in drain contact, which is Ec inside drain
%             
%             % Ec inside the drain (E_ref_D) is assumed to be  20kT below EF_D
%             %Emin = 9; Emax = 13; % maximum energy at which carriers may be found
%             
%             VV=[ obj.VdsStart obj.VdsEnd ];
%             % Calcualte max energy level and minimum levels
%             obj.Emin = min(obj.EF_S, obj.EF_S - max(VV))-20* obj.kBT;
%             obj.Emax = max(obj.EF_S, obj.EF_S - min(VV))+20*obj.kBT; % you have at least add (1.4 - 2 )over max potential val
%             obj.E = (obj.Emin:obj.dE:obj.Emax)';
%             obj.NE = length(obj.E);
%             Energy_states = obj.E ;
            
        end
        
        
        % Characterisation
        function exit=characterise(device, variable_to_fix, fixed_value,plotfigures)
            
%             if device.geometry == "single-gate transistor"
%                 device.sweepVoltage=linspace(0,device.VoltageLimit,device.VoltagePoints);
%             else
                device.sweepVoltage=linspace(-device.VoltageLimit,device.VoltageLimit,device.VoltagePoints);
%             end
            
            % Initialization 
            V = device.sweepVoltage;
            I_TLM_final = [];
            I_TMM_final = [];
            I_Simmon_final = [];

%             slices = linspace(0, -device.heights{1}, slice_count);
            device.i_tlm=zeros(1, device.VoltagePoints);
            
            for index_Vd = 1:length(V)                
                Vg=0;
                switch (variable_to_fix)
                    case 'vgs'
                        Vg=fixed_value;
                        Vd =V(index_Vd);
                        disp("Vds = " + Vd)

                    case 'vds'
                        Vd=fixed_value;
                        Vg =V(index_Vd);
                        disp("Vgs = " + Vg)
                end                
                Vs=0;                
                exit=find_potential(device, [Vg, Vs, Vd]);
                if(exit == -1)
                    return;
                end
                fprintf('Solved for Channel Length of %f\n',  sum(device.widths{2}) );
                [Xpde,Vpde] = takeslice(device, 'y', 0); 
                
%                 if device.geometry == 'single-gate transistor'
%                     Vpde=(device.Phi_B+device.Ef)-Vpde;
%                 end
                
                %figure;plot(Xpde,Vpde)
                figure(5);
                plot(Xpde,-Vpde)
                %[xx,d1,pot_barrier]=Analytical_Pot(Vg,Vd,E)
                VdsWindowLimit(device,0,device.VoltageLimit,device.VoltagePoints);
                [X_mid,E_mid,E_mid_Equil] = CalculateEBD(device,Xpde,Vpde,Vd,1,device.widths{1},...
                                                         device.permitivity{1},device.channel_affinity,device.geometry);
                [meff_mid,w,xx,UB] = device.ImageForceEffect(X_mid,E_mid,E_mid_Equil,Vd, plotfigures,device.geometry);
                [T_TLM,T_TMM,T_Simmon,T_WKB,I_TMM,I_TLM,I_Simmon,I_WKB] = device.TunnelingCurrentCalc(meff_mid,w,UB,Vd,[V(1) V(end)],xx(1:end-1),plotfigures);
                
                device.i_tlm(index_Vd) = device.i_tlm(index_Vd) + I_TLM;
                device.i_tmm(index_Vd) = I_TMM;
                device.i_simmon(index_Vd) = I_Simmon;
                device.i_wkb(index_Vd) = I_WKB;
                
            end
            
            end_val=0
        end
        
        function [X_mid,E_mid,E_mid_Equil] = CalculateEBD(obj,Xpde,Vpde,Vd,pdeFlag,InsulatorsLen,Dielectrics,Affinity,geometry)
            
            
            E_ref_S  = obj.E_ref_S;
            EF_S = obj.EF_S;

            TotLength = sum(obj.widths{1});
            lengthL =InsulatorsLen;
          
            InsDielectric = obj.permitivity{1};
            NumInsulators = length(InsDielectric);
            
            numMeshingPoints = obj.point_count;
            
            %--Get the Percentage for Num of points at each insualtor according to itsEF_S length------------%
            InsulatorRatioVect=floor(numMeshingPoints*(lengthL./TotLength));
            %Find difference in num of points between formed vector and received vector
            pointsDiff=numMeshingPoints-sum(InsulatorRatioVect);
            %If there is a difference, increase last insulator vector by this difference
            InsulatorRatioVect(end)=InsulatorRatioVect(end)+pointsDiff;
            
            x_mid_PDE = Xpde;
            % Initializations
            EL_max=[];
            EL_min =[];
            PrevlengthL=[0 lengthL];
            XLength=[];
            E_mid_Equil=[];
            X_mid=[];
            E_mid=[];
            Vbias = Vd;


            % if Analytical Solution
            if pdeFlag == 0
                
                Vbi = obj.WFL - obj.WFR;
                for i=1:NumInsulators
                    %Calc indices for each insulator start(min) and end(max)
                    EL_max(i)=sum(InsulatorRatioVect(1:i));
                end
                EL_min=EL_max-InsulatorRatioVect+1; %[1 31 81]
                
                for i=1:NumInsulators
                    
                    XLength=x_mid_PDE(EL_min(i):EL_max(i));
                    
                    if(NumInsulators == 1)
                        %Ux = (U(i) -(dPHI+VD)*xx/L-(lamda*L)./(xx.*(L-xx)))'; %barrier including image force
                        %dV = (VD+Vbi).*(XLength./L1);
                        dV(i) = (-Vbias+Vbi)*(lengthL(i)/InsDielectric(i))/sum(lengthL./InsDielectric);
                        
                    else
                        dV(i) = (-Vbias-Vbi)*(lengthL(i)/InsDielectric(i))/sum(lengthL./InsDielectric);
                    end
                    
                    %for First isulator
                    if i == 1
                        
                        U(i) =EF_S+(WFL-EA(1));
                        
                        %If reached to last Insulator
                    elseif i == NumInsulators
                        
                        U(i) = (EF_S+Vbias)+(WFR-EA(end)+dV(i));
                        %Else for any Inside insulator
                    else
                        %Here Still Needs Some Check!
                        U(i) = (EF_S+Vbias)+(EA(i)-EA(i-1))+dV(i);
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    E_mid_Eq=linspace(U(i),U(i),InsulatorRatioVect(i));
                    E_mid_Equil=[E_mid_Equil E_mid_Eq ];
                    %Sum PDE Pot on Equilib Pot
                    disp('-------------------------');
                    
                    E_mid_new=E_mid_Eq-(dV(i)*(XLength-PrevlengthL(i))./lengthL(i));
                    E_mid=[E_mid E_mid_new ];
                    
                    X_mid_new=x_mid_PDE(EL_min(i):EL_max(i));
                    X_mid=[X_mid X_mid_new ];
                    
                end
                
            % if Numerical Solution
            else
                
                % Indices = find(diff(obj.widths{1})~=0);
                % Indices = [Indices length(obj.widths{1})]
                
                % Calculate the correct ratio for each insulator
                Ratio = floor( obj.point_count*InsulatorsLen./TotLength );
                
                % Verify proper Ratio values
                if( sum(Ratio) > obj.point_count )
                    
                    error('Ratio is wrong')

                elseif ( sum(Ratio) < obj.point_count )

                    diff =  obj.point_count - sum(Ratio) ; 
                    Ratio = [ Ratio(1:end-1) Ratio(end)+diff ];

                end

                E_mid_Eq = [];
                E_mid = [];
                X_mid =[];
               
                % Equilibrium Energy Band Diagram
                for i = 1:length(Ratio)

                    Equi_Pot = obj.EF_S+(obj.WFL-Affinity(i));
                    E_mid =linspace(Equi_Pot,Equi_Pot,Ratio(i));
                    E_mid_Equil = [ E_mid_Equil E_mid];
                end
                
                if geometry == "diode"
                    % Calculate Barrier Heights at each inuslator, and Obtain

                    E_mid=E_mid_Equil-Vpde;
                    E_mid(1)=0;
                    E_mid(end)=-Vbias;
                else
                    %E_mid=E_mid_Eq+Vpde;
                    %Equi_Pot = obj.EF_S+(obj.WFL-Affinity(i));
                    E_mid=obj.EF_S+obj.WFL-Vpde;
                    E_mid(1)=0;
                    E_mid(end)=-Vbias;
                end
                X_mid = x_mid_PDE;

            end
        end
        
        function [meff_mid,w,xx,UB] = ImageForceEffect(obj,X_mid,E_mid,E_mid_Equil,VD, plotfigures,geometry)
            NofDivs=length(X_mid);
            InsDielectric=obj.permitivity{1};
            obj.Np = obj.point_count;
            obj.Ns = obj.Np - 1; % no of segments (mesh elements)
            
           if geometry == "diode"
                Vnew=ImageForceFunction(InsDielectric,X_mid,E_mid,NofDivs);
                Vnew(1) = 0;
                Vnew(end) =-VD;    
           else
               % for transistor don`t calculate the image potential
               Vnew = E_mid;
           end

            
            %plot(Vnew)
            %---------------------- Effective mass array ------------------------------%
            % Determine the effective mass and potential distribution
            % meff_array = zeros(Np,1); % Initialize the effective mass array
            % meff is constant and equals to that of Al(0.5)_Ga(0.5)_As
            
            % this part will be further generalized to allow inserting different
            % effective masses for all insulators
            
            meff = 1* obj.m0;
            meff_array = meff*ones(obj.Np,1); % meff at every node (used by FDM)
            meff_mid = (meff_array(1:end-1)+ meff_array(2:end))/2; % meff at every segment (used by TMM)
            
            %--------------------------------------------------------------%
            % potential parameters
            %U0=EbarrierVect;
            %-----------------------------%
            L = sum(obj.widths{1}); % the width of the barrier in nm
            w = (L/obj.Ns); % Tunneling Window mesh spacing in nm
            %xx =0:w:L;% position vector within the barrier         % meff_B is the effective mass in either left or right boundary region in Kg
            xx=linspace(0,L,length(Vnew));        
            
            %U0= E_mid_Equil;%Vnew(1); %EF_S+(WFL-EA(1)); %
            Ub =Vnew;
            %------------------------ Plot Potential Barrier --------------------%
            if plotfigures
                
                fig6=figure(6);
                
                plot(X_mid,E_mid,X_mid,E_mid_Equil,X_mid,Vnew,'markersize',5,'linewidth',2,'XDataSource','real(x)','YDataSource','imag(y)')
                xlabel('{\itx}[nm]','fontsize',14)
                ylabel('{\itU} [eV]','fontsize',14)
                lg=legend('Potential After Bias','Equilibrium Potential','Image Potential','Location','Best')
                %lg.Location= 'southWest'
                linkdata(fig6);
            end
            % fig2=figure(2);
            % plot(xx,Ub,'b-','markersize',5,'linewidth',2,'XDataSource','real(x)','YDataSource','imag(y)')
            % linkdata(fig2)
            % %---- Set The font type and size for the figure ---%
            
            %  %---- Set title for x, y axes ----%
            %  xlabel('{\itx}[nm]','fontsize',14)
            %  ylabel('{\itU} [eV]','fontsize',14)
            %Ub=Ub';
            %UB(1) = 0; UB(end) =VD;
            UB = (Ub(1:end-1)+Ub(2:end))/2; % The potential in each segment is taken as
            UB=UB';

        end
        
        function [T_TLM,T_TMM,T_Simmon,T_WKB,I_TMM,I_TLM,I_Simmon,I_WKB] = TunnelingCurrentCalc(obj,meff_mid,w,UB,Vd,Vd_window,x,plotfigures)
            % Variables Initialization
            T_Simmon=0;
            I_Simmon=0;
            T_WKB=0;
            I_WKB=0;

            T_Simmon=0;
            I_Simmon = 0;
            T_WKB=0;
            I_WKB = 0;

            TunnelingWindow(obj,UB,Vd,Vd_window,'TLM_TMM');
            E_TLM_TMM=obj.E;

            obj.alfa_S=1e-9*sqrt(2*obj.m0*obj.q*(obj.E_ref_S-E_TLM_TMM))/obj.hbar; % wave number in source % nm-1
            obj.alfa_S(obj.alfa_S==0)= obj.eps0*1i; % to avoid division by zero
            obj.alfa_D=1e-9*sqrt(2*obj.m0*obj.q*(obj.E_ref_D-E_TLM_TMM))/obj.hbar; % wave number in drain % nm-1

            %% Default Tunneling Current Calculations for both Diode and Transistor
            [T_TMM,I_TMM] = CurrentTMM(obj,UB,w,meff_mid,obj.alfa_S,obj.alfa_D);
            [T_TLM,I_TLM] = CurrentTLM(obj,meff_mid,w,UB);

            %           [T_WKB,I_WKB]
            % Only for Transistor
            if( obj.geometry == "single-gate transistor" )

                tch=obj.heights{1};
                TunnelingWindow(obj,UB',Vd,0.01,'transistor');
                E_WKB=obj.E;

                %[T_Simmon,I_Simmon] = CurrentSimmmon(obj,UB',x,-Vd_vect);
                [T_WKB,I_WKB] = CurrentWKB(obj,UB',Vd,x);
                I_WKB=I_WKB*tch*1e-3;

            end

            if plotfigures
                %                 fig7 = figure(7)
                %                 % Draw using (r^) which is red arrow
                %                 plot(E_TLM_TMM,T_TMM,'r--','linewidth',1,'XDataSource','real(x)','YDataSource','imag(y)');
                %                 hold on
                %                 %----- Plot the Tunneling TLM two times
                %                 % Draw using (b-.) which is blue dot and dash line
                %                 semilogy(E_TLM_TMM,T_TLM,'b-.','linewidth',1,'XDataSource','real(x)','YDataSource','imag(y)');
                %                 hold on

                fig7 = figure(7)
                % Draw using (r^) which is red arrow
                semilogy(E_TLM_TMM,T_TMM,'r--','linewidth',1,'XDataSource','real(x)','YDataSource','imag(y)');
                hold on
                %----- Plot the Tunneling TLM two times
                % Draw using (b-.) which is blue dot and dash line
                plot(E_TLM_TMM,T_TLM,'b-.','linewidth',1,'XDataSource','real(x)','YDataSource','imag(y)');
                hold on
                savefig('TunnelingProb(0.4V)_ResonantTunneling.fig')

                if( obj.geometry == "single-gate transistor" )
                    plot(obj.E,T_WKB,'v-.','linewidth',1,'XDataSource','real(x)','YDataSource','imag(y)');
                    hold on
                end

                % Draw using (m-.) which is blue dot and dash line
                %plot(obj.E,T_Simmon,'m-.','linewidth',1,'XDataSource','real(x)','YDataSource','imag(y)');
                xlabel('E{\it[ev]}','fontsize',14)
                ylabel('T{\it[E]}','fontsize',14)
                lg=legend('Tunneling TMM','Tunneling TLM','Tunneling WKB','Location','Best')
                linkdata(fig7);
                hold off
            end

        end

        %         function [T_TLM,T_TMM,I_TMM,I_TLM,T_Simmon,I_Simmon] =Current_from_slice(device,electrode_potentials,slice_loc)
        %                 find_potential(device, electrode_potentials);
        %                 if nargin==2
        %                     [Xpde,Vpde] = takeslice(device, 'y', 0);
        %                 else
        %                     [Xpde,Vpde] = takeslice(device, 'y', slice_loc);
        %                 end
        %                 VD=Vpde(1);
        %                 VS=Vpde(end);
        %                 VDS= VD - VS;
        %                 TunnelingWindow(device,VDS);
        %                 [X_mid,E_mid,E_mid_Equil] = CalculateEBD(device,Xpde,Vpde,VDS,1);
        %                 [E,meff_mid,w,~,UB] = device.ImageForceEffect(X_mid,E_mid,E_mid_Equil,VDS, false);
        %                 [T_TLM,T_TMM,T_Simmon,I_TMM,I_TLM,I_Simmon] = device.TunnelingCurrentCalc(E,meff_mid,w,UB,VDS);
        %         end
        function extract_features(device, type, smoothing)
            v = device.sweepVoltage;           
            switch (type)
                case 'TLM'
                    i = device.i_tlm;
                case 'TMM'
                    i = device.i_tlm;
                case 'simmon'
                    i = device.i_simmon;
                case 'wkb'
                    i = device.i_wkb;
            end
            i = smooth(i, smoothing(1))';
            vq =linspace(min(v), max(v), 1000);
            i = interp1(v,i,vq,"pchip");
            v = vq;
            device.volt =v;
            device.curr =i;
            i1 = gradient(i,v);
            i1=smooth(i1,smoothing(2))';
            i2 = gradient(i1,v);
            i2=smooth(i2,smoothing(3))';
            device.resistivity = 1./i1;
            device.responsivity =i2./ (2.*i1);
            device.nonlinearity = abs(i1.*v./i);
            vq = [linspace(min(v), 0, 1000) linspace(0,max(v),1000)];
            iq = interp1(v,i,vq);
            device.assymetry.x=vq;
            device.assymetry.y=zeros(size(vq));
            device.assymetry.y(vq>=0)=iq(vq<=0)./iq(vq>=0);
            device.assymetry.y(vq<=0)=iq(vq>=0)./iq(vq<=0);
            device.assymetry.y = abs(device.assymetry.y);
        end

    end
end
function str = SegName(model,dl,x1,y1,x2,y2)

    wgeom(dl,'GetCoor');
    MaxNum=model.Geometry.NumEdges;
    %str=MaxNum;
    for i=1:MaxNum
        [xstart,ystart] =GetCoor(i,0);
        [xend,yend] =GetCoor(i,1);
        if(all([xstart,ystart] ==[x1,y1]) && all([xend,yend] ==[x2,y2]) || all([xstart,ystart] ==[x2,y2] & [xend,yend] ==[x1,y1]))
            str=i;
            break;
        end

    end
    
end
function Vnew = ImageForceFunction(er,xinit,Vxinit,NofDivs)

    w=length(xinit);
    q=1.6e-19; e0=8.85e-14; % in cm units
    Z=length(NofDivs);
    Vimage(1:NofDivs(Z))=0; Vimage_rev=Vimage;
    % for j=1:Z
    %     if j==1
    %         A=2; B=NofDivs(j)-1;
    %     else
    %         A=NofDivs(j-1); B=NofDivs(j)-1; 
    %     end
    %     Vimage(A:B)=q./(16*pi*er(j)*e0*1e-7*(xinit(A:B)-xinit(1))); % on left side
    %     Vimage_rev(A:B)=q./(16*pi*er(j)*e0*1e-7*(max(xinit)-xinit(A:B))); %on right side
    % end
    %Vnew=Vxinit-Vimage-Vimage_rev;


    % new code for modification in image force % march 2009
    coeff=q./(16*pi*e0*1e-7);
    erXfwd(1:NofDivs(Z))=0; erXrev(1:NofDivs(Z))=0;
    for j=1:Z
        if j==1
           % disp('first image j=1');
            A=2 ;
            B=NofDivs(j)-1;
            erXfwd(A:B)=er(j)*(xinit(A:B)-xinit(1));
        else
             %disp('any other image j>1');
             A=NofDivs(j-1);
             B=NofDivs(j)-1 ;
             erXfwd(A:B)=er(j)*(xinit(A:B)-xinit(A-1))+erXfwd(A-1);
        end
    end

    for j=Z:-1:1
        if j==1
             A=NofDivs(j)-1; B=2;
            erXrev(B:A)=er(j)*(xinit(A+1)-xinit(B:A))+erXrev(A+1);
        else
            A=NofDivs(j)-1; B=NofDivs(j-1);
            erXrev(B:A)=er(j)*(xinit(A+1)-xinit(B:A))+erXrev(A+1);
        end
    end
    Vnew=Vxinit;
    Vnew(2:NofDivs(Z)-1)=Vxinit(2:NofDivs(Z)-1)-coeff.*(1./erXfwd(2:NofDivs(Z)-1)+1./erXrev(2:NofDivs(Z)-1));

    % end new code
    %% commented on jan 2011, to correct for Ni-NiO(1.5)-NbOx(0.5)-NbN 
    % n=1;
    % while Vnew(n+1)<Vnew(1)&&n<w-1  % left side for V<EFL
    %     n=n+1;
    % end
    % Vnew(1:n)=Vnew(1); %left side for V<EFL
    % 
    % m=w;
    % while Vnew(m-1)<Vnew(w)&&m>2
    %     m=m-1;
    % end
    % Vnew(m:w)=Vnew(w);
    % 
    % if n==w-1||m==2  %if image force penetrates complete barrier
    %     state='no barrier with image force, use orginal barrier'
    %     Vnew=Vxinit;
    % end
    %% added in liew of the above comment, to correct -ve Vnew values.
    n=1;
    while Vnew(n+1)<Vnew(1)&&abs(Vnew(n+1)-Vxinit(n+1))>0.1  % left side for V<EFL
    %     Vnew(n+1)-Vxinit(n+1);
        n=n+1;
    end
    Vnew(1:n)=Vnew(1); %left side for V<EFL

    m=w;
    while Vnew(m-1)<Vnew(w)&&abs(Vnew(m-1)-Vxinit(m-1))>0.1
        m=m-1;
    end
    Vnew(m:w)=Vnew(w);

    if n==w-1||m==2  %if image force penetrates complete barrier
        state='no barrier with image force, use orginal barrier';
        disp(state);
        Vnew=Vxinit;
    end

end
function [T_TLM,I_TLM] = CurrentTLM(EBD,meff_mid,w,UB)
    
 

    Ns = EBD.Ns;
    E_ref_D = EBD.E_ref_D;
    meff_B = EBD.m0; 
    NE = EBD.NE ;
    E = EBD.E;
    hbar = EBD.hbar;
   
    kBT = EBD.kBT;
    q = EBD.q;
    EF_D = EBD.EF_D;
    EF_S = EBD.EF_S;
    
    
    NN = Ns+2; % no of points in addition to left and right boundary points
    UU = [0;UB;E_ref_D]; % potential energy array used for TLM
    mm = [meff_B;meff_mid;meff_B]; % effective mass array used for TLM
    % Tr = zeros(1,NE);
    gg = zeros(NN,NE);
    Zo = zeros(NN,NE);


    for i = 1:NE
       gg(:,i) = 1i*(1e-9)*(sqrt(2*q*mm.*(E(i) - UU)/hbar^2)); % nm-1
       Zo(:,i) = gg(:,i)*hbar./(1i*mm);
    end
    %    gg(gg==0)= eps; % to avoid division by zero
    %    Zo(Zo==0)= eps; % to avoid division by zero
    ZL = Zo(NN,:);
    for n = NN-1:-1:2
         ZL=Zo(n,:).*((ZL+Zo(n,:).*tanh(gg(n,:)*w))./(Zo(n,:)+ZL.*tanh(gg(n,:)*w)));
    %        ZL=Zo(n,:).*((ZL.*cosh(gg(n,:)*w)+Zo(n,:).*sinh(gg(n,:)*w))./(Zo(n,:).*cosh((gg(n,:))*w)+ZL.*sinh(gg(n,:)*w)));
    end
    rho=(ZL-Zo(1,:))./(ZL+Zo(1,:));
    Tr = 1-(abs(rho)).^2;
    Tr=Tr';
    %current density in A/cm2
    if size(E,1)>1
        %I = 2*q^2/(pi*hbar)*trapz(E,T.*log((1+exp((E-EF_S)/kBT))./(1+exp((E-EF_D)/kBT))));
         I = (1e-4*4*pi*meff_B*q^3*kBT/((2*pi*hbar)^3))*trapz(E,Tr.*log((1+exp((EF_S-E)/kBT))./(1+exp((EF_D-E)/kBT))));
    else
        I = 0;
    end
    
    T_TLM = Tr;
    I_TLM = I;
    

end
function [T_TMM,I_TMM] = CurrentTMM(EBD,UB,w,meff_mid,alfa_S,alfa_D)

    % I is the current tunneling from drain to source, in A
    % it is calculated from the product f(E) T(E) using Landauer Formula
    % T(E) is the transmission probibility at a certain energy E
    
    Ns = EBD.Ns;
    meff_B = EBD.m0; 
    NE = EBD.NE ;
    E = EBD.E;
    hbar = EBD.hbar;
   
    kBT = EBD.kBT;
    q = EBD.q;
    EF_D = EBD.EF_D;
    EF_S = EBD.EF_S;
    

    alfa = zeros(NE,Ns);

    for i = 1:NE
        alfa(i,:)=1e-9*sqrt(2*q*meff_mid.*(UB-E(i)))/hbar; % nm-1
    end
    % To avoid divide by zero, replace zero enteries in alfa by a small number
    alfa(alfa==0)=EBD.eps0;
    
    [Mpro11,Mpro12,Mpro21,Mpro22] =Mpro(alfa,meff_mid,w,NE,Ns);

    % We divide here meff_B by a number just for normalization
    mr = meff_B/9.10953e-31;

    WW11 = mr * exp(alfa_D*w*Ns)./(2*alfa_S) .* ...
        (alfa_S.*Mpro11(:,1)/mr + Mpro21(:,1) + alfa_D/mr.*(alfa_S/mr.*Mpro12(:,1)+Mpro22(:,1)));
    t = 1./WW11;
    T = (alfa_D./alfa_S).*abs(t).^2; % the transmission probability
    if size(E,1)>1
        %I = 2*q^2/(pi*hbar)*trapz(E,T.*log((1+exp((E-EF_S)/kBT))./(1+exp((E-EF_D)/kBT))));
        I = (1e-4*4*pi*meff_B*q^3*kBT/((2*pi*hbar)^3))*trapz(E,T.*log((1+exp((EF_S-E)/kBT))./(1+exp((EF_D-E)/kBT))));
    else
        I = 0;
    end
    
   T_TMM = T;
   I_TMM = I;

end


function  [T_WKB,I_WKB] = CurrentWKB(obj,pot_barrier,Vd,x)

    Trans=0;
    Ef=obj.Ef;
    phib=obj.Phi_B; 
    m=obj.m0;
    h=obj.hbar;
    kbT=obj.kBT;
    e=obj.q;
    xstep=0.1e-9;
    Vdstep=0.02;
    Vdmin=0.0;
    T_WKB=[];
    %Loop at each Enegy state in Caclualted tunneling Window that
    %defined at this Vds input
    Evect = obj.E;
    for index_data=1:1:length(Evect)

        E = Evect(index_data);
        % if Electron within Tunneling Window calculate the probability
        if E<Ef+phib   
            Tr=WKB_Tunneling(obj,pot_barrier,E,x);
        % otherwise there is no barrier infront of the electron and it
        % should have 100% to pass from left to right.
        else
            Tr=1;
        end
        T_WKB(index_data)=Tr;
        
    end
    figure(3);plot(E,T_WKB)
    f=T_WKB.*(4*pi*m*e*e*e*kbT/(h^3)).*log((1+exp((Ef-Evect)./kbT))./(1+exp((Ef-Evect-Vd)./kbT)));
    I_WKB=trapz(Evect,f)*obj.heights{1}*1e-9*1e-3;
    index_plot=fix((Vd-Vdmin)/Vdstep+1);
    
end

function  [T_Simmon,I_Simmon] = CurrentSimmmon(obj,Potential,xx,Vd_vect)
    
    T_Simmon=[];
    I_Simmon=[];
    
    for index=1:length(Vd_vect)

        Vd = Vd_vect(index);

        Trans=0;
        Ef = obj.EL;
        e = obj.q;
        m = obj.m0;
        %% --------------- Specify Tunneling Window -----------------------%%

        %% In Dr Shaker Code
         Estep = 0.001;
        % Ef=5.3;
        % y=0; 
        % phib=0.1;
        % L=30e-9;
        % tch=7e-9;
        % R= tch/L;
        % dg=10e-9;
        % er=60/10;		%etunnel/egate
        % xstep=0.1e-9;
        % Vdstep=0.02;
        % Vdmin=0.0;
        % Vdmax=1.5;
        % Vg = 2;
        % Vd_vect=Vdmin:Vdstep:Vdmax; 
        % 
        % lamda=L*sqrt((er)*(dg/L)*(R)+(y/L)*(R)-0.5*(y/L)^2);

        % For Simmon Tunneling 
        Emin=Ef-Vd-10*obj.kBT; 
        Emax=Ef+obj.Phi_B;  %the energy of the tunnel electrons only

        % For Therminoic
        % Emin=Ef+phib; 
        % Emax=Ef+phib+15*kbT;  %thermoionic current component only

        %the energy of both the tunnel and thermoionic current componenents  
        % Emin=Ef-Vd-10*kbT; 
        % Emax=Ef++phib+15*kbT;  

        Evect=Emin:Estep:Emax;
        %% -------------------------------------------------------------- %%

        Int2=zeros(size(Evect));
        for index_data=1:1:length(Evect)

            E = Evect(index_data);

            if E<Emax   

                %%[xx,d1,pot_barrier] =Analytical_Pot(Vg,Vd,E);
                L=sum(obj.channel_lengths) * 1e-9;
                xtun = 0:L/1000:L;
                xstep=0.1e-9;
                [minval, index_min] = min(abs(Potential)); 
                ds = xtun(index_min);
                d1=min(ds);
                %xd=0:xstep*1e9:1;
                %xx=xd*d1;

                %% ----------- Step 1 Tunneling Calc ----------%
                intfunction=sqrt(2*m*e*abs(Potential-E)/obj.hbar^2);
                Kphase=d1*1e9*trapz(xx,intfunction);
                Tr=exp(-2*Kphase);

            else
                Tr=1;
            end

            Trans(index_data)=Tr;
            %% ---------------- Step 2 Tunneling Calc ------------%
            Ex = linspace(E,Emax,1e4); 
            dEx=(max(Ex)-E)/length(Ex);
            f1 = 1./(exp((Ex-Ef)./obj.kBT)+1); % Fermi distribution function for metal 1
            f2 = 1./(exp((Ex-Ef+Vd)./obj.kBT)+1); % Fermi distribution function for metal 2
            fdif = f1 - f2;  % Difference between two fermi function after applied bias voltage 
            Int1 = sum(fdif)*dEx; % sum of 400 points for calculate the integral of fermi function(fdif)
            Int2(index_data) = Tr*Estep*Int1;

        end
        I_Simmon(index) = (4*pi*m*e*e*e/(obj.h^3))*sum(Int2);
        T_Simmon(index) = Trans(index_data) ;
    end
    
end

function [pro11,pro12,pro21,pro22] = Mpro(alfa,meff_mid,w,NE,Ns)
    % To evaluate the product of matrices Pi,Pi+1,...PN
    
    % pro11,pro22 are dimensionless
    % pro12 is in nm
    % pro21 is in nm-1

    % Initially we should multiply by identity matrix
    pro11(:,Ns+1)= ones(NE,1);
    pro12(:,Ns+1)= zeros(NE,1);
    pro21(:,Ns+1)= zeros(NE,1);
    pro22(:,Ns+1)= ones(NE,1);

    for i=Ns:-1:1 % loop on different points in the domain starting from right
        M11 = cosh(alfa(:,i)*w);
        M22 = M11;
        m = meff_mid(i)/9.10953e-31; % just for normalization of number range
        M12 = - m*sinh(alfa(:,i)*w)./alfa(:,i);    
        M21 = - alfa(:,i).*sinh(alfa(:,i)*w)/m;

        pro11(:,i)= M11 .* pro11(:,i+1) + M12 .* pro21(:,i+1);
        pro12(:,i)= M11 .* pro12(:,i+1) + M12 .* pro22(:,i+1);
        pro21(:,i)= M21 .* pro11(:,i+1) + M22 .* pro21(:,i+1);
        pro22(:,i)= M21 .* pro12(:,i+1) + M22 .* pro22(:,i+1);
    end
end
function [xx,d1,pot_barrier]=Analytical_Pot(Vg,Vd,E)
    
    global L lamda  phib Ef 
    
    R=tch/L;
    dg=10e-9;
    er=60/10;
    y=0; 
    L=sum(device.channel_lengths)
    lamda=L*sqrt((er)*(dg/L)*(R)+(y/L)*(R)-0.5*(y/L)^2);

    xtun = 0:L/1000:L;
    xstep=0.1e-9;
    
    funtun = device.Ef - (-device.Phi_B+Vg+(  (  Vd-Vg+Vg*cosh(L/lamda)) /(sinh(L/lamda))  )*sinh(xtun/lamda)-Vg*cosh(xtun/lamda) ) - E;
    [minval, index_min] = min(abs(funtun)); 
    ds = xtun(index_min);
    d1=min(ds);
    xd=0:xstep*1e9:1;
    xx=xd*d1;
    pot_barrier=device.Ef-(-device.Phi_B+Vg+((Vd-Vg+Vg*cosh(L/lamda))/(sinh(L/lamda)))*sinh(xx/lamda)-Vg*cosh(xx/lamda));
    
end
