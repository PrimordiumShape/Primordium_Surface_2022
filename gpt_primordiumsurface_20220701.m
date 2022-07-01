function m = gpt_primordiumsurface_20220701( m )
%m = gpt_primordiumsurface_20220701( m )
%   Morphogen interaction function.
%   Written at 2022-07-01 15:24:31.
%   GFtbox revision 20200318, 2020-03-18 10:00.

% The user may edit any part of this function lying between lines that
% begin "%%% USER CODE" and "%%% END OF USER CODE".  Those lines themselves
% delimiters themselves must not be moved, edited, deleted, or added.

    if isempty(m), return; end

    fprintf( 1, '%s found in %s\n', mfilename(), which(mfilename()) );

    setGlobals();
    realtime = m.globalDynamicProps.currenttime;
    dt = m.globalProps.timestep;

%%% USER CODE: INITIALISATION

% In this section you may modify the mesh in any way whatsoever.
    if (Steps(m)==0) && m.globalDynamicProps.doinit % First iteration
        % Zero out a lot of stuff to create a blank slate.  
        % If no morphogens are set in the GUI it may be useful to
        % zero some arrays by uncommenting the following.
        % m.morphogens(:) = 0;
        % m.morphogenclamp(:) = 0;
        % m.mgen_production(:) = 0;
        % m.mgen_absorption(:) = 0;
        % m.seams(:) = false;
        % m.mgen_dilution(:) = false;

        % Set up names for variant models.
        % Useful for running multiple models on a cluster.
        m.userdata.ranges.modelname.range = { 'LEAF', 'pREV_miPRS', 'FLOWER', 'pFIL_PRS', 'pLFY_amiR_TOR' };  % CLUSTER
        m.userdata.ranges.modelname.index = 5;                      
    end
    
    modelname = m.userdata.ranges.modelname.range{m.userdata.ranges.modelname.index};  % CLUSTER
    disp(sprintf('\nRunning %s model %s\n',mfilename, modelname));
    
    switch modelname
        case 'LEAF'
            % Set up the parameters (e.g. mutations) for this model here.
        case 'pREV_miPRS'
            % Set up the parameters (e.g. mutations) for this model here.
        case 'FLOWER'
            % Set up the parameters (e.g. mutations) for this model here.
        case 'pFIL_PRS'
            % Set up the parameters (e.g. mutations) for this model here.
        case 'pLFY_amiR_TOR'
            % Set up the parameters (e.g. mutations) for this model here.
        otherwise
            % If you reach here, you probably forgot a case.
    end

    % visualize different areas in SAM
    m = leaf_plotoptions( m, 'morphogen', {'V_SAM','V_PRIMORDIUM','ID_AD','ID_MID','ID_AB'});
    
    % setup the time of each running step
    m.globalProps.timestep=1;
%%% END OF USER CODE: INITIALISATION

%%% SECTION 1: ACCESSING MORPHOGENS AND TIME.
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.

% Each call of getMgenLevels below returns four results:
% XXX_i is the index of the morphogen called XXX.
% XXX_p is the vector of all of its values.
% XXX_a is its mutation level.
% XXX_l is the "effective" level of the morphogen, i.e. XXX_p*XXX_a.
% In SECTION 3 of the automatically generated code, all of the XXX_p values
% will be copied back into the mesh.

    polariser_i = FindMorphogenRole( m, 'POLARISER' );
    P = m.morphogens(:,polariser_i);
    [kapar_i,kapar_p,kapar_a,kapar_l] = getMgenLevels( m, 'KAPAR' );  %#ok<ASGLU>
    [kaper_i,kaper_p,kaper_a,kaper_l] = getMgenLevels( m, 'KAPER' );  %#ok<ASGLU>
    [kbpar_i,kbpar_p,kbpar_a,kbpar_l] = getMgenLevels( m, 'KBPAR' );  %#ok<ASGLU>
    [kbper_i,kbper_p,kbper_a,kbper_l] = getMgenLevels( m, 'KBPER' );  %#ok<ASGLU>
    [knor_i,knor_p,knor_a,knor_l] = getMgenLevels( m, 'KNOR' );  %#ok<ASGLU>
    [strainret_i,strainret_p,strainret_a,strainret_l] = getMgenLevels( m, 'STRAINRET' );  %#ok<ASGLU>
    [arrest_i,arrest_p,arrest_a,arrest_l] = getMgenLevels( m, 'ARREST' );  %#ok<ASGLU>
    [id_sam_i,id_sam_p,id_sam_a,id_sam_l] = getMgenLevels( m, 'ID_SAM' );  %#ok<ASGLU>
    [f_boundary_i,f_boundary_p,f_boundary_a,f_boundary_l] = getMgenLevels( m, 'F_BOUNDARY' );  %#ok<ASGLU>
    [v_primordium_i,v_primordium_p,v_primordium_a,v_primordium_l] = getMgenLevels( m, 'V_PRIMORDIUM' );  %#ok<ASGLU>
    [f_border_i,f_border_p,f_border_a,f_border_l] = getMgenLevels( m, 'F_BORDER' );  %#ok<ASGLU>
    [id_ad_i,id_ad_p,id_ad_a,id_ad_l] = getMgenLevels( m, 'ID_AD' );  %#ok<ASGLU>
    [id_mid_i,id_mid_p,id_mid_a,id_mid_l] = getMgenLevels( m, 'ID_MID' );  %#ok<ASGLU>
    [id_ab_i,id_ab_p,id_ab_a,id_ab_l] = getMgenLevels( m, 'ID_AB' );  %#ok<ASGLU>
    [id_g_i,id_g_p,id_g_a,id_g_l] = getMgenLevels( m, 'ID_G' );  %#ok<ASGLU>
    [v_sam_i,v_sam_p,v_sam_a,v_sam_l] = getMgenLevels( m, 'V_SAM' );  %#ok<ASGLU>
    [id_org_plus_i,id_org_plus_p,id_org_plus_a,id_org_plus_l] = getMgenLevels( m, 'ID_ORG_PLUS' );  %#ok<ASGLU>
    [id_org_minus_i,id_org_minus_p,id_org_minus_a,id_org_minus_l] = getMgenLevels( m, 'ID_ORG_MINUS' );  %#ok<ASGLU>
    [id_inh_i,id_inh_p,id_inh_a,id_inh_l] = getMgenLevels( m, 'ID_INH' );  %#ok<ASGLU>
    [f_inhibit_i,f_inhibit_p,f_inhibit_a,f_inhibit_l] = getMgenLevels( m, 'F_INHIBIT' );  %#ok<ASGLU>
    [f_polarity_i,f_polarity_p,f_polarity_a,f_polarity_l] = getMgenLevels( m, 'F_POLARITY' );  %#ok<ASGLU>
    [id_nopolar_i,id_nopolar_p,id_nopolar_a,id_nopolar_l] = getMgenLevels( m, 'ID_NOPOLAR' );  %#ok<ASGLU>
    [id_org_neg_i,id_org_neg_p,id_org_neg_a,id_org_neg_l] = getMgenLevels( m, 'ID_ORG_NEG' );  %#ok<ASGLU>

% Mesh type: circle
%            asym: 0
%          centre: 0
%       circumpts: 120
%       coneangle: 0
%         dealign: 0
%       generalFE: 0
%          height: 0
%        innerpts: 0
%          layers: 0
%             new: 1
%      randomness: 0.1
%           rings: 30
%      semicircle: 0
%       thickness: 0
%          xwidth: 2
%          ywidth: 2

%            Morphogen    Diffusion   Decay   Dilution   Mutant
%            --------------------------------------------------
%                KAPAR         ----    ----       ----     ----
%                KAPER         ----    ----       ----     ----
%                KBPAR         ----    ----       ----     ----
%                KBPER         ----    ----       ----     ----
%                 KNOR         ----    ----       ----     ----
%            POLARISER        0.005    ----       ----     ----
%            STRAINRET         ----    ----       ----     ----
%               ARREST         ----    ----       ----     ----
%               ID_SAM         ----    ----       ----     ----
%           F_BOUNDARY         ----    ----       ----     ----
%         V_PRIMORDIUM         ----    ----       ----     ----
%             F_BORDER         ----    ----       ----     ----
%                ID_AD         ----    ----       ----     ----
%               ID_MID         ----    ----       ----     ----
%                ID_AB         ----    ----       ----     ----
%                 ID_G         ----    ----       ----     ----
%                V_SAM         ----    ----       ----     ----
%          ID_ORG_PLUS         ----    ----       ----     ----
%         ID_ORG_MINUS         ----    ----       ----     ----
%               ID_INH         ----    ----       ----     ----
%            F_INHIBIT         ----    ----       ----     ----
%           F_POLARITY         ----    ----       ----     ----
%           ID_NOPOLAR         ----    ----       ----     ----
%           ID_ORG_NEG         ----    ----       ----     ----


%%% USER CODE: MORPHOGEN INTERACTIONS

% In this section you may modify the mesh in any way that does not
% alter the set of nodes.

    if (Steps(m)==0) && m.globalDynamicProps.doinit  % Initialisation code.
        
        % define factor to form SAM with the shape of hemisphere
        id_sam_p = sqrt(m.nodes(:,1).^2 + m.nodes(:,2).^2);
        id_sam_p = ((max(id_sam_p)- id_sam_p))/max(id_sam_p); % normalise
        id_sam_l = id_sam_p * id_sam_a;
        
        % visualize the entire SAM area
        v_sam_p(sqrt(m.nodes(:,1).^2+m.nodes(:,2).^2) < 0.9) = 1;
        
        % define factors to identify and visualize primordium 
        ind_radial = sqrt(m.nodes(:,1).^2 + m.nodes(:,2).^2);
        
        switch modelname
            case {'LEAF','pREV_miPRS'}
                % define factors related with primordium area
                ind_primordium_disc = find((m.nodes(:,1) > 0.4) & (ind_radial <= 0.85));  
                v_primordium_p(ind_primordium_disc) = 1; % define primordium area
                
                ind_inhibit_line = find((v_primordium_p == 1) & (m.nodes(:,1) < 0.45)); 
                f_inhibit_p(ind_inhibit_line) = 1; % define the boundary of SAM and adaxial primordium
                
                % define factors to produce polarity field 
                ind_nopolar_area = find((m.nodes(:,1) < 0.35) | (ind_radial > 0.9)); 
                id_nopolar_p(ind_nopolar_area) = 1; % define the area outside the primordium
               
                ind_polarity = find (v_primordium_p==1 & ((m.nodes(:,1) + 0.77).^2 + (m.nodes(:,2)).^2 > 1.42.^2) & ((m.nodes(:,1) + 0.77).^2 + (m.nodes(:,2)).^2 < 1.46.^2)); 
                f_polarity_p(ind_polarity) = 1;
                id_org_minus_p( f_polarity_p == 1 ) = 1; % define org- of leaf primordium
                            
                ind_boundary_circle = find((v_primordium_p == 1) & ((m.nodes(:,1) < 0.45) | (ind_radial > 0.81))); 
                f_boundary_p(ind_boundary_circle) = 1;
                id_org_plus_p (f_boundary_p == 1) = 1; % define org+ of leaf primordium
                               
                % define three domains of leaf primordium 
                ind_ad_area = find (v_primordium_p==1 & ((m.nodes(:,1) + 0.77).^2 + (m.nodes(:,2)).^2 < 1.38.^2));
                id_ad_p (ind_ad_area) = 1;
                    
                ind_mid_area = find (v_primordium_p==1 & ((m.nodes(:,1) + 0.77).^2 + (m.nodes(:,2)).^2 > 1.38.^2) & ((m.nodes(:,1) + 0.77).^2 + (m.nodes(:,2)).^2 < 1.48.^2));
                id_mid_p (ind_mid_area) = 1;
                
                ind_ab_area = find (v_primordium_p==1 &((m.nodes(:,1) + 0.77).^2 + (m.nodes(:,2)).^2 > 1.48.^2));  
                id_ab_p (ind_ab_area) = 1;  
            case {'FLOWER'}
                % define factors related with primordium area
                ind_primordium_disc = find((m.nodes(:,1) > 0.4) & (ind_radial <= 0.85)); 
                v_primordium_p(ind_primordium_disc) = 1;  % define primordium area
                             
                ind_boundary_circle = find(((m.nodes(:,1) < 0.4) & (m.nodes(:,1) > 0.35) & (ind_radial <= 0.9)) | ((m.nodes(:,1) > 0.4) & (ind_radial > 0.85) & (ind_radial <= 0.9)));
                f_boundary_p(ind_boundary_circle) = 1; % define the boundary of SAM and entire primordium
                
                % define factors to produce polarity field there
                id_org_minus_p(sqrt((m.nodes(:,1)-0.52).^2 + m.nodes(:,2).^2)  < 0.02) = 1;  % define the adaxial org- of flower primordium
                id_org_plus_p (f_boundary_p == 1) = 1;  % define org+ of flower primordium                       
                
                ind_nopolar_area = find((m.nodes(:,1) < 0.3) | (ind_radial > 0.95)); 
                id_nopolar_p(ind_nopolar_area) = 1; % define the area outside the primordium
                              
                % define two domains of flower primordium
                ind_ad_area = find (v_primordium_p==1 & ((m.nodes(:,1)).^2 + (m.nodes(:,2)).^2 < 0.72.^2));
                id_ad_p (ind_ad_area) = 1;
                        
                ind_ab_area = find (v_primordium_p==1 &((m.nodes(:,1)).^2  + (m.nodes(:,2)).^2 > 0.72.^2));  
                id_ab_p (ind_ab_area) = 1;
              case {'pFIL_PRS','pLFY_amiR_TOR'}
                % define factors related with primordium area
                ind_primordium_disc = find((m.nodes(:,1) > 0.4) & (ind_radial <= 0.85)); 
                v_primordium_p(ind_primordium_disc) = 1;  % define primordium area
                             
                ind_boundary_circle = find(((m.nodes(:,1) < 0.4) & (m.nodes(:,1) > 0.35) & (ind_radial <= 0.9)) | ((m.nodes(:,1) > 0.4) & (ind_radial > 0.85) & (ind_radial <= 0.9)));
                f_boundary_p(ind_boundary_circle) = 1; % define the boundary of SAM and entire primordium
                
                % define factors to produce polarity field there
                id_org_minus_p(sqrt((m.nodes(:,1)-0.52).^2 + m.nodes(:,2).^2)  < 0.02) = 1;  % define the adaxial org- of flower primordium
                id_org_neg_p(sqrt((m.nodes(:,1)-0.78).^2 + m.nodes(:,2).^2)  < 0.02) = 1;  % define org- of flower primordium in transgenic plants
                id_org_plus_p (f_boundary_p == 1) = 1;  % define org+ of flower primordium                       
                
                ind_nopolar_area = find((m.nodes(:,1) < 0.3) | (ind_radial > 0.95)); 
                id_nopolar_p(ind_nopolar_area) = 1; % define the area outside the primordium
                              
                % define two domains of flower primordium
                ind_ad_area = find (v_primordium_p==1 & ((m.nodes(:,1)).^2 + (m.nodes(:,2)).^2 < 0.72.^2));
                id_ad_p (ind_ad_area) = 1;
                        
                ind_ab_area = find (v_primordium_p==1 &((m.nodes(:,1)).^2  + (m.nodes(:,2)).^2 > 0.72.^2));  
                id_ab_p (ind_ab_area) = 1;
            otherwise
                % If you reach here, you probably forgot a case.
        end
        
        % fixing x, y and z vertices for the base to prevent base from moving up or down
        f_border_p(ind_radial > 0.9) = 1;
        m = leaf_fix_vertex(m,'vertex',find(f_border_p==1),'dfs'); 
    end
   
 
    % define different phases according to growth time
    if (realtime < 10)  % SAM formation phase
        f_samformation_p = 0.15;
    else
        f_samformation_p = 0;
    end
    
    if (realtime >= 10) % promordium growth phase
        f_development_p = 0.1;
    else
        f_development_p = 0;
    end


    % Code common to all models.
    % @@PRN Polariser Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
    % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
    % @@KRN Growth Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.

    % Code for specific models.
    switch modelname
        case 'LEAF'  
            % @@PRN Polariser Regulatory Network 
            if (Steps(m)==0)
                P(id_org_plus_p == 1) = 1;
                P(id_org_minus_p == 1) = 0;
                P(id_nopolar_p == 1) = 0;
                m.morphogenclamp((id_org_plus_p == 1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_org_minus_p==1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_nopolar_p == 1),polariser_i) = 1; % double(sourcenodes);
            end
            m = leaf_mgen_conductivity( m, 'POLARISER', 0.002);  % diffusion rate
            m = leaf_mgen_absorption( m, 'POLARISER', 0);  % degradation rate
            m = leaf_setproperty(m,'mingradient',0.1); % i.e. threshold for using polariser gradient
            
            % @@GRN Gene Regulatory Network
            d_sym = abs(m.nodes(:,2));
            gradient_1 = 1 .* (1 - d_sym);
            gradient_2 = 10 .* (0.1 - d_sym .^2);

            d_gra = sqrt((m.nodes(:,1).^2 + m.nodes(:,2).^2));
            line_gradient = (d_gra - 0.3).^3;
           
            if (realtime <= 10)
                kapar_p = f_samformation_p .* 0.4 .* id_sam_l;
                kaper_p = kapar_p ; 
            elseif (10 < realtime && realtime <= 20) % stage 1
                id_g_p (find(v_primordium_p==1)) = 2;
                id_g_l = id_g_p .* id_g_a;
            
                id_mid_l = (1) .* (gradient_1 + gradient_2).* id_mid_p .* id_mid_a;
                id_ad_l = (-1.9) .* id_ad_p .* id_ad_a;
                id_ab_l = (-1.5) .* id_ab_p .* id_ab_a;   
                
                id_inh_p(f_inhibit_p==1)= -3;
                id_inh_l = id_inh_p .* id_inh_a;
                
                kapar_p = 0.3 .* f_development_p .* (id_g_l + id_mid_l + id_ad_l + id_ab_l + id_inh_l);
                kaper_p = kapar_p; 
            elseif (20 < realtime && realtime <= 30)  % stage 2
                id_g_p (find(v_primordium_p==1)) = 1.2; 
                id_g_l = id_g_p .* id_g_a; 
                
                id_inh_p (f_inhibit_p==1) = -2;
                id_inh_l = id_inh_p .* id_inh_a;
                
                kapar_p = f_development_p .* (line_gradient .* id_g_l + id_inh_l);
                kaper_p = 0.3 .* f_development_p .* id_g_l ; 
            else  % stage 3
                id_g_p (find(v_primordium_p==1)) = 0.5; 
                id_g_l = id_g_p .* id_g_a; 
                
                id_inh_p (f_inhibit_p==1) = -1;
                id_inh_l = id_inh_p .* id_inh_a;
                
                kapar_p = f_development_p .* (line_gradient .* id_g_l + id_inh_l);
                kaper_p = 0.2 .* f_development_p .* id_g_l ;  
            end

            % @@KRN Growth Regulatory Network     
            kbpar_p = kapar_p; 
            kbper_p = kaper_p; 
            knor_p  = 0; 
        case 'pREV_miPRS'  
            % @@PRN Polariser Regulatory Network 
            if (Steps(m)==0)
                P(id_org_plus_p == 1) = 1;
                P(id_org_minus_p == 1) = 0;
                P(id_nopolar_p == 1) = 0;
                m.morphogenclamp((id_org_plus_p == 1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_org_minus_p==1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_nopolar_p == 1),polariser_i) = 1; % double(sourcenodes);
            end
            m = leaf_mgen_conductivity( m, 'POLARISER', 0.002);  % diffusion rate
            m = leaf_mgen_absorption( m, 'POLARISER', 0);  % degradation rate
            m = leaf_setproperty(m,'mingradient',0.1); % i.e. threshold for using polariser gradient
            
            % @@GRN Gene Regulatory Network       
            if (realtime <= 10)
                kapar_p = f_samformation_p .* 0.4 .* id_sam_l;
                kaper_p = kapar_p ;          
            else
                id_g_p (find(v_primordium_p==1)) = 1;
                id_g_l = id_g_p .* id_g_a;
            
                id_mid_l = (0) .* id_mid_p .* id_mid_a;
                id_ad_l = (0).* id_ad_p .* id_ad_a;
                id_ab_l = (-0.5) .* id_ab_p .* id_ab_a;   
                
                id_inh_p(f_inhibit_p==1)= -3;
                id_inh_l = id_inh_p .* id_inh_a;
                
                kapar_p = 0.5 .* f_development_p .* (id_g_l + id_mid_l + id_ad_l + id_ab_l + id_inh_l);
                kaper_p = 0.2 .* f_development_p .* (id_g_l + id_inh_l + id_ab_l);    
            end

            % @@KRN Growth Regulatory Network     
            kbpar_p = kapar_p; 
            kbper_p = kaper_p; 
            knor_p  = 0;
        case 'FLOWER'  
            % @@PRN Polariser Regulatory Network     
            if (Steps(m)==0) 
                P(id_org_plus_p == 1) = 1;
                P(id_org_minus_p == 1) = 0;
                P(id_nopolar_p == 1) = 0;
                m.morphogenclamp((id_org_plus_p == 1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_org_minus_p==1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_nopolar_p == 1),polariser_i) = 1; % double(sourcenodes);
            end
        
            m = leaf_mgen_conductivity( m, 'POLARISER', 0.002);  % diffusion rate
            m = leaf_mgen_absorption( m, 'POLARISER', 0);  % degradation rate
            m = leaf_setproperty(m,'mingradient',0.1); % i.e. threshold for using polariser gradient

            % @@GRN Gene Regulatory Network                    
            if (realtime <= 10)
                kapar_p = f_samformation_p .* 0.4 .* id_sam_l;
                kaper_p = kapar_p ; 
            else
                id_g_p (v_primordium_p==1) = 0.5;
                id_g_l = id_g_p .* id_g_a;
                
                id_inh_p (f_boundary_p==1) = -1.5;
                id_inh_l = id_inh_p .* id_inh_a;
                
                id_ad_l = (0.2).* id_ad_p .* id_ad_a;
                id_ab_l = (0.1).* id_ab_p .* id_ab_a;   
                
                kapar_p = 0.4 .* f_development_p .* (id_g_l + id_ad_l + id_ab_l+ id_inh_l);
                kaper_p = 0.4 .* f_development_p .* (0.3 .* id_g_l + 1.5 .* id_ad_l + 2.5 .* id_ab_l+ id_inh_l);
            end
            
            % @@KRN Growth Regulatory Network
            kbpar_p = kapar_p; 
            kbper_p = kaper_p;  
            knor_p  = 0;   
        case 'pFIL_PRS'  % Flower transgene line 1                 
            % @@PRN Polariser Regulatory Network                             
            if (Steps(m)==0 ) 
                P(id_org_plus_p == 1) = 1;
                P(id_org_minus_p == 1) = 0;
                P(id_nopolar_p == 1) = 0;
                m.morphogenclamp((id_org_plus_p == 1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_org_minus_p==1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_nopolar_p == 1),polariser_i) = 1; % double(sourcenodes);
            end
            
            if (Steps(m)> 19) % add a second polarity sink point at stage 2 of primordium growth phase 
                P(id_org_plus_p == 1) = 1;
                P(id_org_minus_p == 1) = 0;
                P(id_org_neg_p == 1) = 0;
                P(id_nopolar_p == 1) = 0;
                m.morphogenclamp((id_org_plus_p == 1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_org_minus_p==1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_org_neg_p==1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_nopolar_p == 1),polariser_i) = 1; % double(sourcenodes);
            end
            
            m = leaf_mgen_conductivity( m, 'POLARISER', 0.005);  % diffusion rate
            m = leaf_mgen_absorption( m, 'POLARISER', 0);  % degradation rate
            m = leaf_setproperty(m,'mingradient',0.1); % i.e. threshold for using polariser gradient

            % @@GRN Gene Regulatory Network
            if (realtime <= 10)
                kapar_p = f_samformation_p .* 0.4 .* id_sam_l;
                kaper_p = kapar_p ; 
            elseif (10 < realtime && realtime <= 19)  % stage 1
                id_g_p (v_primordium_p==1) = 0.5;
                id_g_l = id_g_p .* id_g_a;
                
                id_inh_p (f_boundary_p==1) = -1.5;
                id_inh_l = id_inh_p .* id_inh_a;
                
                id_ad_l = (0.2).* id_ad_p .* id_ad_a;
                id_ab_l = (0.1).* id_ab_p .* id_ab_a;   
                
                kapar_p = 0.4 .* f_development_p .* (id_g_l + id_ad_l + id_ab_l+ id_inh_l);
                kaper_p = 0.4 .* f_development_p .* (0.3 .* id_g_l + 1.5 .* id_ad_l + 2.5 .* id_ab_l+ id_inh_l);
            else  % stage 2
                id_g_p (v_primordium_p==1) = 0.3;
                id_g_l = id_g_p .* id_g_a;
                
                id_inh_p (f_boundary_p==1) = -1.5;
                id_inh_l = id_inh_p .* id_inh_a;
                
                id_ad_l = (0.2).* id_ad_p .* id_ad_a;
                id_ab_l = (1.5).* id_ab_p .* id_ab_a;   
                
                kapar_p = 0.4 .* f_development_p .* ( id_g_l + id_ad_l + (1.2 - P.^3) .* id_ab_l+ id_inh_l);
                kaper_p = kapar_p ;
            end
            
            % @@KRN Growth Regulatory Network
            kbpar_p = kapar_p; 
            kbper_p = kaper_p;  
            knor_p  = 0;  
       case 'pLFY_amiR_TOR'  % Flower transgene line 2
            % @@PRN Polariser Regulatory Network 
            if (Steps(m)==0)
                P(id_org_plus_p == 1) = 1;
                P(id_org_minus_p == 1) = 0;
                P(id_org_neg_p == 1) = 0;
                P(id_nopolar_p == 1) = 0;
                m.morphogenclamp((id_org_plus_p == 1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_org_minus_p==1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_org_neg_p==1),polariser_i) = 1; % double(sourcenodes);
                m.morphogenclamp((id_nopolar_p == 1),polariser_i) = 1; % double(sourcenodes);
            end
            
            m = leaf_mgen_conductivity( m, 'POLARISER', 0.005);  % diffusion rate
            m = leaf_mgen_absorption( m, 'POLARISER', 0);  % degradation rate
            m = leaf_setproperty(m,'mingradient',0.1); % i.e. threshold for using polariser gradient

            % @@GRN Gene Regulatory Network: normal 
            if (realtime <= 10)
                kapar_p = f_samformation_p .* 0.4 .* id_sam_l;
                kaper_p = kapar_p ; 
            else
                id_g_p (v_primordium_p==1) = 0.5;
                id_g_l = id_g_p .* id_g_a;
                
                id_inh_p (f_boundary_p==1) = -1;
                id_inh_l = id_inh_p .* id_inh_a;
                
                id_ad_l = (0.1).* id_ad_p .* id_ad_a;
                id_ab_l = (0.4).* id_ab_p .* id_ab_a;  % flatten
%                 id_ab_l = (1).* id_ab_p .* id_ab_a;  % bulge
                
                kapar_p = 0.3 .* f_development_p .* (id_g_l + id_ad_l + (1.3 - P).* id_ab_l + id_inh_l);
                kaper_p = kapar_p ;
            end
            
            % @@KRN Growth Regulatory Network
            kbpar_p = kapar_p; 
            kbper_p = kaper_p;  
            knor_p  = 0;  
        otherwise
            % If this happens, maybe you forgot a model.
    end
%%% END OF USER CODE: MORPHOGEN INTERACTIONS

%%% SECTION 3: INSTALLING MODIFIED VALUES BACK INTO MESH STRUCTURE
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.
    m.morphogens(:,polariser_i) = P;
    m.morphogens(:,kapar_i) = kapar_p;
    m.morphogens(:,kaper_i) = kaper_p;
    m.morphogens(:,kbpar_i) = kbpar_p;
    m.morphogens(:,kbper_i) = kbper_p;
    m.morphogens(:,knor_i) = knor_p;
    m.morphogens(:,strainret_i) = strainret_p;
    m.morphogens(:,arrest_i) = arrest_p;
    m.morphogens(:,id_sam_i) = id_sam_p;
    m.morphogens(:,f_boundary_i) = f_boundary_p;
    m.morphogens(:,v_primordium_i) = v_primordium_p;
    m.morphogens(:,f_border_i) = f_border_p;
    m.morphogens(:,id_ad_i) = id_ad_p;
    m.morphogens(:,id_mid_i) = id_mid_p;
    m.morphogens(:,id_ab_i) = id_ab_p;
    m.morphogens(:,id_g_i) = id_g_p;
    m.morphogens(:,v_sam_i) = v_sam_p;
    m.morphogens(:,id_org_plus_i) = id_org_plus_p;
    m.morphogens(:,id_org_minus_i) = id_org_minus_p;
    m.morphogens(:,id_inh_i) = id_inh_p;
    m.morphogens(:,f_inhibit_i) = f_inhibit_p;
    m.morphogens(:,f_polarity_i) = f_polarity_p;
    m.morphogens(:,id_nopolar_i) = id_nopolar_p;
    m.morphogens(:,id_org_neg_i) = id_org_neg_p;

%%% USER CODE: FINALISATION

% In this section you may modify the mesh in any way whatsoever.

    % If needed force FE to subdivide (increase number FE's) here
    % if realtime==280+dt
         % m = leaf_subdivide( m, 'morphogen','id_vent',...
         %       'min',0.5,'max',1,...
         %       'mode','mid','levels','all');
    % end
% Cut the mesh along the seams (see above)
    % if m.userdata.CutOpen==1
    %    m=leaf_dissect(m);
    %    m.userdata.CutOpen=2;        
    %    Relax accumulated stresses slowly i.e. 0.95 to 0.999
    %    m = leaf_setproperty( m, 'freezing', 0.999 );
    % end
%%% END OF USER CODE: FINALISATION

end


%%% USER CODE: SUBFUNCTIONS

function m = local_setproperties( m )
% This function is called at time zero in the INITIALISATION section of the
% interaction function.  It provides commands to set each of the properties
% that are contained in m.globalProps.  Uncomment whichever ones you would
% like to set yourself, and put in whatever value you want.
%
% Some of these properties are for internal use only and should never be
% set by the user.  At some point these will be moved into a different
% component of m, but for the present, just don't change anything unless
% you know what it is you're changing.

%    m = leaf_setproperty( m, 'trinodesvalid', true );
%    m = leaf_setproperty( m, 'prismnodesvalid', true );
%    m = leaf_setproperty( m, 'thresholdsq', 0.007521 );
%    m = leaf_setproperty( m, 'lengthscale', 2.000000 );
%    m = leaf_setproperty( m, 'initialArea', 3.140157 );
%    m = leaf_setproperty( m, 'bendunitlength', 1.772049 );
%    m = leaf_setproperty( m, 'thicknessRelative', 0.100000 );
%    m = leaf_setproperty( m, 'thicknessArea', 1.000000 );
%    m = leaf_setproperty( m, 'hybridMesh', false );
%    m = leaf_setproperty( m, 'thicknessMode', 'physical' );
%    m = leaf_setproperty( m, 'activeGrowth', 1.000000 );
%    m = leaf_setproperty( m, 'displayedGrowth', 1.000000 );
%    m = leaf_setproperty( m, 'displayedMulti', [] );
%    m = leaf_setproperty( m, 'allowNegativeGrowth', true );
%    m = leaf_setproperty( m, 'usePrevDispAsEstimate', true );
%    m = leaf_setproperty( m, 'perturbInitGrowthEstimate', 0.000010 );
%    m = leaf_setproperty( m, 'perturbRelGrowthEstimate', 0.010000 );
%    m = leaf_setproperty( m, 'perturbDiffusionEstimate', 0.000100 );
%    m = leaf_setproperty( m, 'resetRand', false );
%    m = leaf_setproperty( m, 'mingradient', 0.000000 );
%    m = leaf_setproperty( m, 'relativepolgrad', false );
%    m = leaf_setproperty( m, 'usefrozengradient', true );
%    m = leaf_setproperty( m, 'userpolarisation', false );
%    m = leaf_setproperty( m, 'twosidedpolarisation', false );
%    m = leaf_setproperty( m, 'splitmargin', 1.400000 );
%    m = leaf_setproperty( m, 'splitmorphogen', '' );
%    m = leaf_setproperty( m, 'thresholdmgen', 0.500000 );
%    m = leaf_setproperty( m, 'bulkmodulus', 1.000000 );
%    m = leaf_setproperty( m, 'unitbulkmodulus', true );
%    m = leaf_setproperty( m, 'poissonsRatio', 0.300000 );
%    m = leaf_setproperty( m, 'starttime', 0.000000 );
%    m = leaf_setproperty( m, 'timestep', 0.010000 );
%    m = leaf_setproperty( m, 'timeunitname', '' );
%    m = leaf_setproperty( m, 'distunitname', 'mm' );
%    m = leaf_setproperty( m, 'scalebarvalue', 0.000000 );
%    m = leaf_setproperty( m, 'validateMesh', true );
%    m = leaf_setproperty( m, 'rectifyverticals', false );
%    m = leaf_setproperty( m, 'allowSplitLongFEM', true );
%    m = leaf_setproperty( m, 'allowSplitThinFEM', false );
%    m = leaf_setproperty( m, 'splitthinness', 10.000000 );
%    m = leaf_setproperty( m, 'longSplitThresholdPower', 0.000000 );
%    m = leaf_setproperty( m, 'allowSplitBentFEM', false );
%    m = leaf_setproperty( m, 'allowSplitBio', true );
%    m = leaf_setproperty( m, 'allowFlipEdges', false );
%    m = leaf_setproperty( m, 'allowElideEdges', true );
%    m = leaf_setproperty( m, 'mincellangle', 0.200000 );
%    m = leaf_setproperty( m, 'mincellrelarea', 0.040000 );
%    m = leaf_setproperty( m, 'alwaysFlat', 0.000000 );
%    m = leaf_setproperty( m, 'flattenforceconvex', true );
%    m = leaf_setproperty( m, 'flatten', false );
%    m = leaf_setproperty( m, 'flattenratio', 1.000000 );
%    m = leaf_setproperty( m, 'useGrowthTensors', false );
%    m = leaf_setproperty( m, 'useMorphogens', true );
%    m = leaf_setproperty( m, 'plasticGrowth', false );
%    m = leaf_setproperty( m, 'maxFEcells', 0 );
%    m = leaf_setproperty( m, 'inittotalcells', 0 );
%    m = leaf_setproperty( m, 'bioApresplitproc', '' );
%    m = leaf_setproperty( m, 'bioApostsplitproc', '' );
%    m = leaf_setproperty( m, 'maxBioAcells', 0 );
%    m = leaf_setproperty( m, 'biosplitarea', 0.000000 );
%    m = leaf_setproperty( m, 'biosplitarrestmgen', 'ARREST' );
%    m = leaf_setproperty( m, 'biosplitarrestmgenthreshold', 0.990000 );
%    m = leaf_setproperty( m, 'bioMinEdgeLength', 0.000000 );
%    m = leaf_setproperty( m, 'bioSpacePullInRatio', 0.100000 );
%    m = leaf_setproperty( m, 'colors', (6 values) );
%    m = leaf_setproperty( m, 'colorvariation', 0.050000 );
%    m = leaf_setproperty( m, 'colorparams', (12 values) );
%    m = leaf_setproperty( m, 'biocolormode', 'auto' );
%    m = leaf_setproperty( m, 'userpostiterateproc', [] );
%    m = leaf_setproperty( m, 'canceldrift', false );
%    m = leaf_setproperty( m, 'mgen_interaction', '' );
%    m = leaf_setproperty( m, 'mgen_interactionName', 'gpt_circle_20200821' );
%    m = leaf_setproperty( m, 'allowInteraction', true );
%    m = leaf_setproperty( m, 'interactionValid', true );
%    m = leaf_setproperty( m, 'gaussInfo', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'D', (36 values) );
%    m = leaf_setproperty( m, 'C', (36 values) );
%    m = leaf_setproperty( m, 'G', (6 values) );
%    m = leaf_setproperty( m, 'solver', 'cgs' );
%    m = leaf_setproperty( m, 'solverprecision', 'double' );
%    m = leaf_setproperty( m, 'solvertolerance', 0.001000 );
%    m = leaf_setproperty( m, 'solvertolerancemethod', 'max' );
%    m = leaf_setproperty( m, 'diffusiontolerance', 0.000010 );
%    m = leaf_setproperty( m, 'allowsparse', true );
%    m = leaf_setproperty( m, 'maxsolvetime', 1000.000000 );
%    m = leaf_setproperty( m, 'cgiters', 0 );
%    m = leaf_setproperty( m, 'simsteps', 0 );
%    m = leaf_setproperty( m, 'stepsperrender', 0 );
%    m = leaf_setproperty( m, 'growthEnabled', true );
%    m = leaf_setproperty( m, 'diffusionEnabled', true );
%    m = leaf_setproperty( m, 'flashmovie', false );
%    m = leaf_setproperty( m, 'makemovie', false );
%    m = leaf_setproperty( m, 'moviefile', '' );
%    m = leaf_setproperty( m, 'codec', 'Motion JPEG AVI' );
%    m = leaf_setproperty( m, 'autonamemovie', true );
%    m = leaf_setproperty( m, 'overwritemovie', false );
%    m = leaf_setproperty( m, 'framesize', [] );
%    m = leaf_setproperty( m, 'mov', [] );
%    m = leaf_setproperty( m, 'boingNeeded', false );
%    m = leaf_setproperty( m, 'defaultinterp', 'min' );
%    m = leaf_setproperty( m, 'readonly', false );
%    m = leaf_setproperty( m, 'projectdir', 'C:\Users\APG\GFtbox_Projects' );
%    m = leaf_setproperty( m, 'modelname', 'GPT_Circle_20200821' );
%    m = leaf_setproperty( m, 'allowsave', true );
%    m = leaf_setproperty( m, 'addedToPath', false );
%    m = leaf_setproperty( m, 'bendsplit', 0.300000 );
%    m = leaf_setproperty( m, 'usepolfreezebc', false );
%    m = leaf_setproperty( m, 'dorsaltop', true );
%    m = leaf_setproperty( m, 'defaultazimuth', -45.000000 );
%    m = leaf_setproperty( m, 'defaultelevation', 33.750000 );
%    m = leaf_setproperty( m, 'defaultroll', 0.000000 );
%    m = leaf_setproperty( m, 'defaultViewParams', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'comment', '' );
%    m = leaf_setproperty( m, 'legendTemplate', '%T: %q\n%m' );
%    m = leaf_setproperty( m, 'bioAsplitcells', true );
%    m = leaf_setproperty( m, 'bioApullin', 0.142857 );
%    m = leaf_setproperty( m, 'bioAfakepull', 0.202073 );
%    m = leaf_setproperty( m, 'viewrotationstart', -45.000000 );
%    m = leaf_setproperty( m, 'viewrotationperiod', 0.000000 );
%    m = leaf_setproperty( m, 'interactive', false );
%    m = leaf_setproperty( m, 'coderevision', 0 );
%    m = leaf_setproperty( m, 'coderevisiondate', '' );
%    m = leaf_setproperty( m, 'modelrevision', 0 );
%    m = leaf_setproperty( m, 'modelrevisiondate', '' );
%    m = leaf_setproperty( m, 'savedrunname', '' );
%    m = leaf_setproperty( m, 'savedrundesc', '' );
%    m = leaf_setproperty( m, 'vxgrad', (108 values) );
end

% Here you may write any functions of your own, that you want to call from
% the interaction function, but never need to call from outside it.
% Remember that they do not have access to any variables except those
% that you pass as parameters, and cannot change anything except by
% returning new values as results.
% Whichever section they are called from, they must respect the same
% restrictions on what modifications they are allowed to make to the mesh.

% For example:

% function m = do_something( m )
%   % Change m in some way.
% end

% Call it from the main body of the interaction function like this:
%       m = do_something( m );
