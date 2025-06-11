%1. load your SDSU model using the function load()
load("/Users/nhinguyen/Desktop/Z.mobilis/Models/SDSU.mat")
%2. Change the name of the just loaded SDSU model to "model_sdsu"
model_sdsu=SDSU;

%3. Load 478 model using the function load()
load('i_ZM478.mat')
model_478= iZM4_478;
%4. Change the name of the just loaded 478 model to "model_478"

%5. Load the excel file with your metabolites corrected (unification of metabolites of 478
%   model IDs to SDSU IDs) using [~,~,metData]=xlsread()
cd("/Users/nhinguyen/Desktop/Z.mobilis")
[~,metData,~]= xlsread("ModelsMet2_unified.xlsx")
%6. Change the metabolite IDs in 478_model with the new annotation (SDSU
%   annotation) in the excel file using the following code as template:

OG_mets=metData(:,2);
new_mets=metData(:,3);
problematic_mets={};
for i=1:length(OG_mets)
    
    new_met=new_mets{i,1};
    if isempty(new_met)
        continue
    end
    
    OG_met=OG_mets{i,1};
    OG_met_ID_478=findMetIDs(model_478,OG_met);
    if OG_met_ID_478 == 0
        problematic_mets=vertcat(problematic_mets,OG_met);
        continue
    else
        model_478.mets{OG_met_ID_478}=new_met;
    end
end

%save('/Users/nhinguyen/Desktop/Z.mobilis/Models/iZM_478_new_march_17_2025', 'model_478')
%save('/Users/nhinguyen/Desktop/Z.mobilis/Models/iZM_478_new_march_19_2025', 'model_478')

save('/Users/nhinguyen/Desktop/Z.mobilis/Models/iZM_478_new_march_21_2025', 'model_478')

%% Adding rnxs from iZM478 model to SDSU model

%1. Load SDSU model
load("/Users/nhinguyen/Desktop/Z.mobilis/Models/SDSU.mat")
%2. Change the name to model_sdsu
model_sdsu=SDSU;
%3. Load iZM478_march_19_2025 model (the new one, 17 is the one we did but
% i want to try my model to see)
load("/Users/nhinguyen/Desktop/Z.mobilis/Models/iZM_478_new_march_21_2025.mat")
%4. Change the name to model_478
% model_478= iZM_478_new_march_19_2025; No need, it already model_478

%6. load the excel file using [~,~,data]=xlsread("")
cd("/Users/nhinguyen/Desktop/Z.mobilis")
[~,~,metData]= xlsread("SDSU_ID_gene_478.xlsx")
%7. Create a new variable with the rxns IDS
rxns_IDs=metData(:,1);

%8. Use the function [NewModel]=takeRxnstoModel2(model_sdsu,model_478,rxns_IDs);
NewModel = takeRxnstoModel2(model_sdsu, model_478, rxns_IDs);
model_sdsu=NewModel;

%1. Save the new SDSU model

%save('/Users/nhinguyen/Desktop/Z.mobilis/Models/SDSUmodel_march_28_2025', 'model_sdsu');
%% Add the genes of the new reactions to the SDSU model


findRxnIDs(model_sdsu,rxn) %For finding the index of the reactions in the model_sdsu model
model_sdsu.grRules % This is the field of the model in which you need to add the grRules from the excel file



rxn_index=findRxnIDs(model_sdsu,rxn); %rxn is the reaction ID that comes from the excel file, for example '2DGLCNR'
model_sdsu.grRules{rxn_index}=new_grRule % this new_grRule comes from the "gene" column of the excel file called SDSU_ID_gene_478

%Notes
% Add an if that check if the rxn_index is empty. If is empty creates a new
% variable called problematicRxns in which you store the rxns that were not
% found. 

%Purpose: adding gene rules (grTules) to model_sdsu using excel file
%SDSU_ID_gene_478.xlsx

%from excel
%1. take the abbreviation 
%2. take corresponding gene
%3. find the reaction in model_sdsu
%4. updates the grRules for the reaction with new gene 


load('/Users/nhinguyen/Desktop/Z.mobilis/Models/SDSUmodel_march_28_2025')
model=model_sdsu;
%step 1: read excel file
[~,geneData,~]= xlsread('SDSU_ID_gene_478.xlsx');

%put the adding new met and reaction here
% Adding problematic reactions manually
%1. Add problematicreactions using the addReaction2 function. You can find
%the structure of adding reaction in the takeRcnstoModel2
%Adding the metabolites
model_sdsu=addMetabolite(model_sdsu,'rnaasn_c','TRNA Asn C10H17O10PR2','C10H17O10PR2','','','','',0,0);%Do not modify the last 0
model_sdsu=addMetabolite(model_sdsu,'glntrna_c', 'Gln-tRNA Gln C15H23N5O11PR2', 'C15H23N5O11PR2', '', '', '', '', 0, 0);%Do the same with the other metabolite (asntrna_c)

%Adding the reactions for ASNTRS
rxnName={'ASNTRS','Asparaginyl-tRNA synthetase'};
metaboliteList={'asn__L_c','atp_c','trnaasn_c','amp_c','asntrna_c','ppi_c'};
stoichCoeffList=[-1 -1 -1 1 1 1]
revFlag=false;
lowerBound=0;
upperBound=1000;
objCoeff=0; %no modify
subSystem='tRNA Charging';
grRule=''; % no modify
geneNameList='';% no modify
systNameList='';% no modify
checkDuplicate=true;% no modify
confScores={1};% no modify
rxnECNumbers=''; % no modify
rxnKEGGID='';% no modify
rxnReferences='';% no modify
rxnNotes='';% no modify
model_sdsu.csense(length(model_sdsu.mets),1)='E';
model_sdsu=addReaction2(model_sdsu,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate,confScores);


%ZAdd grRules from both problemtaic reactions to the model
%MO0782 or ZMO0784 
index=findRxnIDs(model_sdsu,'ASNTRS');
model_sdsu.grRules{index}= 'ZMO0782 or ZMO0784' ;
model_sdsu.grRules{findRxnIDs(model_sdsu,'ASNTRS')}= 'ZMO0782 or ZMO0784' ;

%Adding the reactions for GLNTRS
rxnName={'GLNTRS','Glutaminyl-tRNA synthetase'};
metaboliteList={'atp_c','gln__L_c','trnagln_c','amp_c','glntrna_c','ppi_c'};
stoichCoeffList=[-1 -1 -1 1 1 1]
revFlag=false;
lowerBound=0;
upperBound=1000;
objCoeff=0; %no modify
subSystem='tRNA Charging';
grRule=''; % no modify
geneNameList='';% no modify
systNameList='';% no modify
checkDuplicate=true;% no modify
confScores={1};% no modify
rxnECNumbers=''; % no modify
rxnKEGGID='';% no modify
rxnReferences='';% no modify
rxnNotes='';% no modify
model_sdsu.csense(length(model_sdsu.mets),1)='E';
model_sdsu=addReaction2(model_sdsu,rxnName,metaboliteList,stoichCoeffList,revFlag,lowerBound,upperBound,objCoeff,subSystem,grRule,geneNameList,systNameList,checkDuplicate,confScores);

%ZAdd grRules from both problemtaic reactions to the model
%ZMO0782 or ZMO0783 or ZMO0784
index=findRxnIDs(model_sdsu,'GLNTRS');
model_sdsu.grRules{index}= 'ZMO0782 or ZMO0783 or ZMO0784' ;
model_sdsu.grRules{findRxnIDs(model_sdsu,'GLNTRS')}= 'ZMO0782 or ZMO0783 or ZMO0784' ;


%step 2: grab the columns
rxns = geneData(:,1) % grab abbreviatiom
new_grRules =geneData(:,2) % grab Gene

%track reaction that cant found
problematicRxns={};

%loop through each reaction
for i=1:length(rxns)
    rxn=rxns{i};  %abbreviation from excel
    new_grRule = new_grRules{i}; %Gene frome excel
    
    % Find the index of the reaction in the model
    rxn_index = findRxnIDs(model_sdsu, rxn);
    
    
    if rxn_index==0 
        problematicRxns=vertcat(problematicRxns,rxn);
        % save the reaction ID that is not found in model
    else
        model_sdsu.grRules{rxn_index} = new_grRule;
        % update the grRule in to model_sdsu
    end
end 
save('/Users/nhinguyen/Desktop/Z.mobilis/Models/SDSUmodel_april_14_2025','model_sdsu')

%% Make the model grow with minimal media and under anaerobic conditions 
% 1. Load the model SDSUmodel_april_14_2025
% 2. Try function SetMedia.m. Call the function with the correct inputs

data = load('SDSUmodel_april_14_2025.mat');
ZM_model = data.model_sdsu %model is a problematic name in MATLAB (because of Simulink), so I give it a new name= ZM_model

EX_DM_rxns=ZM_model.rxns(find(findExcRxns(ZM_model)));
for i=1:length(EX_DM_rxns)
    rxn=EX_DM_rxns{i};
    if contains(rxn,'EX_')
        ZM_model.lb(findRxnIDs(ZM_model,rxn))=0;
    end
end

% Gaby will send me exchangeSingleModel function
% Use [a,b]=exchangeSingleModel function to check what exhcnage reactions are
% opend and if the SetMedia function works properly. Focus just in the b
% output
% If when trying the function exchangeSi..., I have an error, open the lb
% of oxygen to -1000

% Gaby will send the media formulation of zymomonas used experimentally. 
% Task: get the basic molecules of these media in their exchange reactions
% and annotate them in the excel file with the media with lb=-1000 ub=1000
% For example:
% Compound in medium: NH4NO3
%   EX_nh4_e
%   EX_no3_e
% You might use bigg to get oriented. Hint: Go from bigger molecules to
% smaller ones.
%Hint: You will add EX_h2_e and EX_h_e and EX_h2o_e and EX_co2_e and
%EX_o2_e
% Look to the environmental medium that is used by the model iML1515. For
% this you have to use the exchangeSingleModel function

% GOAL :to test if a metabolic model of Zymomonas mobilis 
% can grow using a special nutrient mixture (ZMM medium) without oxygen (anaerobic conditions).
% 1. Load the model 





% 2. Load medium and apply it to the model ( media ingredient)
%[~, ~, media] = xlsread('/Users/nhinguyen/Desktop/Z.mobilis/ZM_medium_anaerobic.xlsx'); % read excel file that include list of nutrient grow in medium  
Newmodel = SetMedia('/Users/nhinguyen/Desktop/Z.mobilis/ZM_medium_anaerobic.xlsx', ZM_model); % function setmedia tells the model which compounds are available by opening the correct exchange reactions.

% close all reactions
EX_DM_rxns = Newmodel.rxns(find(findExcRxns(Newmodel)));
% findExcRxns(ZM_model) = finds all exchange reactions in this model that
% start withg EX_ and find() will give index position
% ZM_model.rxns(...) gets the reaction IDs (names like 'EX_o2_e',
% 'EX_nh4_e') from the index position
%stores the list of exchange reaction names in EX_DM_rxns
for i = 1:length(EX_DM_rxns)  % loop thru each exchange 
    rxn = EX_DM_rxns{i}; % get the name of current rxns
    if contains(rxn, 'EX_') % check to seee start with EX_
        Newmodel.lb(findRxnIDs(Newmodel, rxn)) = 0;
        % Finds the index of this reaction in the model.
        %Sets its lower bound (lb) to 0 = the reaction can no longer take in anything from the environment (no uptake allowed).
    end
end

% then make sure the oxygen exchange reaction is closed (no oxygen uptake)
% this is to ensure anaerobic condition
%EX_o2_e = 0

%close oxy
Newmodel.lb(findRxnIDs(Newmodel, 'EX_o2_e')) = 0;
% find index number of the oxygen exc in the rxn then set it = 0


% now get the model grow anaerobically :
%make sure the exchange reactions for all media components are open from the ZMM_1.28.25.rtf file


% open 
% List of required exchange reactions from ZMM
required_exchanges = {
    'EX_glc__D_e'      % glucose
    'EX_nh4_e'         % ammonium (from (NH4)2SO4)
    'EX_so4_e'         % sulfate (from (NH4)2SO4, MgSO4, FeSO4)
    'EX_pi_e'          % phosphate (from KH2PO4, K2HPO4)
    'EX_mg2_e'         % magnesium (MgSO4)
    'EX_fe2_e'         % ferrous iron (FeSO4)
    'EX_ca2_e'         % calcium (CaCl2)
    'EX_pnto__R_e'     % pantothenate
    'EX_na1_e'         % sodium (NaCl, NaMoO4)
    'EX_mobd_e'        % molybdate (NaMoO4)
    'EX_cl_e'          % chloride (NaCl, CaCl2)
    'EX_h2o_e'         % water
    'EX_h_e'           % proton
    'EX_co2_e'         % carbon dioxide
    
};
% put in excel / to run setmedia

% from ZMM zymo minimal medium 
% each of this must be open
% Open the correct exchange reactions manually (in case SetMedia missed any)

% loop thru each exc rxn
% find the index in the model
% set lower = -1000 and upper =1000
% if not foun then print warning
for i = 1:length(required_exchanges)
    rxnID = findRxnIDs(Newmodel, required_exchanges{i});
    if rxnID > 0
        Newmodel.lb(rxnID) = -1000;
        Newmodel.ub(rxnID) = 1000;
    else
        fprintf('Reaction %s not found in the model.\n', required_exchanges{i});
        % %s = insert string , \n = new line , will display reaction
        % EX_djmd_D_e not found in the model
    end
end

% Fix csense field if it's incorrectly formatted
% the csense field is missing or incorrectly formatted, the function optimizeCbModel crash with errors.
% csense = constraint sense tells what kind of constraint each row of the
% matrix represent

Newmodel.csense = repmat('E', size(Newmodel.S, 1), 1);  % char vector only

% Set biomass reaction as obj 
Newmodel.c(:) = 0;
biomass_rxn_id = findRxnIDs(Newmodel, 'BIOMASS_Ec_iML1515_WT_75p37M');
if biomass_rxn_id == 0
    error('Biomass reaction not found');
else
    Newmodel.c(biomass_rxn_id) = 1;
end

solution = optimizeCbModel(Newmodel);
fprintf('Biomass flux under anaerobic minimal media: %f\n', solution.f);

save('Newmodel_after_media.mat_04_29', 'Newmodel')









