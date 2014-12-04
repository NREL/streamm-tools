.. _db_interface:

.. index:: query, SQL, database, OPV

*********************************************
Interface to NREL database
*********************************************

One of the projects at NREL that uses the STREAMM toolkit is the
Organic Photovoltaics (OPV) modeling effort. High-throughput density
functional theory (DFT) calculations are used to screen candidate
molecules so as to let synthetic chemists focus on only the most
promising materials for energy applications. The simulation results
from ~10,000's of calculations are stored in the NREL OPV database
The STREAMM toolkit is used not only to generate the simulation files
needed to perform the DFT calculations, but also to access and
effectively mine the massive amounts of data generated on the
Peregrine compute cluster. The python modules and scripts that enable
efficient access to the OPV database is described below.


Description of OPV data
=======================================

=============================    ======================
Column                           Type
=============================    ======================
project_id                       integer
project_name                     character varying(150)
oligomer_id                      integer
tag                              character varying(150)
oligomer_number                  integer
backbone_type_name               character varying
molecule_type                    character varying(150)       
bblock_list                      character varying[]          
isomer_type                      character varying(10)        
conformer_type                   character varying(10)        
oligomer_created_on              timestamp without time zone  
oligomer_result_id               integer                      
basis                            character varying(150)       
oligomer_result_date_added       timestamp without time zone  
oligomer_result_created_on       timestamp without time zone  
date_calculation_finished        timestamp without time zone  
oligomer_result_last_updated     timestamp without time zone  
homo                             numeric(45,20)               
lumo                             numeric(45,20)               
optical_lumo                     numeric(45,20)               
homo_contour                     numeric(45,20)               
lumo_contour                     numeric(45,20)               
gap                              numeric(45,20)               
logp                             numeric(45,20)               
mu_x                             numeric(45,20)               
mu_y                             numeric(45,20)               
mu_z                             numeric(45,20)               
mu_tot                           numeric(45,20)               
total_energy                     numeric(45,20)               
mu2_x                            numeric                      
mu2_y                            numeric                      
mu2_z                            numeric                      
mu2_tot                          numeric                      
dipole_difference                numeric                      
date_cation_calc_finished        timestamp without time zone  
vert_cat_en                      numeric(45,20)               
neutral_reorganization_energy    numeric(45,20)               
cation_reorganization_energy     numeric(45,20)               
vertical_ionization_energy       numeric(45,20)               
adiabatic_ionization_energy      numeric(45,20)               
relaxed_neut_en                  numeric(45,20)               
relaxed_cat_en                   numeric(45,20)               
mol                              text
=============================    ======================


vw_structure

=============================    ======================
Column                           Type
=============================    ======================
project_id                       integer                      
project_name                     character varying(150)       
structure_id                     integer                      
tag                              character varying(150)       
backbone_type_name               character varying            
max_oligomer_number              integer                      
isomer_type                      character varying(10)        
conformer_type                   character varying(10)        
structure_created_on             timestamp without time zone  
structure_result_id              integer                      
basis                            character varying(150)       
structure_result_date_added      timestamp without time zone  
structure_result_created_on      timestamp without time zone  
delta_homo                       numeric(45,20)               
delta_lumo                       numeric(45,20)               
delta_optical_lumo               numeric(45,20)               
homo_extrapolated                numeric(45,20)               
lumo_extrapolated                numeric(45,20)               
gap_extrapolated                 numeric(45,20)               
optical_lumo_extrapolated        numeric(45,20)
=============================    ======================


Desciption of opvSQL python module
======================================================

.. automodule:: opvSQL
   :members:
   :show-inheritance:


Examples of using STREAMM classes with database
======================================================

simulationGaussian example
