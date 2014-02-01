import FWCore.ParameterSet.Config as cms

branchesToSkim = cms.vstring()
branchesToSkim.extend((
 'bc_chx',
 'bc_hitdetid',
 'bc_hybrid_n',
 'bc_islbar_n',
 'bc_islend_n',
 'bc_nhits',
 'bc_s4',
 'bc_sieip',
 'bc_sipip',
 'bc_type',
 'bs_sigmaZ0Error',
 'bs_x0Error',
 'bs_y0Error',
 'bs_z0Error',
 'conv_MVALikelihood',
 'conv_ch1ch2',
 'conv_detatrksatecal',
 'conv_distofminapproach',
 'conv_dphitrksatecal',
 'conv_dphitrksatvtx',
 'conv_dxy',
 'conv_dz',
 'conv_eleind',
 'conv_lxy',
 'conv_lz',
 'conv_nHitsBeforeVtx',
 'conv_nHitsMax',
 'conv_nSharedHits',
 'conv_p4',
 'conv_pair_momentum',
 'conv_paircotthetasep',
 'conv_pairinvmass',
 'conv_quality',
 'conv_tk1_d0',
 'conv_tk1_dz',
 'conv_tk1_dzerr',
 'conv_tk1_etaerr',
 'conv_tk1_lambdaerr',
 'conv_tk1_nh',
 'conv_tk1_phierr',
 'conv_tk1_pin',
 'conv_tk1_pout',
 'conv_tk1_pterr',
 'conv_tk1_thetaerr',
 'conv_tk2_d0',
 'conv_tk2_dz',
 'conv_tk2_dzerr',
 'conv_tk2_etaerr',
 'conv_tk2_lambdaerr',
 'conv_tk2_nh',
 'conv_tk2_phierr',
 'conv_tk2_pin',
 'conv_tk2_pout',
 'conv_tk2_pterr',
 'conv_tk2_thetaerr',
 'conv_type',
 'conv_vtxProb',
 'conv_vtx_xErr',
 'conv_vtx_yErr',
 'conv_vtx_zErr',
 'conv_zofprimvtxfromtrks',
 'ct_n',
 'ct_p4',
 'ct_emEnergy',
 'ct_emL1',
 'ct_hadEnergy',
 'ct_hadL1',
 'ct_outerEnergy',
 'ct_size',
 'ecalhit_flag',
 'ecalhit_ieta',
 'ecalhit_iphi',
 'ecalhit_ix',
 'ecalhit_iy',
 'ecalhit_time',
 'ecalhit_type',
 'ecalhit_zside',
 'el_std_bchits',
 'el_std_calib_energy',
 'el_std_calib_energyerr',
 'el_std_catbased',
 'el_std_ch_gsf',
 'el_std_ch_scpix',
 'el_std_chi2',
 'el_std_class',
 'el_std_conv_vtxProb',
 'el_std_corr_energy',
 'el_std_corr_energyerr',
 'el_std_crack',
 'el_std_d0',
 'el_std_detaeleout',
 'el_std_detaout',
 'el_std_dphiout',
 'el_std_e1x5',
 'el_std_e2x5',
 'el_std_e5x5',
 'el_std_eleopout',
 'el_std_epf',
 'el_std_eseedopout',
 'el_std_eseedpf',
 'el_std_eseffsixix',
 'el_std_eseffsiyiy',
 'el_std_gsfchi2',
 'el_std_hcalbciso03',
 'el_std_hcalbciso04',
 'el_std_hcalsolidiso03',
 'el_std_hcalsolidiso04',
 'el_std_hoebc',
 'el_std_hoebcd1',
 'el_std_hoebcd2',
 'el_std_hoed1',
 'el_std_hoed2',
 'el_std_hp_1pxb',
 'el_std_hp_1pxf',
 'el_std_hp_expout',
 'el_std_ip3d',
 'el_std_ip3d_err',
 'el_std_ip3d_sig',
 'el_std_ip_ctf',
 'el_std_kfchi2',
 'el_std_kfhits',
 'el_std_losthits',
 'el_std_momcalo',
 'el_std_momout',
 'el_std_momvtx',
 'el_std_momvtxconst',
 'el_std_must',
 'el_std_mustnc',
 'el_std_nambtk',
 'el_std_nbrem',
 'el_std_nbrempf',
 'el_std_passcutpresel',
 'el_std_passmvapresel',
 'el_std_poscalo',
 'el_std_pout',
 'el_std_psenergy',
 'el_std_psenergypf',
 'el_std_psly1',
 'el_std_psly2',
 'el_std_psnstriply1',
 'el_std_psnstriply2',
 'el_std_r9',
 'el_std_sc_time',
 'el_std_sieiesc',
 'el_std_schits',
 'el_std_sipip',
 'el_std_tkiso04',
 'el_std_validhits',
 'el_std_z0',
 'filter_names_HLT1',
 'filter_pass',
 'genjet_algo1_aux',
 'genjet_algo1_em',
 'genjet_algo1_had',
 'genjet_algo1_inv',
 'genjet_algo2_aux',
 'genjet_algo2_em',
 'genjet_algo2_had',
 'genjet_algo2_inv',
 'genjet_algo3_aux',
 'genjet_algo3_em',
 'genjet_algo3_had',
 'genjet_algo3_inv',
 'gp_vtx',
 'gsf_tk_charge',
 'gsf_tk_chi2',
 'gsf_tk_d0',
 'gsf_tk_d0err',
 'gsf_tk_dof',
 'gsf_tk_dz',
 'gsf_tk_dzerr',
 'gsf_tk_etaerr',
 'gsf_tk_fbrem',
 'gsf_tk_hp_expin',
 'gsf_tk_hp_expout',
 'gsf_tk_hp_nlost',
 'gsf_tk_hp_nvalid',
 'gsf_tk_hp_nvalidpix',
 'gsf_tk_n',
 'gsf_tk_nhits',
 'gsf_tk_nlosthit',
 'gsf_tk_p4',
 'gsf_tk_phierr',
 'gsf_tk_pin',
 'gsf_tk_pinmode',
 'gsf_tk_pout',
 'gsf_tk_poutmode',
 'gsf_tk_pterr',
 'gsf_tk_qoverpinerr',
 'gsf_tk_qoverpouterr',
 'gsf_tk_shared',
 'gsf_tk_tkind',
 'gsf_tk_tpind',
 'gsf_tk_vtx_pos',
 'gv_nTkHi',
 'gv_nTkLo',
 'gv_p3',
 'gv_sumPtHi',
 'gv_sumPtLo',
 'hlt_candpath2',
 'jet_algoPF1_R',
 'jet_algoPF1_Rchg',
 'jet_algoPF1_Rchg_QC',
 'jet_algoPF1_Rneutral',
 'jet_algoPF1_axis1',
 'jet_algoPF1_axis1_QC',
 'jet_algoPF1_axis2',
 'jet_algoPF1_calotwind',
 'jet_algoPF1_cutbased_mva',
 'jet_algoPF1_cutbased_mva_ext',
 'jet_algoPF1_frac06',
 'jet_algoPF1_frac07',
 'jet_algoPF1_nCharged_ptCut_QC',
 'jet_algoPF1_pull',
 'jet_algoPF1_pull_QC',
 'jet_algoPF1_pull_dphi',
 'jet_algoPF1_pull_dy',
 'jet_algoPF1_rmsCand',
 'jet_algoPF1_rmsCand_QC',
 'jet_algoPF1_tana',
 'jet_algoPF1_tana_QC',
 'jet_algoPF3_R',
 'jet_algoPF3_Rchg',
 'jet_algoPF3_Rchg_QC',
 'jet_algoPF3_Rneutral',
 'jet_algoPF3_axis1',
 'jet_algoPF3_axis1_QC',
 'jet_algoPF3_axis2',
 'jet_algoPF3_axis2_QC',
 'jet_algoPF3_calotwind',
 'jet_algoPF3_cutbased_mva',
 'jet_algoPF3_cutbased_mva_ext',
 'jet_algoPF3_frac06',
 'jet_algoPF3_frac07',
 'jet_algoPF3_nCharged_QC',
 'jet_algoPF3_nCharged_ptCut_QC',
 'jet_algoPF3_nNeutrals_ptCut',
 'jet_algoPF3_ptD_QC',
 'jet_algoPF3_pull',
 'jet_algoPF3_pull_QC',
 'jet_algoPF3_pull_dphi',
 'jet_algoPF3_pull_dy',
 'jet_algoPF3_rmsCand',
 'jet_algoPF3_rmsCand_QC',
 'jet_algoPF3_tana',
 'jet_algoPF3_tana_QC',
 'l1_labels',
 'l1bits_phy',
 'l1bits_tec',
 'l1cenjet_et',
 'l1cenjet_eta',
 'l1cenjet_n',
 'l1cenjet_phi',
 'l1forjet_et',
 'l1forjet_eta',
 'l1forjet_n',
 'l1forjet_phi',
 'l1met_et',
 'l1met_phi',
 'l1mu_et',
 'l1mu_eta',
 'l1mu_n',
 'l1mu_phi',
 'l1taujet_et',
 'l1taujet_eta',
 'l1taujet_n',
 'l1taujet_phi',
 'lpt_drmatch',
 'lpt_duplicate',
 'lpt_el_n',
 'lpt_emu_n',
 'lpt_ind',
 'lpt_indgen',
 'lpt_mu_n',
 'lpt_n',
 'lpt_p4',
 'lpt_pdgid',
 'lpt_pho_n',
 'met_met',
 'met_met_crossed',
 'met_met_jet',
 'met_met_mip',
 'met_met_nocalo',
 'met_met_s9',
 'met_phi',
 'met_phi_crossed',
 'met_phi_jet',
 'met_phi_mip',
 'met_phi_nocalo',
 'met_phi_s9',
 'met_phi_tcmet',
 'mu_glo_d0',
 'mu_glo_d0err',
 'mu_glo_dz',
 'mu_glo_dzerr',
 'mu_glo_em',
 'mu_glo_emS9',
 'mu_glo_had',
 'mu_glo_hadS9',
 'mu_glo_ho',
 'mu_glo_hoS9',
 'mu_glo_innerhits',
 'mu_glo_losthits',
 'mu_glo_momvtx',
 'mu_glo_posecal',
 'mu_glo_poshcal',
 'mu_glo_posvtx',
 'mu_glo_staind',
 'mu_glo_tkpterr',
 'mu_rhoCorr',
 'pfcand_ecalenergy',
 'pfcand_hcalenergy',
 'pfcand_ispu',
 'pfcand_momerr',
 'pfcand_mva_e_mu',
 'pfcand_mva_e_pi',
 'pfcand_mva_gamma_nh',
 'pfcand_mva_nothing_gamma',
 'pfcand_mva_nothing_nh',
 'pfcand_mva_pi_mu',
 'pfcand_overlappho',
 'pfcand_ps1energy',
 'pfcand_ps2energy',
 'pfcand_rawecalenergy',
 'pfcand_rawhcalenergy',
 'pfcand_vz',
 'pho_IsConvOutIn',
 'pho_PfEleVeto',
 'pho_barrel',
 'pho_bchits',
 'pho_cep',
 'pho_cep_global',
 'pho_conv_MVALikelihood',
 'pho_conv_ch1ch2',
 'pho_conv_chi2',
 'pho_conv_chi2_probability',
 'pho_conv_detatrksatecal',
 'pho_conv_distofminapproach',
 'pho_conv_dphitrksatecal',
 'pho_conv_dphitrksatvtx',
 'pho_conv_ntracks',
 'pho_conv_pair_momentum',
 'pho_conv_paircotthetasep',
 'pho_conv_pairinvmass',
 'pho_conv_tk1_d0',
 'pho_conv_tk1_dz',
 'pho_conv_tk1_dzerr',
 'pho_conv_tk1_nh',
 'pho_conv_tk1_pin',
 'pho_conv_tk1_pout',
 'pho_conv_tk2_d0',
 'pho_conv_tk2_dz',
 'pho_conv_tk2_dzerr',
 'pho_conv_tk2_nh',
 'pho_conv_tk2_pin',
 'pho_conv_tk2_pout',
 'pho_conv_validvtx',
 'pho_conv_vertexcorrected_p4',
 'pho_conv_zofprimvtxfromtrks',
 'pho_e1x3',
 'pho_h1oe',
 'pho_h1oe_bc',
 'pho_h2oe',
 'pho_h2oe_bc',
 'pho_hasConvPf',
 'pho_hasSLConvPf',
 'pho_hcal1sumetconedr03',
 'pho_hcal1sumetconedr04',
 'pho_hcal2sumetconedr03',
 'pho_hcal2sumetconedr04',
 'pho_hcalbc1sumetconedr03',
 'pho_hcalbc1sumetconedr04',
 'pho_hcalbc2sumetconedr03',
 'pho_hcalbc2sumetconedr04',
 'pho_hcalbcsumetconedr03',
 'pho_hcalbcsumetconedr04',
 'pho_isEBEEGap',
 'pho_isEBEtaGap',
 'pho_isEBGap',
 'pho_isEBPhiGap',
 'pho_isEEDeeGap',
 'pho_isEEGap',
 'pho_isEERingGap',
 'pho_lambdadivcov',
 'pho_lambdadivcov_global',
 'pho_lambdaratio_global',
 'pho_maxoraw',
 'pho_mustEtout',
 'pho_ntrkhollowconedr04',
 'pho_ntrksolidconedr03',
 'pho_ntrksolidconedr04',
 'pho_pfClusECorr',
 'pho_pfMatch',
 'pho_pfconvVtxZ',
 'pho_pfconvVtxZErr',
 'pho_pi0disc',
 'pho_r19',
 'pho_r1x5',
 'pho_r2x5',
 'pho_schits',
 'pho_see',
 'pho_seed_time',
 'pho_smaj',
 'pho_trksumpthollowconedr04',
 'pho_trksumptsolidconedr03',
 'pho_trksumptsolidconedr04',
 'pho_zernike20',
 'pho_zernike42',
 'rho_algo2',
 'rho_algo3',
 'sc_2xN',
 'sc_5xN',
 'sc_barrel',
 'sc_bccrackcorr',
 'sc_bclocalcorr',
 'sc_brem',
 'sc_hybrid_n',
 'sc_islbar_bcind',
 'sc_islbar_bcseedind',
 'sc_islbar_n',
 'sc_islbar_nbc',
 'sc_islbar_p4',
 'sc_islbar_raw',
 'sc_islbar_seedenergy',
 'sc_islbar_xyz',
 'sc_islend_n',
 'sc_r9',
 'sc_sieie',
 'selector_bits',
 'tk_algo',
 'tk_charge',
 'tk_chi2',
 'tk_d0',
 'tk_d0err',
 'tk_dof',
 'tk_dz',
 'tk_dzerr',
 'tk_etaerr',
 'tk_hp_expin',
 'tk_hp_expout',
 'tk_hp_nlost',
 'tk_hp_nvalid',
 'tk_hp_nvalidpix',
 'tk_nhits',
 'tk_nlosthit',
 'tk_phierr',
 'tk_qoverperr',
 'tk_tpind',
 'vtx_nobs_dxdydz',
 'vtx_nobs_n',
 'vtx_nobs_ndof',
 'vtx_nobs_ntks',
 'vtx_nobs_scalarpt',
 'vtx_nobs_tkind',
 'vtx_nobs_tkweight',
 'vtx_nobs_vectorp3',
 'vtx_nobs_x2dof',
 'vtx_nobs_xyz',
 'vtx_std_ndof',
 'vtx_std_scalarpt',
 'vtx_std_vectorp3'))

