devtools::load_all()


## Simple simulation -----------------------------------------------------------
# cnodata = CNOlist("data/ToyModel_PB/ToyModelPB.csv")
# 
# pkn_dataGen= readSIF("data/ToyModel_PB/ToyModelPB_TRUE.sif")
# plotModel(pkn_dataGen,cnodata)
# 
# 
# model = preprocessing(cnodata,pkn_dataGen,expansion = F, compression = F)
# 
# 
# ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
# 								   LB_tau = 0, UB_n = 5, UB_k = 5, UB_tau = 10, default_n = 3,
# 								   default_k = 0.5, default_tau = 0.01, opt_n = TRUE, opt_k = TRUE,
# 								   opt_tau = TRUE, random = TRUE)
# 
# 
# init_param = c("map3k7_n_nik"=1,
# 			   "map3k7_k_nik"=0.0905,
# 			   "tau_nik"=7.1398,
# 			   "map3k7_n_mkk4"=1.09,
# 			   "map3k7_k_mkk4"=0.95,
# 			   "map3k1_n_mkk4"=3.1526,
# 			   "map3k1_k_mkk4"=0.4721,
# 			   "tau_mkk4"=5.4167,
# 			   "traf2_n_ask1"=3.4873,
# 			   "traf2_k_ask1"=0.15,
# 			   "tau_ask1"=5.2731,
# 			   "traf2_n_map3k7"=1,
# 			   "traf2_k_map3k7"=0.0918,
# 			   "tau_map3k7"=0.1,
# 			   "ask1_n_mkk7"=3.5966,
# 			   "ask1_k_mkk7"=0.7772,
# 			   "map3k1_n_mkk7"=3.9788,
# 			   "map3k1_k_mkk7"=0.4026,
# 			   "tau_mkk7"=10,
# 			   "tnfa_n_tnfr"=1.0188,
# 			   "tnfa_k_tnfr"=0.1173,
# 			   "tau_tnfr"=9.1369,
# 			   "egf_n_egfr"=4.892,
# 			   "egf_k_egfr"=0.0905,
# 			   "tau_egfr"=1.0937,
# 			   "erk_n_ph"=4.611,
# 			   "erk_k_ph"=0.1027,
# 			   "tau_ph"=0.106,
# 			   "nfkb_n_ex"=4.9909,
# 			   "nfkb_k_ex"=0.4268,
# 			   "tau_ex"=0.6146,
# 			   "raf1_n_mek"=1.4875,
# 			   "raf1_k_mek"=0.5195,
# 			   "tau_mek"=9.5228,
# 			   "sos_n_ras"=1.454,
# 			   "sos_k_ras"=0.7119,
# 			   "tau_ras"=2.8385,
# 			   "tnfr_n_traf2"=4.1622,
# 			   "tnfr_k_traf2"=0.5476,
# 			   "tau_traf2"=9.9362,
# 			   "nik_n_ikk"=4.997,
# 			   "nik_k_ikk"=0.1029,
# 			   "tau_ikk"=5.7463,
# 			   "pi3k_n_akt"=1.0726,
# 			   "pi3k_k_akt"=0.9486,
# 			   "tau_akt"=9.5917,
# 			   "egfr_n_pi3k"=4.4151,
# 			   "egfr_k_pi3k"=0.0944,
# 			   "tau_pi3k"=7.7292,
# 			   "ex_n_ikb"=4.9991,
# 			   "ex_k_ikb"=0.6852,
# 			   "ikk_n_ikb"=1.5451,
# 			   "ikk_k_ikb"=0.95,
# 			   "tau_ikb"=0.6778,
# 			   "ikb_n_nfkb"=4.9997,
# 			   "ikb_k_nfkb"=0.2329,
# 			   "tau_nfkb"=0.9984,
# 			   "jnk_n_cjun"=1.0759,
# 			   "jnk_k_cjun"=0.9496,
# 			   "tau_cjun"=9.6513,
# 			   "mkk7_n_jnk"=4.8747,
# 			   "mkk7_k_jnk"=0.0913,
# 			   "tau_jnk"=1.3215,
# 			   "ras_n_map3k1"=1.1997,
# 			   "ras_k_map3k1"=0.8298,
# 			   "tau_map3k1"=9.995,
# 			   "mek_n_erk"=2.5354,
# 			   "mek_k_erk"=0.7099,
# 			   "tau_erk"=10,
# 			   "ras_n_raf1"=1,
# 			   "ras_k_raf1"=0.2725,
# 			   "tau_raf1"=0.8552,
# 			   "egfr_n_sos"=1.0012,
# 			   "egfr_k_sos"=0.9375,
# 			   "ph_n_sos"=1.4651,
# 			   "ph_k_sos"=0.72,
# 			   "tau_sos"=9.3196,
# 			   "mkk4_n_p38"=1.0032,
# 			   "mkk4_k_p38"=0.9329,
# 			   "tau_p38"=0.1,
# 			   "akt_n_gsk3"=1.0608,
# 			   "akt_k_gsk3"=0.09,
# 			   "tau_gsk3"=0.2602,
# 			   "cjun_n_ap1"=4.7985,
# 			   "cjun_k_ap1"=0.2653,
# 			   "tau_ap1"=0.2045)
# 
# parMatchID = match(ode_parameters$parNames, names(init_param))
# # match(names(init_param),ode_parameters$parNames)
# if(any(is.na(parMatchID))) stop(" some parameters were not matched in the initial parameter set")
# if(!all(ode_parameters$parNames == names(init_param)[parMatchID])) stop("some parametes were not found")
# ode_parameters$parValues = init_param[parMatchID]
# simulatedData=plotLBodeFitness(cnodata, model, transfer_function=3, ode_parameters=ode_parameters, reltol = 1e-6, atol = 1e-3, maxStepSize = 0.001)
# 

## EDGE KO cases ---------------------------------------------------------------
print("****************   edge ko   *******************************")
cnodata = CNOlist("data/ToyModel_PB/ToyModelPB.csv")

pkn_dataGen= readSIF("data/ToyModel_PB/ToyModelPB_TRUE.sif")
# plotModel(pkn_dataGen,cnodata)


model = preprocessing(cnodata,pkn_dataGen,expansion = F, compression = F)


model_KD = model

counter = which(model_KD$reacID =="!akt=gsk3")
model_KD$reacID <- model_KD$reacID[-counter]
model_KD$interMat <- model_KD$interMat[,-counter]
model_KD$notMat <- model_KD$notMat[,-counter]

# plotModel(model_KD,cnodata)


ode_parameters=createLBodeContPars(model_KD, LB_n = 1, LB_k = 0,
								   LB_tau = 0, UB_n = 5, UB_k = 5, UB_tau = 10, default_n = 3,
								   default_k = 0.5, default_tau = 0.01, opt_n = TRUE, opt_k = TRUE,
								   opt_tau = TRUE, random = TRUE)


init_param = c("map3k7_n_nik"=1,
			   "map3k7_k_nik"=0.0905,
			   "tau_nik"=7.1398,
			   "map3k7_n_mkk4"=1.09,
			   "map3k7_k_mkk4"=0.95,
			   "map3k1_n_mkk4"=3.1526,
			   "map3k1_k_mkk4"=0.4721,
			   "tau_mkk4"=5.4167,
			   "traf2_n_ask1"=3.4873,
			   "traf2_k_ask1"=0.15,
			   "tau_ask1"=5.2731,
			   "traf2_n_map3k7"=1,
			   "traf2_k_map3k7"=0.0918,
			   "tau_map3k7"=0.1,
			   "ask1_n_mkk7"=3.5966,
			   "ask1_k_mkk7"=0.7772,
			   "map3k1_n_mkk7"=3.9788,
			   "map3k1_k_mkk7"=0.4026,
			   "tau_mkk7"=10,
			   "tnfa_n_tnfr"=1.0188,
			   "tnfa_k_tnfr"=0.1173,
			   "tau_tnfr"=9.1369,
			   "egf_n_egfr"=4.892,
			   "egf_k_egfr"=0.0905,
			   "tau_egfr"=1.0937,
			   "erk_n_ph"=4.611,
			   "erk_k_ph"=0.1027,
			   "tau_ph"=0.106,
			   "nfkb_n_ex"=4.9909,
			   "nfkb_k_ex"=0.4268,
			   "tau_ex"=0.6146,
			   "raf1_n_mek"=1.4875,
			   "raf1_k_mek"=0.5195,
			   "tau_mek"=9.5228,
			   "sos_n_ras"=1.454,
			   "sos_k_ras"=0.7119,
			   "tau_ras"=2.8385,
			   "tnfr_n_traf2"=4.1622,
			   "tnfr_k_traf2"=0.5476,
			   "tau_traf2"=9.9362,
			   "nik_n_ikk"=4.997,
			   "nik_k_ikk"=0.1029,
			   "tau_ikk"=5.7463,
			   "pi3k_n_akt"=1.0726,
			   "pi3k_k_akt"=0.9486,
			   "tau_akt"=9.5917,
			   "egfr_n_pi3k"=4.4151,
			   "egfr_k_pi3k"=0.0944,
			   "tau_pi3k"=7.7292,
			   "ex_n_ikb"=4.9991,
			   "ex_k_ikb"=0.6852,
			   "ikk_n_ikb"=1.5451,
			   "ikk_k_ikb"=0.95,
			   "tau_ikb"=0.6778,
			   "ikb_n_nfkb"=4.9997,
			   "ikb_k_nfkb"=0.2329,
			   "tau_nfkb"=0.9984,
			   "jnk_n_cjun"=1.0759,
			   "jnk_k_cjun"=0.9496,
			   "tau_cjun"=9.6513,
			   "mkk7_n_jnk"=4.8747,
			   "mkk7_k_jnk"=0.0913,
			   "tau_jnk"=1.3215,
			   "ras_n_map3k1"=1.1997,
			   "ras_k_map3k1"=0.8298,
			   "tau_map3k1"=9.995,
			   "mek_n_erk"=2.5354,
			   "mek_k_erk"=0.7099,
			   "tau_erk"=10,
			   "ras_n_raf1"=1,
			   "ras_k_raf1"=0.2725,
			   "tau_raf1"=0.8552,
			   "egfr_n_sos"=1.0012,
			   "egfr_k_sos"=0.9375,
			   "ph_n_sos"=1.4651,
			   "ph_k_sos"=0.72,
			   "tau_sos"=9.3196,
			   "mkk4_n_p38"=1.0032,
			   "mkk4_k_p38"=0.9329,
			   "tau_p38"=0.1,
			   "akt_n_gsk3"=1.0608,
			   "akt_k_gsk3"=0.09,
			   "tau_gsk3"=0.2602,
			   "cjun_n_ap1"=4.7985,
			   "cjun_k_ap1"=0.2653,
			   "tau_ap1"=0.2045)

parMatchID = match(ode_parameters$parNames, names(init_param))
# match(names(init_param),ode_parameters$parNames) 
if(any(is.na(parMatchID))) stop(" some parameters were not matched in the initial parameter set")
if(!all(ode_parameters$parNames == names(init_param)[parMatchID])) stop("some parametes were not found")
ode_parameters$parValues = init_param[parMatchID]
simulatedData=plotLBodeFitness(cnodata, model_KD, transfer_function=3, ode_parameters=ode_parameters, reltol = 1e-6, atol = 1e-3, maxStepSize = 0.001)



# KO other pi3k->AKT
# 
# model_KD = model
# counter = which(model_KD$reacID =="pi3k=akt")
# model_KD$reacID <- model_KD$reacID[-counter]
# model_KD$interMat <- model_KD$interMat[,-counter]
# model_KD$notMat <- model_KD$notMat[,-counter]
# 
# plotModel(model_KD,cnodata)
# 
# 
# ode_parameters=createLBodeContPars(model_KD, LB_n = 1, LB_k = 0,
# 								   LB_tau = 0, UB_n = 5, UB_k = 5, UB_tau = 10, default_n = 3,
# 								   default_k = 0.5, default_tau = 0.01, opt_n = TRUE, opt_k = TRUE,
# 								   opt_tau = TRUE, random = TRUE)
# 
# parMatchID = match(ode_parameters$parNames, names(init_param))
# # match(names(init_param),ode_parameters$parNames) 
# if(any(is.na(parMatchID))) stop(" some parameters were not matched in the initial parameter set")
# if(!all(ode_parameters$parNames == names(init_param)[parMatchID])) stop("some parametes were not found")
# ode_parameters$parValues = init_param[parMatchID]
# simulatedData=plotLBodeFitness(cnodata, model_KD, transfer_function=3, ode_parameters=ode_parameters, reltol = 1e-6, atol = 1e-3, maxStepSize = 0.001)
# 
# 
# 
