context("initial conditions for logicODEs")
library(CNORode)

load_test_case = function(){
	
	cnodata = CNOlist(system.file("testdata/ToyModelPB_csv",package="CNORode"))
	pkn_dataGen= readSIF(system.file("testdata/ToyModelPB_TRUE_sif",package="CNORode"))
	model = preprocessing(cnodata,pkn_dataGen,expansion = F, compression = F)
	
	ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
									   LB_tau = 0, UB_n = 5, UB_k = 5, UB_tau = 10, default_n = 3,
									   default_k = 0.5, default_tau = 0.01, opt_n = TRUE, opt_k = TRUE,
									   opt_tau = TRUE, random = FALSE)
	
	
	res = list()
	res$cnodata = cnodata
	res$model = model
	res$ode_parameters = ode_parameters
	return(res)
}

test_that("create initial condition matrix", {
	
	test_case = load_test_case()
	model = test_case$model
	cnodata = test_case$cnodata
	ode_parameters = test_case$ode_parameters
	
	
	initValue1 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.2)
	expect_true(is.matrix(initValue1))
	
	initValue2 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.2,useMeasurements = F)
	expect_true(all(initValue2 == 0.2))
	
	expect_error(createLBodeInitialConditions(model = model,data = cnodata, initialValue = data.frame(0.0)))
	
	expect_error(createLBodeInitialConditions(model = model,data = cnodata, initialValue = list(0.0)))

})


test_that("simulate with initial condition matrix", {
	
	test_case = load_test_case()
	model = test_case$model
	cnodata = test_case$cnodata
	ode_parameters = test_case$ode_parameters
	
	initValue1 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.0)
	S1 = getLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue1)
	S2 = getLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters)
	expect_true(all(S1[[1]] == S2[[1]],na.rm = T))
	
	initValue1 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.0)
	S1 = getLBodeDataSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue1)
	S2 = getLBodeDataSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters)
	expect_true(all(S1[[1]] == S2[[1]],na.rm = T))
	
	initValue2 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.5)
	S1 = getLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue2)
	S2 = getLBodeModelSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters)
	expect_false(all(S1[[1]] == S2[[1]],na.rm = T))
	
	initValue2 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.5)
	S1 = getLBodeDataSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters,initialValueMatrix = initValue2)
	S2 = getLBodeDataSim(cnolist = cnodata,model = model,ode_parameters = ode_parameters)
	expect_true(all(S1[[1]] == S2[[1]],na.rm = T))
	expect_false(all(S1[[2]] == S2[[2]],na.rm = T))
	
})


test_that("call objective function", {
	
	test_case = load_test_case()
	model = test_case$model
	cnodata = test_case$cnodata
	ode_parameters = test_case$ode_parameters
	
	cnolist = compatCNOlist(cnodata)
	indices <- indexFinder(cnolist,model,verbose=FALSE)
	
	
	initValue1 = createLBodeInitialConditions(model = model,data = cnodata, initialValue = 0.0)
	
	objFun = getLBodeContObjFunction(cnolist=cnolist,
							model=model,
							ode_parameters=ode_parameters,
							indices=indices,
							time=1,
							verbose=F,
							transfer_function=4,
							reltol=1e-6,
							atol=1e-6,
							maxStepSize=2,
							maxNumSteps=1e5,
							maxErrTestsFails=50,
							nan_fac=0,
							lambda_tau=0,
							lambda_k=0,
							bootstrap=F,
							SSpenalty_fac=0,
							SScontrolPenalty_fac=0,
							boot_seed=0,
							initialValueMatrix = initValue1)
	
	expect_true( (objFun(x = ode_parameters$parValues) > 134) && (objFun(x = ode_parameters$parValues) < 135))
})

