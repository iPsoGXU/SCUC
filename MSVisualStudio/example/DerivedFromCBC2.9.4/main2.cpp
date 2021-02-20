#include "CbcConfig.h"
#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinTime.hpp"
#include "CglPreProcess.hpp"
#include "CbcSOS.hpp"
#include "CbcBranchLotsize.hpp"
#include "CbcSimpleInteger.hpp"
#include "CglProbing.hpp"
#include "CoinFinite.hpp"
#include "CbcLinked.hpp"
#include "CbcSolverHeuristics.hpp"
#include "CbcOrClpParam.hpp"
#include "CbcSolver.hpp"
#include "CglKnapsackCover.hpp"
#include "CglGMI.hpp"
#include "CglGomory.hpp"
#include "CglRedSplit.hpp"
#include "CglRedSplit2.hpp"
#include "CglClique.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglFlowCover.hpp"
#include "CglTwomir.hpp"
#include "CglLandP.hpp"
#include "CglResidualCapacity.hpp"
#include "CglZeroHalf.hpp"
#include "CbcCutGenerator.hpp"
#include "CoinMpsIO.hpp"
#include "CbcSolverExpandKnapsack.hpp"
#include "CbcMipStartIO.hpp"
#include "CbcSimpleIntegerPseudoCost.hpp"
#include "CbcCompareDefault.hpp"
#include "CbcStrategy.hpp"
#include "CbcBranchDefaultDecision.hpp"
#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>

#define IN_BRANCH_AND_BOUND (0x01000000|262144|128|1024|2048)
static void putBackOtherSolutions(CbcModel* presolvedModel, CbcModel* model,
	CglPreProcess* preProcess)
{
	int numberSolutions = presolvedModel->numberSavedSolutions();
	int numberColumns = presolvedModel->getNumCols();
	if (numberSolutions > 1) {
		model->deleteSolutions();
		double* bestSolution = CoinCopyOfArray(presolvedModel->bestSolution(), numberColumns);
		//double cutoff = presolvedModel->getCutoff();
		double objectiveValue = presolvedModel->getObjValue();
		//model->createSpaceForSavedSolutions(numberSolutions-1);
		for (int iSolution = numberSolutions - 1; iSolution >= 0; iSolution--) {
			presolvedModel->setCutoff(COIN_DBL_MAX);
			presolvedModel->solver()->setColSolution(presolvedModel->savedSolution(iSolution));
			//presolvedModel->savedSolutionObjective(iSolution));
			preProcess->postProcess(*presolvedModel->solver(), false);
			model->setBestSolution(preProcess->originalModel()->getColSolution(), model->solver()->getNumCols(),
				presolvedModel->savedSolutionObjective(iSolution));
		}
		presolvedModel->setBestObjectiveValue(objectiveValue);
		presolvedModel->solver()->setColSolution(bestSolution);
		//presolvedModel->setBestSolution(bestSolution,numberColumns,objectiveValue);
	}
}

int main()
{
	// 读取数据，实例化模型
	std::string dataFileName = "uc10.mps";
	OsiClpSolverInterface initialSolver;
	int numMpsReadErrors = initialSolver.readMps(dataFileName.c_str(), "");
	if (numMpsReadErrors != 0)
	{
		printf("%d errors reading MPS file\n", numMpsReadErrors);
		return numMpsReadErrors;
	}
	CbcModel model(initialSolver);

	// 基本参数，并设置默认值
	CbcSolverUsefulData parameterData;
	CbcMain0(model, parameterData);
	parameterData.noPrinting_ = false;
	parameterData.useSignalHandler_ = true;
	CbcOrClpParam* parameters = parameterData.parameters_;
	int numberParameters = parameterData.numberParameters_;
	double totalTime = parameterData.totalTime_; // 记录运行总时间
	bool noPrinting = parameterData.noPrinting_;
	bool useSignalHandler = parameterData.useSignalHandler_;
	int doIdiot = -1;
	int outputFormat = 2;
	int substitution = 3;
	int dualize = 3;
	int preSolve = 5;
	int doSprint = -1;
	int testOsiParameters = -1;
	int nodeStrategy = 0;
	bool useStrategy = false;
	bool strongChanged = false;
	bool pumpChanged = false;
	bool biLinearProblem = false;
	int returnMode = 1;
	int numberChanged = 0;

	// Say no resolve after cuts 
	model.setResolveAfterTakeOffCuts(false);
	model.setNumberBeforeTrust(10);
	model.setNumberStrong(5);
	OsiSolverInterface* solver = model.solver();
	OsiClpSolverInterface* clpSolver = dynamic_cast<OsiClpSolverInterface*> (solver);
	ClpSimplex* lpSolver = clpSolver->getModelPtr();
	lpSolver->setPerturbation(50);
	lpSolver->messageHandler()->setPrefix(false);
	clpSolver->messageHandler()->setLogLevel(1);
	model.messageHandler()->setLogLevel(1);
	lpSolver->setLogLevel(1);

	int slog = whichParam(CLP_PARAM_INT_SOLVERLOGLEVEL, numberParameters, parameters);
	int log = whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters, parameters);
	

	// 统计用变量
	double statistics_seconds = 0.0, statistics_obj = 0.0;
	double statistics_sys_seconds = 0.0, statistics_elapsed_seconds = 0.0;
	double statistics_continuous = 0.0, statistics_tighter = 0.0;
	double statistics_cut_time = 0.0;
	int statistics_nodes = 0, statistics_iterations = 0;
	int statistics_nrows = 0, statistics_ncols = 0;
	int statistics_nprocessedrows = 0, statistics_nprocessedcols = 0;
	std::string statistics_result;
	int* statistics_number_cuts = NULL;
	const char** statistics_name_generators = NULL;
	int statistics_number_generators = 0;


	OsiClpSolverInterface* originalSolver = dynamic_cast<OsiClpSolverInterface*> (model.solver());

	// 定义cuts备用，action代表了cut的默认使用方法
	CglGomory gomoryGen;
	// try larger limit
	gomoryGen.setLimitAtRoot(1000);
	gomoryGen.setLimit(50);
	gomoryGen.setAwayAtRoot(0.005);
	// set default action (0=off,1=on,2=root)
	int gomoryAction = 3;

	CglProbing probingGen;
	probingGen.setUsingObjective(1);
	probingGen.setMaxPass(1);
	probingGen.setMaxPassRoot(1);
	// Number of unsatisfied variables to look at
	probingGen.setMaxProbe(10);
	probingGen.setMaxProbeRoot(50);
	// How far to follow the consequences
	probingGen.setMaxLook(10);
	probingGen.setMaxLookRoot(50);
	probingGen.setMaxLookRoot(10);
	// Only look at rows with fewer than this number of elements
	probingGen.setMaxElements(200);
	probingGen.setMaxElementsRoot(300);
	probingGen.setRowCuts(3);
	// set default action (0=off,1=on,2=root)
	int probingAction = 3;

	CglKnapsackCover knapsackGen;
	//knapsackGen.switchOnExpensive();
	//knapsackGen.setMaxInKnapsack(100);
	// set default action (0=off,1=on,2=root)
	int knapsackAction = 3;

	CglRedSplit redsplitGen;
	//redsplitGen.setLimit(100);
	// set default action (0=off,1=on,2=root)
	// Off as seems to give some bad cuts
	int redsplitAction = 0;

	CglRedSplit2 redsplit2Gen;
	//redsplit2Gen.setLimit(100);
	// set default action (0=off,1=on,2=root)
	// Off
	int redsplit2Action = 0;

	CglGMI GMIGen;
	//GMIGen.setLimit(100);
	// set default action (0=off,1=on,2=root)
	int GMIAction = 0;

	CglFakeClique cliqueGen(NULL, false);
	//CglClique cliqueGen(false,true);
	cliqueGen.setStarCliqueReport(false);
	cliqueGen.setRowCliqueReport(false);
	cliqueGen.setMinViolation(0.1);
	// set default action (0=off,1=on,2=root)
	int cliqueAction = 3;

	// maxaggr,multiply,criterion(1-3)
	CglMixedIntegerRounding2 mixedGen(1, true, 1);
	// set default action (0=off,1=on,2=root)
	int mixedAction = 3;
	mixedGen.setDoPreproc(1); // safer (and better)

	CglFlowCover flowGen;
	// set default action (0=off,1=on,2=root)
	int flowAction = 3;

	CglTwomir twomirGen;
	twomirGen.setMaxElements(250);
	twomirGen.setAwayAtRoot(0.005);
	twomirGen.setAway(0.01);
	// set default action (0=off,1=on,2=root)
	int twomirAction = 3;

	CglLandP landpGen;
	landpGen.validator().setMinViolation(1.0e-4);
	// set default action (0=off,1=on,2=root)
	int landpAction = 0;

	CglResidualCapacity residualCapacityGen;
	residualCapacityGen.setDoPreproc(1); // always preprocess
	// set default action (0=off,1=on,2=root)
	int residualCapacityAction = 0;

	CglZeroHalf zerohalfGen;
	int zerohalfAction = 0;

	bool dominatedCuts = false;

	int initialPumpTune = -1;
	int iParam;
	iParam = whichParam(CBC_PARAM_INT_DIVEOPT, numberParameters, parameters);
	parameters[iParam].setIntValue(2);
	iParam = whichParam(CBC_PARAM_INT_FPUMPITS, numberParameters, parameters);
	parameters[iParam].setIntValue(30);
	iParam = whichParam(CBC_PARAM_INT_FPUMPTUNE, numberParameters, parameters);
	parameters[iParam].setIntValue(1005043);
	initialPumpTune = 1005043;
	iParam = whichParam(CLP_PARAM_INT_PROCESSTUNE, numberParameters, parameters);
	parameters[iParam].setIntValue(6);
	iParam = whichParam(CBC_PARAM_STR_DIVINGC, numberParameters, parameters);
	parameters[iParam].setCurrentOption("on");
	iParam = whichParam(CBC_PARAM_STR_RINS, numberParameters, parameters);
	parameters[iParam].setCurrentOption("on");
	iParam = whichParam(CBC_PARAM_STR_PROBINGCUTS, numberParameters, parameters);
	parameters[iParam].setCurrentOption("on");

	// 开始计时
	double time1 = CoinCpuTime();
	double time1Elapsed = CoinGetTimeOfDay();

	int doScaling = 4;
	int logLevel = parameters[slog].intValue();
	// 0 normal, 1 from ampl or MIQP etc (2 allows cuts)
	int complicatedInteger = 0;
	{
		OsiSolverInterface* solver = model.solver();
		OsiClpSolverInterface* si = dynamic_cast<OsiClpSolverInterface*>(solver);
		assert(si != NULL);
		si->getModelPtr()->scaling(doScaling);
		ClpSimplex* lpSolver = si->getModelPtr();
		statistics_nrows = si->getNumRows();
		statistics_ncols = si->getNumCols();
		statistics_nprocessedrows = si->getNumRows();
		statistics_nprocessedcols = si->getNumCols();
		if (logLevel <= 1)
			si->setHintParam(OsiDoReducePrint, true, OsiHintTry);
		si->setSpecialOptions(0x40000000);
	}

	// 开始预处理 
	OsiSolverInterface* saveSolver = NULL;
	CglPreProcess process;
	CbcModel* babModel = new CbcModel(model);
	OsiSolverInterface* solver3 = clpSolver->clone();
	babModel->assignSolver(solver3);
	OsiClpSolverInterface* clpSolver2 = dynamic_cast<OsiClpSolverInterface*> (babModel->solver());
	if (clpSolver2->messageHandler()->logLevel())
		clpSolver2->messageHandler()->setLogLevel(1);
	if (logLevel > -1)
		clpSolver2->messageHandler()->setLogLevel(logLevel);
	lpSolver = clpSolver2->getModelPtr();
	if (lpSolver->factorizationFrequency() == 200)
	{
		// User did not touch preset
		int numberRows = lpSolver->numberRows();
		const int cutoff1 = 10000;
		const int cutoff2 = 100000;
		const int base = 75;
		const int freq0 = 50;
		const int freq1 = 200;
		const int freq2 = 400;
		const int maximum = 1000;
		int frequency;
		if (numberRows < cutoff1)
			frequency = base + numberRows / freq0;
		else if (numberRows < cutoff2)
			frequency = base + cutoff1 / freq0 + (numberRows - cutoff1) / freq1;
		else
			frequency = base + cutoff1 / freq0 + (cutoff2 - cutoff1) / freq1 + (numberRows - cutoff2) / freq2;
		lpSolver->setFactorizationFrequency(CoinMin(maximum, frequency));
	}
	double time2 = CoinCpuTime();
	totalTime += time2 - time1;
	double timeLeft = babModel->getMaximumSeconds();
	int numberOriginalColumns = babModel->solver()->getNumCols();

	int preProcess = 4;
	if (preProcess)
	{
		// see whether to switch off preprocessing
		// only allow SOS and integer
		OsiObject** objects = babModel->objects();
		int numberObjects = babModel->numberObjects();
		for (int iObj = 0; iObj < numberObjects; iObj++)
		{
			CbcSOS* objSOS = dynamic_cast <CbcSOS*>(objects[iObj]);
			CbcSimpleInteger* objSimpleInteger = dynamic_cast <CbcSimpleInteger*>(objects[iObj]);
			if (!objSimpleInteger && !objSOS)
			{
				// find all integers anyway
				babModel->findIntegers(true);
				preProcess = 0;
				break;
			}
		}
	}
	int numberSOS = 0;
	int doSOS = 1;
	int* sosStart = NULL;
	char* sosType = NULL;
	int* sosIndices = NULL;
	double* sosReference = NULL;

	int integerStatus = -1; // 表示解是否为整数
	int truncateRows = -1;
	int* newPriorities = NULL;
	double* truncatedRhsLower = NULL;
	double* truncatedRhsUpper = NULL;
	int truncateColumns = COIN_INT_MAX;
	bool integersOK;
	int cutPass = -1234567;
	CbcModel* currentBranchModel = NULL;
	std::vector< std::pair< std::string, double > > mipStart;
	std::vector< std::pair< std::string, double > > mipStartBefore;
	double* solutionIn = NULL;
	int useSolution = -1;
	int* prioritiesIn = NULL;
	double* pseudoDown = NULL;
	double* pseudoUp = NULL;
	int* priorities = NULL;
	int* branchDirection = NULL;
	int* cut = NULL;
	int* sosPriority = NULL;
	int printOptions = 0;

	// 预处理的整个过程
	if (preProcess)
	{
		saveSolver = babModel->solver()->clone();
		/* Do not try and produce equality cliques and
		   do up to 10 passes */
		OsiSolverInterface* solver2;
		int tunePreProcess = 0;  // 控制是否使用更强的probing
		{
			// Tell solver we are in Branch and Cut
			saveSolver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo);
			// Default set of cut generators
			CglProbing generator1;
			generator1.setUsingObjective(1);
			generator1.setMaxPass(1);
			generator1.setMaxPassRoot(1);
			generator1.setMaxProbeRoot(CoinMin(3000, saveSolver->getNumCols()));
			generator1.setMaxElements(100);
			generator1.setMaxElementsRoot(200);
			generator1.setMaxLookRoot(50);
			if (saveSolver->getNumCols() > 3000)
				generator1.setMaxProbeRoot(123);
			generator1.setRowCuts(3);
			
			// Add in generators
			if ((model.moreSpecialOptions() & 65536) == 0)
				process.addCutGenerator(&generator1);
			int translate[] = { 9999, 0, 0, -3, 2, 3, -2, 9999, 4, 5 };
			process.passInMessageHandler(babModel->messageHandler());
			/* model may not have created objects
				 If none then create
			*/
			if (!model.numberObjects() && true)
				model.findIntegers(true);
			if (model.numberObjects())
			{
				OsiObject** oldObjects = babModel->objects();
				int numberOldObjects = babModel->numberObjects();
				if (!numberOldObjects)
				{
					oldObjects = model.objects();
					numberOldObjects = model.numberObjects();
				}
				// SOS
				int numberColumns = saveSolver->getNumCols();
				char* prohibited = new char[numberColumns];
				memset(prohibited, 0, numberColumns);
				int numberProhibited = 0;
				for (int iObj = 0; iObj < numberOldObjects; iObj++)
				{
					CbcSOS* obj = dynamic_cast <CbcSOS*>(oldObjects[iObj]);
					if (obj)
					{
						int n = obj->numberMembers();
						const int* which = obj->members();
						for (int i = 0; i < n; i++)
						{
							int iColumn = which[i];
							prohibited[iColumn] = 1;
							numberProhibited++;
						}
					}
					CbcLotsize* obj2 = dynamic_cast <CbcLotsize*>(oldObjects[iObj]);
					if (obj2)
					{
						int iColumn = obj2->columnNumber();
						prohibited[iColumn] = 1;
						numberProhibited++;
					}
				}
				if (numberProhibited)
					process.passInProhibited(prohibited, numberColumns);
				delete[] prohibited;
			}
			int numberPasses = 10;

			CglPreProcess* cbcPreProcessPointer = NULL;
			{
				OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (saveSolver);
				osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
				int savePerturbation = osiclp->getModelPtr()->perturbation();
				cbcPreProcessPointer = &process;
				int saveOptions = osiclp->getModelPtr()->moreSpecialOptions();
				solver2 = process.preProcessNonDefault(*saveSolver, translate[preProcess], numberPasses,
					tunePreProcess);
				/*solver2->writeMps("after");
				  saveSolver->writeMps("before");*/
				osiclp->getModelPtr()->setPerturbation(savePerturbation);
				osiclp->getModelPtr()->setMoreSpecialOptions(saveOptions);
			}

			integersOK = false; // We need to redo if CbcObjects exist
			// Tell solver we are not in Branch and Cut
			saveSolver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
			if (solver2)
				solver2->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
		}

		if (!solver2)
		{
			// say infeasible for solution
			integerStatus = 6;
			delete saveSolver;
			saveSolver = NULL;
			model.setProblemStatus(0);
			model.setSecondaryStatus(1);
			babModel->setProblemStatus(0);
			babModel->setSecondaryStatus(1);
			return -999;
		}
		else
		{
			statistics_nprocessedrows = solver2->getNumRows();
			statistics_nprocessedcols = solver2->getNumCols();
			model.setProblemStatus(-1);
			babModel->setProblemStatus(-1);
		}

		{
			// look at new integers
			int numberOriginalColumns = process.originalModel()->getNumCols();
			const int* originalColumns = process.originalColumns();
			OsiClpSolverInterface* osiclp2 = dynamic_cast<OsiClpSolverInterface*> (solver2);
			int numberColumns = osiclp2->getNumCols();
			OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (saveSolver);
			for (int i = 0; i < numberColumns; i++)
			{
				int iColumn = originalColumns[i];
				if (iColumn < numberOriginalColumns)
				{
					if (osiclp2->isInteger(i) && !osiclp->isInteger(iColumn))
						osiclp2->setOptionalInteger(i); // say optional
				}
			}
		}
		// we have to keep solver2 so pass clone
		solver2 = solver2->clone();
		// see if extra variables wanted
		int threshold =
			parameters[whichParam(CBC_PARAM_INT_EXTRA_VARIABLES, numberParameters, parameters)].intValue();
		int more2 = parameters[whichParam(CBC_PARAM_INT_MOREMOREMIPOPTIONS, numberParameters, parameters)].intValue();

		babModel->assignSolver(solver2);
		babModel->setOriginalColumns(process.originalColumns(),
			truncateColumns);
		babModel->initialSolve();
		babModel->setMaximumSeconds(timeLeft - (CoinCpuTime() - time2));
	}

	CglStored storedAmpl;
	CoinModel* coinModel = NULL;
	CoinModel saveCoinModel;
	CoinModel saveTightenedModel;
	int* whichColumn = NULL;
	int* knapsackStart = NULL;
	int* knapsackRow = NULL;
	int numberKnapsack = 0;

	int testOsiOptions = parameters[whichParam(CBC_PARAM_INT_TESTOSI, numberParameters, parameters)].intValue();

	int useCosts = 0;

	// 设置启发式
	doHeuristics(babModel, 1, parameters,
		numberParameters, noPrinting, initialPumpTune);

	int experimentFlag = parameters[whichParam(CBC_PARAM_INT_EXPERIMENT, numberParameters,
		parameters)].intValue();
	int strategyFlag = parameters[whichParam(CBC_PARAM_INT_STRATEGY, numberParameters,
		parameters)].intValue();
	int bothFlags = CoinMax(CoinMin(experimentFlag, 1), strategyFlag);
	// add cut generators if wanted
	int switches[30];
	int accuracyFlag[30];
	char doAtEnd[30];
	memset(doAtEnd, 0, 30);
	int numberGenerators = 0;
	int translate[] = { -100, -1, -99, -98, 1, -1098, -999, 1, 1, 1, -1 };
	int maximumSlowPasses =
		parameters[whichParam(CBC_PARAM_INT_MAX_SLOW_CUTS,
			numberParameters, parameters)].intValue();

	// 开始设置割平面
	// probing
	if (probingAction) { //probingAction=3
		int numberColumns = babModel->solver()->getNumCols();
		if (probingAction > 7) {
			probingGen.setMaxElements(numberColumns);
			probingGen.setMaxElementsRoot(numberColumns);
		}
		probingGen.setMaxProbeRoot(CoinMin(2000, numberColumns));
		probingGen.setMaxProbeRoot(123);
		probingGen.setMaxProbe(123);
		probingGen.setMaxLookRoot(20);
		if (probingAction == 7 || probingAction == 9)
			probingGen.setRowCuts(-3); // strengthening etc just at root
		if (probingAction == 8 || probingAction == 9) {
			// Number of unsatisfied variables to look at
			probingGen.setMaxProbeRoot(numberColumns);
			probingGen.setMaxProbe(numberColumns);
			// How far to follow the consequences
			probingGen.setMaxLook(50);
			probingGen.setMaxLookRoot(50);
		}
		if (probingAction == 10) {
			probingGen.setMaxPassRoot(2);
			probingGen.setMaxProbeRoot(numberColumns);
			probingGen.setMaxLookRoot(numberColumns);
		}
		// If 5 then force on
		int iAction = translate[probingAction];
		if (probingAction == 5)
			iAction = 1;
		babModel->addCutGenerator(&probingGen, iAction, "Probing");
		accuracyFlag[numberGenerators] = 5;
		switches[numberGenerators++] = 0;
	}
	// gomory
	if (gomoryAction && (complicatedInteger != 1 ||   //gomoryAction=3
		(gomoryAction == 1 || gomoryAction >= 4))) {
		// try larger limit
		int numberColumns = babModel->getNumCols();
		if (gomoryAction == 7) {
			gomoryAction = 4;
			gomoryGen.setLimitAtRoot(numberColumns);
			gomoryGen.setLimit(numberColumns);
		}
		else if (gomoryAction == 8) {
			gomoryAction = 3;
			gomoryGen.setLimitAtRoot(numberColumns);
			gomoryGen.setLimit(200);
		}
		else if (numberColumns > 5000) {
			//#define MORE_CUTS2
#ifdef MORE_CUTS2
									// try larger limit
			gomoryGen.setLimitAtRoot(numberColumns);
			gomoryGen.setLimit(200);
#else
			gomoryGen.setLimitAtRoot(2000);
			//gomoryGen.setLimit(200);
#endif
		}
		else {
#ifdef MORE_CUTS2
			// try larger limit
			gomoryGen.setLimitAtRoot(numberColumns);
			gomoryGen.setLimit(200);
#endif
		}
		int cutLength =
			parameters[whichParam(CBC_PARAM_INT_CUTLENGTH, numberParameters, parameters)].intValue();
		// 不执行
		if (cutLength != -1) {
			gomoryGen.setLimitAtRoot(cutLength);
			if (cutLength < 10000000) {
				gomoryGen.setLimit(cutLength);
			}
			else {
				gomoryGen.setLimit(cutLength % 10000000);
			}
		}
		int laGomory = parameters[whichParam(CBC_PARAM_STR_LAGOMORYCUTS, numberParameters, parameters)].currentOptionAsInteger();
		int gType = translate[gomoryAction];
		if (!laGomory) {
			// Normal
			babModel->addCutGenerator(&gomoryGen, translate[gomoryAction], "Gomory");
			accuracyFlag[numberGenerators] = 3;
			switches[numberGenerators++] = 0;
		}
		else {
			laGomory--;
			int type = (laGomory % 3) + 1;
			int when = laGomory / 3;
			char atEnd = (when < 2) ? 1 : 0;
			int gomoryTypeMajor = 10;
			if (when < 3) {
				// normal as well
				babModel->addCutGenerator(&gomoryGen, gType, "Gomory");
				accuracyFlag[numberGenerators] = 3;
				switches[numberGenerators++] = 0;
				if (when == 2)
					gomoryTypeMajor = 20;
			}
			else {
				when--; // so on
				gomoryTypeMajor = 20;
			}
			if (!when)
				gType = -99; // root
			gomoryGen.passInOriginalSolver(babModel->solver());
			if ((type & 1) != 0) {
				// clean
				gomoryGen.setGomoryType(gomoryTypeMajor + 1);
				babModel->addCutGenerator(&gomoryGen, gType, "GomoryL1");
				accuracyFlag[numberGenerators] = 3;
				doAtEnd[numberGenerators] = atEnd;
				if (atEnd) {
					babModel->cutGenerator(numberGenerators)->setMaximumTries(99999999);
					babModel->cutGenerator(numberGenerators)->setHowOften(1);
				}
				switches[numberGenerators++] = 0;
			}
			if ((type & 2) != 0) {
				// simple
				gomoryGen.setGomoryType(gomoryTypeMajor + 2);
				babModel->addCutGenerator(&gomoryGen, gType, "GomoryL2");
				accuracyFlag[numberGenerators] = 3;
				doAtEnd[numberGenerators] = atEnd;
				if (atEnd) {
					babModel->cutGenerator(numberGenerators)->setMaximumTries(99999999);
					babModel->cutGenerator(numberGenerators)->setHowOften(1);
				}
				switches[numberGenerators++] = 0;
			}
		}
	}
	// Knapsack
	if (knapsackAction) {//3
		babModel->addCutGenerator(&knapsackGen, translate[knapsackAction], "Knapsack");
		accuracyFlag[numberGenerators] = 1;
		switches[numberGenerators++] = -2;
	}
	// redsplit
	if (redsplitAction && !complicatedInteger) {//0
		babModel->addCutGenerator(&redsplitGen, translate[redsplitAction], "Reduce-and-split");
		accuracyFlag[numberGenerators] = 5;
		// slow ? - just do a few times
		if (redsplitAction != 1) {
			babModel->cutGenerator(numberGenerators)->setMaximumTries(maximumSlowPasses);
			babModel->cutGenerator(numberGenerators)->setHowOften(10);
		}
		switches[numberGenerators++] = 1;
	}
	// redsplit2
	if (redsplit2Action && !complicatedInteger) {
		int maxLength = 256;
		if (redsplit2Action > 2) {
			redsplit2Action -= 2;
			maxLength = COIN_INT_MAX;
		}
		CglRedSplit2Param& parameters = redsplit2Gen.getParam();
		parameters.setMaxNonzeroesTab(maxLength);
		babModel->addCutGenerator(&redsplit2Gen, translate[redsplit2Action], "Reduce-and-split(2)");
		accuracyFlag[numberGenerators] = 5;
		// slow ? - just do a few times
		if (redsplit2Action != 1) {
			babModel->cutGenerator(numberGenerators)->setHowOften(maximumSlowPasses);
			babModel->cutGenerator(numberGenerators)->setMaximumTries(maximumSlowPasses);
			babModel->cutGenerator(numberGenerators)->setHowOften(5);
		}

		switches[numberGenerators++] = 1;
	}
	// GMI
	if (GMIAction && !complicatedInteger) {
		if (GMIAction > 5) {
			// long
			GMIAction -= 5;
			CglGMIParam& parameters = GMIGen.getParam();
			parameters.setMaxSupportRel(1.0);
		}
		babModel->addCutGenerator(&GMIGen, translate[GMIAction], "Gomory(2)");
		if (GMIAction == 5) {
			// just at end and root
			GMIAction = 2;
			doAtEnd[numberGenerators] = 1;
			babModel->cutGenerator(numberGenerators)->setMaximumTries(99999999);
			babModel->cutGenerator(numberGenerators)->setHowOften(1);
		}
		accuracyFlag[numberGenerators] = 5;
		switches[numberGenerators++] = 0;
	}
	// clique
	if (cliqueAction) {
		babModel->addCutGenerator(&cliqueGen, translate[cliqueAction], "Clique");
		accuracyFlag[numberGenerators] = 0;
		switches[numberGenerators++] = 0;
	}
	// MixedIntegerRounding2
	if (mixedAction) {
		babModel->addCutGenerator(&mixedGen, translate[mixedAction], "MixedIntegerRounding2");
		accuracyFlag[numberGenerators] = 2;
		switches[numberGenerators++] = 0;
	}
	// FlowCover
	if (flowAction) {
		babModel->addCutGenerator(&flowGen, translate[flowAction], "FlowCover");
		accuracyFlag[numberGenerators] = 2;
		switches[numberGenerators++] = 0;
	}
	// TwoMirCuts
	if (twomirAction && (complicatedInteger != 1 ||
		(twomirAction == 1 || twomirAction >= 4))) {
		// try larger limit
		int numberColumns = babModel->getNumCols();
		if (twomirAction == 7) {
			twomirAction = 4;
			twomirGen.setMaxElements(numberColumns);
		}
		else if (numberColumns > 5000 && twomirAction == 4) {
			twomirGen.setMaxElements(2000);
		}
		int laTwomir = parameters[whichParam(CBC_PARAM_STR_LATWOMIRCUTS, numberParameters, parameters)].currentOptionAsInteger();
		int twomirType = translate[twomirAction];
		if (!laTwomir) {
			// Normal
			babModel->addCutGenerator(&twomirGen, translate[twomirAction], "TwoMirCuts");
			accuracyFlag[numberGenerators] = 4;
			switches[numberGenerators++] = 1;
		}
		else {
			laTwomir--;
			int type = (laTwomir % 3) + 1;
			int when = laTwomir / 3;
			char atEnd = (when < 2) ? 1 : 0;
			int twomirTypeMajor = 10;
			if (when < 3) {
				// normal as well
				babModel->addCutGenerator(&twomirGen, translate[twomirAction], "TwoMirCuts");
				accuracyFlag[numberGenerators] = 4;
				switches[numberGenerators++] = 1;
				if (when == 2)
					twomirTypeMajor = 10;
			}
			else {
				when--; // so on
				twomirTypeMajor = 20;
			}
			if (!when)
				twomirType = -99; // root
			twomirGen.passInOriginalSolver(babModel->solver());
			if ((type & 1) != 0) {
				// clean
				twomirGen.setTwomirType(twomirTypeMajor + 1);
				babModel->addCutGenerator(&twomirGen, twomirType, "TwoMirCutsL1");
				accuracyFlag[numberGenerators] = 4;
				doAtEnd[numberGenerators] = atEnd;
				switches[numberGenerators++] = atEnd ? 0 : 1;
			}
			if ((type & 2) != 0) {
				// simple
				twomirGen.setTwomirType(twomirTypeMajor + 2);
				babModel->addCutGenerator(&twomirGen, twomirType, "TwoMirCutsL2");
				accuracyFlag[numberGenerators] = 4;
				doAtEnd[numberGenerators] = atEnd;
				switches[numberGenerators++] = atEnd ? 0 : 1;
			}
		}
	}
	// LiftAndProject
	if (landpAction) {//0
		babModel->addCutGenerator(&landpGen, translate[landpAction], "LiftAndProject");
		accuracyFlag[numberGenerators] = 5;
		// slow ? - just do a few times
		if (landpAction != 1) {
			babModel->cutGenerator(numberGenerators)->setMaximumTries(maximumSlowPasses);
			babModel->cutGenerator(numberGenerators)->setHowOften(10);
		}
		switches[numberGenerators++] = 1;
	}
	// ResidualCapacity
	if (residualCapacityAction) {
		babModel->addCutGenerator(&residualCapacityGen, translate[residualCapacityAction], "ResidualCapacity");
		accuracyFlag[numberGenerators] = 5;
		switches[numberGenerators++] = 1;
	}
	// ZeroHalf
	if (zerohalfAction) {
		if (zerohalfAction > 4) {
			//zerohalfAction -=4;
			zerohalfGen.setFlags(1);
		}
		babModel->addCutGenerator(&zerohalfGen, translate[zerohalfAction], "ZeroHalf");
		accuracyFlag[numberGenerators] = 5;
		babModel->cutGenerator(numberGenerators)->
			setNeedsRefresh(true);
		switches[numberGenerators++] = 2;
	}

	// Say we want timings
	numberGenerators = babModel->numberCutGenerators();
	int iGenerator;
	int cutDepth =
		parameters[whichParam(CBC_PARAM_INT_CUTDEPTH, numberParameters, parameters)].intValue();
	for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
		CbcCutGenerator* generator = babModel->cutGenerator(iGenerator);
		int howOften = generator->howOften();
		if (howOften == -98 || howOften == -99 || generator->maximumTries() > 0)
			generator->setSwitchOffIfLessThan(switches[iGenerator]);
		// Use if any at root as more likely later and fairly cheap
		//if (switches[iGenerator]==-2)
		//generator->setWhetherToUse(true);
		generator->setInaccuracy(accuracyFlag[iGenerator]);
		if (doAtEnd[iGenerator]) {
			generator->setWhetherCallAtEnd(true);
			//generator->setMustCallAgain(true);
		}
		generator->setTiming(true);
		if (cutDepth >= 0)
			generator->setWhatDepth(cutDepth);
	}
	if (cutPass != -1234567)
		babModel->setMaximumCutPassesAtRoot(cutPass);

	if (!noPrinting) {
		int iLevel = parameters[log].intValue();
		if (iLevel < 0) {
			if (iLevel > -10) {
				babModel->setPrintingMode(1);
			}
			else {
				babModel->setPrintingMode(2);
				iLevel += 10;
				parameters[log].setIntValue(iLevel);
			}
			iLevel = -iLevel;
		}
		babModel->messageHandler()->setLogLevel(iLevel);
		if (babModel->getNumCols() > 2000 || babModel->getNumRows() > 1500 ||
			babModel->messageHandler()->logLevel() > 1)
			babModel->setPrintFrequency(100);
	}

	babModel->solver()->setIntParam(OsiMaxNumIterationHotStart,
		parameters[whichParam(CBC_PARAM_INT_MAXHOTITS, numberParameters, parameters)].intValue());
	OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (babModel->solver());
	// go faster stripes
	if ((osiclp->getNumRows() < 300 && osiclp->getNumCols() < 500)) {
		osiclp->setupForRepeatedUse(2, parameters[slog].intValue());
		if (bothFlags >= 1) {
			ClpSimplex* lp = osiclp->getModelPtr();
			int specialOptions = lp->specialOptions();
			lp->setSpecialOptions(specialOptions | (2048 + 4096));
		}
	}
	else {
		osiclp->setupForRepeatedUse(0, parameters[slog].intValue());
	}

	double increment = babModel->getCutoffIncrement();;
	int* changed = NULL;

	babModel->setCutoffIncrement(CoinMax(babModel->getCutoffIncrement(), increment));
	// Turn this off if you get problems
	// Used to be automatically set
	int mipOptions = parameters[whichParam(CBC_PARAM_INT_MIPOPTIONS, numberParameters, parameters)].intValue() % 10000;
	osiclp->setSpecialOptions(mipOptions);

	// probably faster to use a basis to get integer solutions
	babModel->setSpecialOptions(babModel->specialOptions() | 2);
	currentBranchModel = babModel;
	if (strategyFlag == 1) {
		// try reduced model
		babModel->setSpecialOptions(babModel->specialOptions() | 512);
	}
	{
		int extra1 = parameters[whichParam(CBC_PARAM_INT_EXTRA1, numberParameters, parameters)].intValue();
		if (extra1 != -1) {
			if (extra1 < 0) {
				if (extra1 == -7777)
					extra1 = -1;
				else if (extra1 == -8888)
					extra1 = 1;
				babModel->setWhenCuts(-extra1);
			}
			else if (extra1 < 19000) {
				babModel->setSearchStrategy(extra1);
				printf("XXXXX searchStrategy %d\n", extra1);
			}
			else {
				int n = extra1 - 20000;
				if (!n)
					n--;
				babModel->setNumberAnalyzeIterations(n);
				printf("XXXXX analyze %d\n", extra1);
			}
		}
		else if (bothFlags >= 1)
			babModel->setWhenCuts(999998);
	}

	const int* originalColumns = preProcess ? process.originalColumns() : NULL;
	OsiSolverInterface* testOsiSolver = (testOsiOptions >= 0) ? babModel->solver() : NULL;
	if (!testOsiSolver) {
		// *************************************************************
		// CbcObjects
		if (preProcess && (process.numberSOS() || babModel->numberObjects())) {
			int numberSOS = process.numberSOS();
			int numberIntegers = babModel->numberIntegers();
			/* model may not have created objects
			   If none then create
			*/
			if (!numberIntegers || !babModel->numberObjects()) {
				int type = (pseudoUp) ? 1 : 0;
				babModel->findIntegers(true, type);
				numberIntegers = babModel->numberIntegers();
				integersOK = true;
			}
			OsiObject** oldObjects = babModel->objects();
			// Do sets and priorities
			OsiObject** objects = new OsiObject * [numberSOS];
			// set old objects to have low priority
			int numberOldObjects = babModel->numberObjects();
			int numberColumns = babModel->getNumCols();
			// backward pointer to new variables
			// extend arrays in case SOS
			assert(originalColumns);
			int n = CoinMin(truncateColumns, numberColumns);
			n = originalColumns[n - 1] + 1;
			n = CoinMax(n, CoinMax(numberColumns, numberOriginalColumns));
			int* newColumn = new int[n];
			int i;
			for (i = 0; i < numberOriginalColumns; i++)
				newColumn[i] = -1;
			for (i = 0; i < CoinMin(truncateColumns, numberColumns); i++)
				newColumn[originalColumns[i]] = i;
			if (!integersOK) {
				// Change column numbers etc
				int n = 0;
				for (int iObj = 0; iObj < numberOldObjects; iObj++) {
					int iColumn = oldObjects[iObj]->columnNumber();
					if (iColumn < 0 || iColumn >= numberOriginalColumns) {
						oldObjects[n++] = oldObjects[iObj];
					}
					else {
						iColumn = newColumn[iColumn];
						if (iColumn >= 0) {
							CbcSimpleInteger* obj =
								dynamic_cast <CbcSimpleInteger*>(oldObjects[iObj]);
							if (obj) {
								obj->setColumnNumber(iColumn);
							}
							else {
								// only other case allowed is lotsizing
								CbcLotsize* obj2 =
									dynamic_cast <CbcLotsize*>(oldObjects[iObj]);
								assert(obj2);
								obj2->setModelSequence(iColumn);
							}
							oldObjects[n++] = oldObjects[iObj];
						}
						else {
							delete oldObjects[iObj];
						}
					}
				}
				babModel->setNumberObjects(n);
				numberOldObjects = n;
				babModel->zapIntegerInformation();
			}
			int nMissing = 0;
			for (int iObj = 0; iObj < numberOldObjects; iObj++) {
				if (process.numberSOS())
					oldObjects[iObj]->setPriority(numberColumns + 1);
				int iColumn = oldObjects[iObj]->columnNumber();
				if (iColumn < 0 || iColumn >= numberOriginalColumns) {
					CbcSOS* obj =
						dynamic_cast <CbcSOS*>(oldObjects[iObj]);
					if (obj) {
						int n = obj->numberMembers();
						int* which = obj->mutableMembers();
						double* weights = obj->mutableWeights();
						int nn = 0;
						for (i = 0; i < n; i++) {
							int iColumn = which[i];
							int jColumn = newColumn[iColumn];
							if (jColumn >= 0) {
								which[nn] = jColumn;
								weights[nn++] = weights[i];
							}
							else {
								nMissing++;
							}
						}
						obj->setNumberMembers(nn);
					}
					continue;
				}
				if (originalColumns)
					iColumn = originalColumns[iColumn];
				if (branchDirection) {
					CbcSimpleInteger* obj =
						dynamic_cast <CbcSimpleInteger*>(oldObjects[iObj]);
					if (obj) {
						obj->setPreferredWay(branchDirection[iColumn]);
					}
					else {
						CbcObject* obj =
							dynamic_cast <CbcObject*>(oldObjects[iObj]);
						assert(obj);
						obj->setPreferredWay(branchDirection[iColumn]);
					}
				}
				if (pseudoUp) {
					CbcSimpleIntegerPseudoCost* obj1a =
						dynamic_cast <CbcSimpleIntegerPseudoCost*>(oldObjects[iObj]);
					assert(obj1a);
					if (pseudoDown[iColumn] > 0.0)
						obj1a->setDownPseudoCost(pseudoDown[iColumn]);
					if (pseudoUp[iColumn] > 0.0)
						obj1a->setUpPseudoCost(pseudoUp[iColumn]);
				}
			}
			if (nMissing) {
				printf("%d SOS variables vanished due to pre processing? - check validity?", nMissing);
			}
			delete[] newColumn;
			const int* starts = process.startSOS();
			const int* which = process.whichSOS();
			const int* type = process.typeSOS();
			const double* weight = process.weightSOS();
			int iSOS;
			for (iSOS = 0; iSOS < numberSOS; iSOS++) {
				int iStart = starts[iSOS];
				int n = starts[iSOS + 1] - iStart;
				//#define MAKE_SOS_CLIQUES
#ifndef MAKE_SOS_CLIQUES
				objects[iSOS] = new CbcSOS(babModel, n, which + iStart, weight + iStart,
					iSOS, type[iSOS]);
#else
				objects[iSOS] =
					new CbcClique(babModel_, 1, n, which + iStart,
						NULL, -iSOS - 1);
#endif
				// branch on long sets first
				objects[iSOS]->setPriority(numberColumns - n);
			}
			if (numberSOS)
				babModel->addObjects(numberSOS, objects);
			for (iSOS = 0; iSOS < numberSOS; iSOS++)
				delete objects[iSOS];
			delete[] objects;
		}
		else if (priorities || branchDirection || pseudoDown || pseudoUp || numberSOS) {
			// do anyway for priorities etc
			int numberIntegers = babModel->numberIntegers();
			/* model may not have created objects
			   If none then create
			*/
			if (!numberIntegers || !babModel->numberObjects()) {
				int type = (pseudoUp) ? 1 : 0;
				babModel->findIntegers(true, type);
			}
			if (numberSOS) {
				// Do sets and priorities
				OsiObject** objects = new OsiObject * [numberSOS];
				int iSOS;
				if (originalColumns) {
					// redo sequence numbers
					int numberColumns = babModel->getNumCols();
					int nOld = originalColumns[numberColumns - 1] + 1;
					int* back = new int[nOld];
					int i;
					for (i = 0; i < nOld; i++)
						back[i] = -1;
					for (i = 0; i < numberColumns; i++)
						back[originalColumns[i]] = i;
					// Really need better checks
					int nMissing = 0;
					int n = sosStart[numberSOS];
					for (i = 0; i < n; i++) {
						int iColumn = sosIndices[i];
						int jColumn = back[iColumn];
						if (jColumn >= 0)
							sosIndices[i] = jColumn;
						else
							nMissing++;
					}
					delete[] back;
					if (nMissing) {
						printf("%d SOS variables vanished due to pre processing? - check validity?", nMissing);
					}
				}
				for (iSOS = 0; iSOS < numberSOS; iSOS++) {
					int iStart = sosStart[iSOS];
					int n = sosStart[iSOS + 1] - iStart;
					objects[iSOS] = new CbcSOS(babModel, n, sosIndices + iStart, sosReference + iStart,
						iSOS, sosType[iSOS]);
					if (sosPriority)
						objects[iSOS]->setPriority(sosPriority[iSOS]);
					else if (!prioritiesIn)
						objects[iSOS]->setPriority(10);  // rather than 1000
				}
				// delete any existing SOS objects
				int numberObjects = babModel->numberObjects();
				OsiObject** oldObjects = babModel->objects();
				int nNew = 0;
				for (int i = 0; i < numberObjects; i++) {
					OsiObject* objThis = oldObjects[i];
					CbcSOS* obj1 =
						dynamic_cast <CbcSOS*>(objThis);
					OsiSOS* obj2 =
						dynamic_cast <OsiSOS*>(objThis);
					if (!obj1 && !obj2) {
						oldObjects[nNew++] = objThis;
					}
					else {
						delete objThis;
					}
				}
				babModel->setNumberObjects(nNew);
				babModel->addObjects(numberSOS, objects);
				for (iSOS = 0; iSOS < numberSOS; iSOS++)
					delete objects[iSOS];
				delete[] objects;
			}
		}
		OsiObject** objects = babModel->objects();
		int numberObjects = babModel->numberObjects();
		for (int iObj = 0; iObj < numberObjects; iObj++) {
			// skip sos
			CbcSOS* objSOS =
				dynamic_cast <CbcSOS*>(objects[iObj]);
			if (objSOS)
				continue;
#ifdef MAKE_SOS_CLIQUES
			// skip cliques
			CbcClique* objClique =
				dynamic_cast <CbcClique*>(objects[iObj]);
			if (objClique)
				continue;
#endif
			int iColumn = objects[iObj]->columnNumber();
			assert(iColumn >= 0);
			if (originalColumns)
				iColumn = originalColumns[iColumn];
			if (branchDirection) {
				CbcSimpleInteger* obj =
					dynamic_cast <CbcSimpleInteger*>(objects[iObj]);
				if (obj) {
					obj->setPreferredWay(branchDirection[iColumn]);
				}
				else {
					CbcObject* obj =
						dynamic_cast <CbcObject*>(objects[iObj]);
					assert(obj);
					obj->setPreferredWay(branchDirection[iColumn]);
				}
			}
			if (priorities) {
				int iPriority = priorities[iColumn];
				if (iPriority > 0)
					objects[iObj]->setPriority(iPriority);
			}
			if (pseudoUp && pseudoUp[iColumn]) {
				CbcSimpleIntegerPseudoCost* obj1a =
					dynamic_cast <CbcSimpleIntegerPseudoCost*>(objects[iObj]);
				assert(obj1a);
				if (pseudoDown[iColumn] > 0.0)
					obj1a->setDownPseudoCost(pseudoDown[iColumn]);
				if (pseudoUp[iColumn] > 0.0)
					obj1a->setUpPseudoCost(pseudoUp[iColumn]);
			}
		}
		// *************************************************************
	}
	else {
		// *************************************************************
		// OsiObjects
		// Find if none
		int numberIntegers = testOsiSolver->getNumIntegers();
		/* model may not have created objects
		   If none then create
		*/
		if (!numberIntegers || !testOsiSolver->numberObjects()) {
			//int type = (pseudoUp) ? 1 : 0;
			testOsiSolver->findIntegers(false);
			numberIntegers = testOsiSolver->getNumIntegers();
		}
		if (preProcess && process.numberSOS()) {
			int numberSOS = process.numberSOS();
			OsiObject** oldObjects = testOsiSolver->objects();
			// Do sets and priorities
			OsiObject** objects = new OsiObject * [numberSOS];
			// set old objects to have low priority
			int numberOldObjects = testOsiSolver->numberObjects();
			int numberColumns = testOsiSolver->getNumCols();
			for (int iObj = 0; iObj < numberOldObjects; iObj++) {
				oldObjects[iObj]->setPriority(numberColumns + 1);
				int iColumn = oldObjects[iObj]->columnNumber();
				assert(iColumn >= 0);
				if (iColumn >= numberOriginalColumns)
					continue;
				if (originalColumns)
					iColumn = originalColumns[iColumn];
				if (branchDirection) {
					OsiSimpleInteger* obj =
						dynamic_cast <OsiSimpleInteger*>(oldObjects[iObj]);
					if (obj) {
						obj->setPreferredWay(branchDirection[iColumn]);
					}
					else {
						OsiObject2* obj =
							dynamic_cast <OsiObject2*>(oldObjects[iObj]);
						if (obj)
							obj->setPreferredWay(branchDirection[iColumn]);
					}
				}
				if (pseudoUp) {
					abort();
				}
			}
			const int* starts = process.startSOS();
			const int* which = process.whichSOS();
			const int* type = process.typeSOS();
			const double* weight = process.weightSOS();
			int iSOS;
			for (iSOS = 0; iSOS < numberSOS; iSOS++) {
				int iStart = starts[iSOS];
				int n = starts[iSOS + 1] - iStart;
				objects[iSOS] = new OsiSOS(testOsiSolver, n, which + iStart, weight + iStart,
					type[iSOS]);
				// branch on long sets first
				objects[iSOS]->setPriority(numberColumns - n);
			}
			testOsiSolver->addObjects(numberSOS, objects);
			for (iSOS = 0; iSOS < numberSOS; iSOS++)
				delete objects[iSOS];
			delete[] objects;
		}
		else if (priorities || branchDirection || pseudoDown || pseudoUp || numberSOS) {
			if (numberSOS) {
				// Do sets and priorities
				OsiObject** objects = new OsiObject * [numberSOS];
				int iSOS;
				if (originalColumns) {
					// redo sequence numbers
					int numberColumns = testOsiSolver->getNumCols();
					int nOld = originalColumns[numberColumns - 1] + 1;
					int* back = new int[nOld];
					int i;
					for (i = 0; i < nOld; i++)
						back[i] = -1;
					for (i = 0; i < numberColumns; i++)
						back[originalColumns[i]] = i;
					// Really need better checks
					int nMissing = 0;
					int n = sosStart[numberSOS];
					for (i = 0; i < n; i++) {
						int iColumn = sosIndices[i];
						int jColumn = back[iColumn];
						if (jColumn >= 0)
							sosIndices[i] = jColumn;
						else
							nMissing++;
					}
					delete[] back;
					if (nMissing) {
						printf("%d SOS variables vanished due to pre processing? - check validity?", nMissing);
					}
				}
				for (iSOS = 0; iSOS < numberSOS; iSOS++) {
					int iStart = sosStart[iSOS];
					int n = sosStart[iSOS + 1] - iStart;
					objects[iSOS] = new OsiSOS(testOsiSolver, n, sosIndices + iStart, sosReference + iStart,
						sosType[iSOS]);
					if (sosPriority)
						objects[iSOS]->setPriority(sosPriority[iSOS]);
					else if (!prioritiesIn)
						objects[iSOS]->setPriority(10);  // rather than 1000
				}
				// delete any existing SOS objects
				int numberObjects = testOsiSolver->numberObjects();
				OsiObject** oldObjects = testOsiSolver->objects();
				int nNew = 0;
				for (int i = 0; i < numberObjects; i++) {
					OsiObject* objThis = oldObjects[i];
					OsiSOS* obj1 =
						dynamic_cast <OsiSOS*>(objThis);
					OsiSOS* obj2 =
						dynamic_cast <OsiSOS*>(objThis);
					if (!obj1 && !obj2) {
						oldObjects[nNew++] = objThis;
					}
					else {
						delete objThis;
					}
				}
				testOsiSolver->setNumberObjects(nNew);
				testOsiSolver->addObjects(numberSOS, objects);
				for (iSOS = 0; iSOS < numberSOS; iSOS++)
					delete objects[iSOS];
				delete[] objects;
			}
		}
		OsiObject** objects = testOsiSolver->objects();
		int numberObjects = testOsiSolver->numberObjects();
		int logLevel = parameters[log].intValue();
		for (int iObj = 0; iObj < numberObjects; iObj++) {
			// skip sos
			OsiSOS* objSOS =
				dynamic_cast <OsiSOS*>(objects[iObj]);
			if (objSOS) {
				if (logLevel > 2)
					printf("Set %d is SOS - priority %d\n", iObj, objSOS->priority());
				continue;
			}
			int iColumn = objects[iObj]->columnNumber();
			if (iColumn >= 0) {
				if (originalColumns)
					iColumn = originalColumns[iColumn];
				if (branchDirection) {
					OsiSimpleInteger* obj =
						dynamic_cast <OsiSimpleInteger*>(objects[iObj]);
					if (obj) {
						obj->setPreferredWay(branchDirection[iColumn]);
					}
					else {
						OsiObject2* obj =
							dynamic_cast <OsiObject2*>(objects[iObj]);
						if (obj)
							obj->setPreferredWay(branchDirection[iColumn]);
					}
				}
				if (priorities) {
					int iPriority = priorities[iColumn];
					if (iPriority > 0)
						objects[iObj]->setPriority(iPriority);
				}
				if (logLevel > 2)
					printf("Obj %d is int? - priority %d\n", iObj, objects[iObj]->priority());
				if (pseudoUp && pseudoUp[iColumn]) {
					abort();
				}
			}
		}
		// *************************************************************
	}
	int statistics = (printOptions > 0) ? printOptions : 0;
	free(priorities);
	priorities = NULL;
	free(branchDirection);
	branchDirection = NULL;
	free(pseudoDown);
	pseudoDown = NULL;
	free(pseudoUp);
	pseudoUp = NULL;
	free(solutionIn);
	solutionIn = NULL;
	free(prioritiesIn);
	prioritiesIn = NULL;
	free(sosStart);
	sosStart = NULL;
	free(sosIndices);
	sosIndices = NULL;
	free(sosType);
	sosType = NULL;
	free(sosReference);
	sosReference = NULL;
	free(cut);
	cut = NULL;
	free(sosPriority);
	sosPriority = NULL;

	if (!babModel->numberStrong() && babModel->numberBeforeTrust() > 0)
		babModel->setNumberBeforeTrust(0);

	// If defaults then increase trust for small models
	if (!strongChanged) {
		int numberColumns = babModel->getNumCols();
		if (numberColumns <= 50)
			babModel->setNumberBeforeTrust(1000);
		else if (numberColumns <= 100)
			babModel->setNumberBeforeTrust(100);
		else if (numberColumns <= 300)
			babModel->setNumberBeforeTrust(50);
	}

	osiclp = dynamic_cast<OsiClpSolverInterface*> (babModel->solver());
	lpSolver = osiclp->getModelPtr();
	int hotits = parameters[whichParam(CBC_PARAM_INT_MAXHOTITS, numberParameters, parameters)].intValue();
	if (hotits > 100) {
		osiclp->setSpecialOptions(osiclp->specialOptions() & ~32);
		osiclp->setIntParam(OsiMaxNumIterationHotStart, hotits);
	}
	else {
		osiclp->setIntParam(OsiMaxNumIterationHotStart, hotits);
	}

	if ((experimentFlag >= 1 || strategyFlag >= 1) && babModel->fastNodeDepth() == -1) {
		if (babModel->solver()->getNumCols() +
			babModel->solver()->getNumRows() < 500)
			babModel->setFastNodeDepth(-12);
	}
	else if (babModel->fastNodeDepth() == -999) {
		babModel->setFastNodeDepth(-1);
	}
	int heurOptions = parameters[whichParam(CBC_PARAM_INT_HOPTIONS, numberParameters, parameters)].intValue();
	if (heurOptions > 100)
		babModel->setSpecialOptions(babModel->specialOptions() | 8192);

	int denseCode = parameters[whichParam(CBC_PARAM_INT_DENSE, numberParameters, parameters)].intValue();
	int smallCode = parameters[whichParam(CBC_PARAM_INT_SMALLFACT, numberParameters, parameters)].intValue();
	if (bothFlags >= 1) {
		if (denseCode < 0)
			denseCode = 40;
		if (smallCode < 0 && !lpSolver->factorization()->isDenseOrSmall())
			smallCode = 40;
	}
	if (denseCode > 0) {
		lpSolver->factorization()->setGoDenseThreshold(denseCode);
		assert(osiclp == babModel->solver());
		osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
	}
	if (smallCode > 0 && smallCode > denseCode)
		lpSolver->factorization()->setGoSmallThreshold(smallCode);
	//if (denseCode>=lpSolver->numberRows()) {
	//lpSolver->factorization()->goDense();
	//}
	if (logLevel <= 1)
		babModel->solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);

	int specialOptions = parameters[whichParam(CBC_PARAM_INT_STRONG_STRATEGY, numberParameters, parameters)].intValue();
	if (specialOptions >= 0)
		babModel->setStrongStrategy(specialOptions);
	int jParam = whichParam(CBC_PARAM_STR_CUTOFF_CONSTRAINT,
		numberParameters, parameters);
	if (parameters[jParam].currentOptionAsInteger()) {
		babModel->setCutoffAsConstraint(true);
		int moreOptions = babModel->moreSpecialOptions();
		if (parameters[jParam].currentOptionAsInteger() == 4)
			babModel->setMoreSpecialOptions(moreOptions | 4194304);
	}
	int multipleRoot = parameters[whichParam(CBC_PARAM_INT_MULTIPLEROOTS, numberParameters, parameters)].intValue();
	if (multipleRoot < 10000) {
		//babModel->setMultipleRootTries(multipleRoot);
		babModel->setMultipleRootTries(0);
	}
	else {
		// will be doing repeated solves and saves
		int numberGoes = multipleRoot / 10000;
		multipleRoot -= 10000 * numberGoes;
		int moreOptions = babModel->moreSpecialOptions();
		if (numberGoes < 100) {
			remove("global.cuts");
			remove("global.fix");
			moreOptions |= (67108864 | 134217728);
		}
		else {
			moreOptions |= 67108864 * (numberGoes / 100);
			numberGoes = numberGoes % 100;
		}
		babModel->setMultipleRootTries(multipleRoot);
		babModel->setMoreSpecialOptions(moreOptions);
		int numberColumns = babModel->getNumCols();
		double* bestValues = new double[numberGoes];
		double** bestSolutions = new double* [numberGoes];
		int* which = new int[numberGoes];
		int numberSolutions = 0;
		printf("Starting %d passes each with %d solvers",
			numberGoes, multipleRoot % 10);
		for (int iGo = 0; iGo < numberGoes; iGo++) {
			printf("Starting pass %d", iGo + 1);
			CbcModel tempModel = *babModel;
			tempModel.setMaximumNodes(0);
			// switch off cuts if none generated
			int numberGenerators = tempModel.numberCutGenerators();
			for (int iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
				CbcCutGenerator* generator = tempModel.cutGenerator(iGenerator);
				generator->setSwitchOffIfLessThan(1);
			}
			// random
			tempModel.setRandomSeed(tempModel.getRandomSeed() + 100000000 * (iGo + 1 + 5 * numberGoes));
			for (int i = 0; i < tempModel.numberHeuristics(); i++)
				tempModel.heuristic(i)->setSeed(tempModel.heuristic(i)->getSeed() + 100000000 * iGo);
#ifndef CBC_OTHER_SOLVER
			OsiClpSolverInterface* solver = dynamic_cast<OsiClpSolverInterface*> (tempModel.solver());
			ClpSimplex* simplex = solver->getModelPtr();
			int solverSeed = simplex->randomNumberGenerator()->getSeed();
			simplex->setRandomSeed(solverSeed + 100000000 * (iGo + 1));
#endif
			tempModel.branchAndBound();
			if (tempModel.bestSolution()) {
				bestSolutions[numberSolutions] =
					CoinCopyOfArray(tempModel.bestSolution(),
						numberColumns);
				bestValues[numberSolutions] = -tempModel.getMinimizationObjValue();
				which[numberSolutions] = numberSolutions;
				numberSolutions++;
			}
		}
		// allow solutions
		double sense = babModel->solver()->getObjSense();;
		CoinSort_2(bestValues, bestValues + numberSolutions, which);
		babModel->setMoreSpecialOptions(moreOptions & (~16777216));
		for (int i = 0; i < numberSolutions; i++) {
			int k = which[i];
			if (bestValues[i] < babModel->getCutoff()) {
				babModel->setBestSolution(bestSolutions[k], numberColumns,
					-bestValues[i] * sense, true);
				babModel->incrementUsed(bestSolutions[k]);
			}
			delete[] bestSolutions[k];
		}
		babModel->setMoreSpecialOptions(moreOptions);
		if (numberSolutions)
			printf("Ending major passes - best solution %g", -bestValues[numberSolutions - 1]);
		else
			printf("Ending major passes - no solution found");
		delete[] which;
		delete[] bestValues;
		delete[] bestSolutions;
	}
	if (biLinearProblem)
		babModel->setSpecialOptions(babModel->specialOptions() & (~(512 | 32768)));
	babModel->setMoreSpecialOptions2(parameters[whichParam(CBC_PARAM_INT_MOREMOREMIPOPTIONS, numberParameters, parameters)].intValue());
	
	// 分支定界
	babModel->branchAndBound(statistics);
	if (truncateColumns < babModel->solver()->getNumCols()) {
		OsiSolverInterface* solverX = babModel->solver();
		int numberColumns = solverX->getNumCols();
		int numberRows = solverX->getNumRows();
		int numberDelete = numberColumns - truncateColumns;
		int* delStuff = new int[numberDelete];
		for (int i = 0; i < numberDelete; i++)
			delStuff[i] = i + truncateColumns;
		solverX->deleteCols(numberDelete, delStuff);
		numberDelete = numberRows - truncateRows;
		for (int i = 0; i < numberDelete; i++)
			delStuff[i] = i + truncateRows;
		solverX->deleteRows(numberDelete, delStuff);
		delete[] delStuff;
		if (truncatedRhsLower) {
			numberRows = solverX->getNumRows();
			for (int i = 0; i < numberRows; i++) {
				solverX->setRowLower(i, truncatedRhsLower[i]);
				solverX->setRowUpper(i, truncatedRhsUpper[i]);
			}
			delete[] truncatedRhsLower;
			delete[] truncatedRhsUpper;
		}
	}

	int numberSolutions = babModel->numberSavedSolutions();
	if (numberSolutions > 1) {
		for (int iSolution = numberSolutions - 1; iSolution >= 0; iSolution--) {
			model.setBestSolution(babModel->savedSolution(iSolution),
				model.solver()->getNumCols(),
				babModel->savedSolutionObjective(iSolution));
		}
	}
	currentBranchModel = NULL;
	osiclp = dynamic_cast<OsiClpSolverInterface*> (babModel->solver());
	statistics_cut_time = 0.0;
	if (!noPrinting) {
		// Print more statistics
		printf("Cuts at root node changed objective from %g to %g",
			babModel->getContinuousObjective(), babModel->rootObjectiveAfterCuts());

		numberGenerators = babModel->numberCutGenerators();
		statistics_number_cuts = new int[numberGenerators];;
		statistics_number_generators = numberGenerators;
		statistics_name_generators = new const char* [numberGenerators];
		char timing[30];
		for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
			CbcCutGenerator* generator = babModel->cutGenerator(iGenerator);
			statistics_name_generators[iGenerator] =
				generator->cutGeneratorName();
			statistics_number_cuts[iGenerator] = generator->numberCutsInTotal();
			printf("%s was tried %d times and created %d cuts of which %d were active after adding rounds of cuts",
				generator->cutGeneratorName(),
				generator->numberTimesEntered(),
				generator->numberCutsInTotal() +
				generator->numberColumnCuts(),
				generator->numberCutsActive());
			if (generator->timing()) {
				sprintf(timing, " (%.3f seconds)", generator->timeInCutGenerator());
				statistics_cut_time += generator->timeInCutGenerator();
			}
			CglStored* stored = dynamic_cast<CglStored*>(generator->generator());
			if (stored && !generator->numberCutsInTotal())
				continue;
#ifndef CLP_INVESTIGATE
			CglImplication* implication = dynamic_cast<CglImplication*>(generator->generator());
			if (implication && !generator->numberCutsInTotal())
				continue;
#endif
		}
	}
	// adjust time to allow for children on some systems
	time2 = CoinCpuTime() + CoinCpuTimeJustChildren();
	totalTime += time2 - time1;
	// For best solution
	double* bestSolution = NULL;
	// Say in integer
	if (babModel->status()) {
		// treat as stopped
		integerStatus = 3;
	}
	else {
		if (babModel->isProvenOptimal()) {
			integerStatus = 0;
		}
		else {
			// infeasible
			integerStatus = 6;
			delete saveSolver;
			saveSolver = NULL;
		}
	}
	if (babModel->getMinimizationObjValue() < 1.0e50) {
		// post process
		int n;
		if (preProcess) {
			n = saveSolver->getNumCols();
			bestSolution = new double[n];
#ifndef CBC_OTHER_SOLVER
			OsiClpSolverInterface* clpSolver = dynamic_cast<OsiClpSolverInterface*> (babModel->solver());
#else
			OsiCpxSolverInterface* clpSolver = dynamic_cast<OsiCpxSolverInterface*> (babModel_->solver());
#endif
			// Save bounds on processed model
			const int* originalColumns = process.originalColumns();
			int numberColumns2 = clpSolver->getNumCols();
			double* solution2 = new double[n];
			double* lower2 = new double[n];
			double* upper2 = new double[n];
			for (int i = 0; i < n; i++) {
				solution2[i] = COIN_DBL_MAX;
				lower2[i] = COIN_DBL_MAX;
				upper2[i] = -COIN_DBL_MAX;
			}
			const double* columnLower = clpSolver->getColLower();
			const double* columnUpper = clpSolver->getColUpper();
			const double* solution = babModel->bestSolution();
			for (int i = 0; i < numberColumns2; i++) {
				int jColumn = originalColumns[i];
				if (jColumn < n) {
					solution2[jColumn] = solution[i];
					lower2[jColumn] = columnLower[i];
					upper2[jColumn] = columnUpper[i];
				}
			}
#ifndef CBC_OTHER_SOLVER
			ClpSimplex* lpSolver = clpSolver->getModelPtr();
			lpSolver->setSpecialOptions(lpSolver->specialOptions() | IN_BRANCH_AND_BOUND); // say is Cbc (and in branch and bound)
#endif
			// put back any saved solutions
			putBackOtherSolutions(babModel, &model, &process);
			process.postProcess(*babModel->solver());
#ifdef COIN_DEVELOP
			if (model_.bestSolution() && fabs(model_.getMinimizationObjValue() -
				babModel_->getMinimizationObjValue()) < 1.0e-8) {
				const double* b1 = model_.bestSolution();
				const double* b2 = saveSolver->getColSolution();
				const double* columnLower = saveSolver->getColLower();
				const double* columnUpper = saveSolver->getColUpper();
				for (int i = 0; i < n; i++) {
					if (fabs(b1[i] - b2[i]) > 1.0e-7) {
						printf("%d %g %g %g %g\n", i, b1[i], b2[i],
							columnLower[i], columnUpper[i]);
					}
				}
		}
#endif
			bool tightenB = false;
			{
				int n = babModel->numberObjects();
				for (int i = 0; i < n; i++) {
					const OsiObject* obj = babModel->object(i);
					if (!dynamic_cast<const CbcSimpleInteger*>(obj)) {
						tightenB = true;
						break;
					}
				}
			}
			// Solution now back in saveSolver
			// Double check bounds
			columnLower = saveSolver->getColLower();
			columnUpper = saveSolver->getColUpper();
			solution = saveSolver->getColSolution();
			int numberChanged = 0;
			for (int i = 0; i < n; i++) {
				if (!saveSolver->isInteger(i) && !tightenB)
					continue;
				if (lower2[i] != COIN_DBL_MAX) {
					if (lower2[i] != columnLower[i] ||
						upper2[i] != columnUpper[i]) {
						if (lower2[i] < columnLower[i] ||
							upper2[i] > columnUpper[i]) {
#ifdef COIN_DEVELOP
							printf("odd bounds tighter");
							printf("%d bab bounds %g %g now %g %g\n",
								i, lower2[i], upper2[i], columnLower[i],
								columnUpper[i]);
#endif
					}
						else {
#ifdef COIN_DEVELOP
							printf("%d bab bounds %g %g now %g %g\n",
								i, lower2[i], upper2[i], columnLower[i],
								columnUpper[i]);
#endif
							numberChanged++;
							saveSolver->setColLower(i, lower2[i]);
							saveSolver->setColUpper(i, upper2[i]);
				}
			}
	}
}
			delete[] solution2;
			delete[] lower2;
			delete[] upper2;
			if (numberChanged) {
				printf("%d bounds tightened after postprocessing\n",
					numberChanged);
			}
			saveSolver->resolve();
			if (!saveSolver->isProvenOptimal()) {
				// try all slack
				CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel->solver()->getEmptyWarmStart());
				saveSolver->setWarmStart(basis);
				delete basis;
				saveSolver->initialSolve();
#ifdef COIN_DEVELOP
				saveSolver->writeMps("inf2");
#endif
				OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (saveSolver);
				if (osiclp)
					osiclp->getModelPtr()->checkUnscaledSolution();
			}
			assert(saveSolver->isProvenOptimal());
#ifndef CBC_OTHER_SOLVER
			// and original solver
			originalSolver->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
			assert(n >= originalSolver->getNumCols());
			n = originalSolver->getNumCols();
			originalSolver->setColLower(saveSolver->getColLower());
			originalSolver->setColUpper(saveSolver->getColUpper());
			// basis
			CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel->solver()->getWarmStart());
			originalSolver->setBasis(*basis);
			delete basis;
			originalSolver->resolve();
			if (!originalSolver->isProvenOptimal()) {
				// try all slack
				CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel->solver()->getEmptyWarmStart());
				originalSolver->setBasis(*basis);
				delete basis;
				originalSolver->initialSolve();
				OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (originalSolver);
				if (osiclp)
					osiclp->getModelPtr()->checkUnscaledSolution();
			}
			assert(originalSolver->isProvenOptimal());
#endif
			babModel->assignSolver(saveSolver);
			memcpy(bestSolution, babModel->solver()->getColSolution(), n * sizeof(double));
		}
		else {
			n = babModel->solver()->getNumCols();
			bestSolution = new double[n];
			memcpy(bestSolution, babModel->solver()->getColSolution(), n * sizeof(double));
		}
		if (returnMode == 1 && model.numberSavedSolutions() < 2) {
			model.deleteSolutions();
			model.setBestSolution(bestSolution, n, babModel->getMinimizationObjValue());
		}
		babModel->deleteSolutions();
		babModel->setBestSolution(bestSolution, n, babModel->getMinimizationObjValue());
#ifndef CBC_OTHER_SOLVER
		// and put back in very original solver
		{
			ClpSimplex* original = originalSolver->getModelPtr();
			double* lower = original->columnLower();
			double* upper = original->columnUpper();
			double* solution = original->primalColumnSolution();
			int n = original->numberColumns();
			//assert (!n||n==babModel_->solver()->getNumCols());
			for (int i = 0; i < n; i++) {
				solution[i] = bestSolution[i];
				if (originalSolver->isInteger(i)) {
					lower[i] = solution[i];
					upper[i] = solution[i];
				}
			}
			// basis
			CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel->solver()->getWarmStart());
			originalSolver->setBasis(*basis);
			delete basis;
			originalSolver->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
			originalSolver->resolve();
			if (!originalSolver->isProvenOptimal()) {
				// try all slack
				CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel->solver()->getEmptyWarmStart());
				originalSolver->setBasis(*basis);
				delete basis;
				originalSolver->initialSolve();
				OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (originalSolver);
				if (osiclp)
					osiclp->getModelPtr()->checkUnscaledSolution();
#ifdef CLP_INVESTIGATE
				if (!originalSolver->isProvenOptimal()) {
					if (saveSolver) {
						printf("saveSolver and originalSolver matrices saved\n");
						saveSolver->writeMps("infA");
					}
					else {
						printf("originalSolver matrix saved\n");
						originalSolver->writeMps("infB");
					}
			}
#endif
		}
			assert(originalSolver->isProvenOptimal());
		}
#endif
	}
	else if (model.bestSolution() && model.getMinimizationObjValue() < 1.0e50 && preProcess) {
		printf("Restoring heuristic best solution of %g", model.getMinimizationObjValue());
		int n = saveSolver->getNumCols();
		bestSolution = new double[n];
		// Put solution now back in saveSolver
		saveSolver->setColSolution(model.bestSolution());
		babModel->assignSolver(saveSolver);
		saveSolver = NULL;
		babModel->setMinimizationObjValue(model.getMinimizationObjValue());
		memcpy(bestSolution, babModel->solver()->getColSolution(), n * sizeof(double));
#ifndef CBC_OTHER_SOLVER
		// and put back in very original solver
		{
			ClpSimplex* original = originalSolver->getModelPtr();
			double* lower = original->columnLower();
			double* upper = original->columnUpper();
			double* solution = original->primalColumnSolution();
			int n = original->numberColumns();
			//assert (!n||n==babModel_->solver()->getNumCols());
			for (int i = 0; i < n; i++) {
				solution[i] = bestSolution[i];
				if (originalSolver->isInteger(i)) {
					lower[i] = solution[i];
					upper[i] = solution[i];
				}
			}
			// basis
			CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel->solver()->getWarmStart());
			originalSolver->setBasis(*basis);
			delete basis;
		}
#endif
	}

	lpSolver = clpSolver->getModelPtr();
	if (numberChanged) {
		for (int i = 0; i < numberChanged; i++) {
			int iColumn = changed[i];
			clpSolver->setContinuous(iColumn);
		}
		delete[] changed;
	}

	int n = lpSolver->getNumCols();
	if (bestSolution) {
		memcpy(lpSolver->primalColumnSolution(), bestSolution, n * sizeof(double));
		// now see what that does to row solution
		int numberRows = lpSolver->numberRows();
		double* rowSolution = lpSolver->primalRowSolution();
		memset(rowSolution, 0, numberRows * sizeof(double));
		lpSolver->clpMatrix()->times(1.0, bestSolution, rowSolution);
		lpSolver->setObjectiveValue(babModel->getObjValue());
	}

	delete saveSolver;
	delete[] bestSolution;
	std::string statusName[] = { "", "Stopped on ", "Run abandoned", "", "", "User ctrl-c" };
	std::string minor[] = { "Optimal solution found", "Linear relaxation infeasible", "Optimal solution found (within gap tolerance)", "node limit", "time limit", "user ctrl-c", "solution limit", "Linear relaxation unbounded", "Problem proven infeasible" };
	int iStat = babModel->status();
	int iStat2 = babModel->secondaryStatus();
	if (!iStat && !iStat2 && !bestSolution)
		iStat2 = 8;
	if (!iStat && iStat2 == 1 && bestSolution)
		iStat2 = 0; // solution and search completed
	statistics_seconds = time2 - time1;
	statistics_sys_seconds = CoinSysTime();
	statistics_elapsed_seconds = CoinWallclockTime();
	statistics_obj = babModel->getObjValue();
	statistics_continuous = babModel->getContinuousObjective();
	statistics_tighter = babModel->rootObjectiveAfterCuts();
	statistics_nodes = babModel->getNodeCount();
	statistics_iterations = babModel->getIterationCount();;
	statistics_result = statusName[iStat];
	if (!noPrinting) {
		printf("\nResult - %s%s\n",
			statusName[iStat].c_str(),
			minor[iStat2].c_str());
		if (babModel->bestSolution()) {
			printf("Objective value:                %.8f\n",
				babModel->getObjValue());
		}
		else {
			printf("No feasible solution found\n");
		}
		if (iStat2 >= 2 && iStat2 <= 6) {
			printf("Lower bound:                    %.3f\n",
				babModel->getBestPossibleObjValue());
			if (babModel->bestSolution()) {
				"Gap:                            %.2f\n",
					(babModel->getObjValue() - babModel->getBestPossibleObjValue()) /
					fabs(babModel->getBestPossibleObjValue());
			}
		}
		printf("Enumerated nodes:               %d\n",
			babModel->getNodeCount());
		printf("Total iterations:               %d\n",
			babModel->getIterationCount());
		printf("Time (CPU seconds):             %.2f\n",
			CoinCpuTime() - time1);
		printf("Time (Wallclock seconds):       %.2f\n",
			CoinGetTimeOfDay() - time1Elapsed);
	}
	return 0;
}