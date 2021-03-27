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
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicRINS.hpp"

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
	double timeStart = CoinCpuTime();
	double totalTime = 0.0;
	OsiClpSolverInterface initialSolver;
	CbcModel model_(initialSolver);
	/*******************cbcMain0中的一些设置**********************/
	{
		OsiSolverInterface* solver = model_.solver();
		OsiClpSolverInterface* clpSolver = dynamic_cast<OsiClpSolverInterface*> (solver);
		ClpSimplex* lpSolver = clpSolver->getModelPtr();
		lpSolver->setPerturbation(50);
		lpSolver->messageHandler()->setPrefix(false);
		clpSolver->messageHandler()->setLogLevel(1);
		model_.messageHandler()->setLogLevel(1);
		lpSolver->setLogLevel(1);
		model_.setNumberBeforeTrust(10);
		model_.setNumberStrong(5);
	}

	CbcModel* babModel_ = NULL;
	OsiClpSolverInterface* originalSolver = dynamic_cast<OsiClpSolverInterface*> (model_.solver());
	model_.setResolveAfterTakeOffCuts(false);  // Say no resolve after cuts
	OsiSolverInterface* solver = model_.solver();
	OsiClpSolverInterface* clpSolver = dynamic_cast<OsiClpSolverInterface*> (solver);
	ClpSimplex* lpSolver = clpSolver->getModelPtr();

	/**********************读取数据**************************/
	{
		std::string fileName = "uc10.mps";
		ClpSimplex* lpSolver = clpSolver->getModelPtr();
		int status = clpSolver->readMps(fileName.c_str());
		if (status != 0)
		{
			printf("%d errors reading data file\n", status);
			return status;
		}
		int numberColumns = lpSolver->numberColumns();
		for (int i = 0; i < numberColumns; i++)
		{
			if (lpSolver->isInteger(i))
				clpSolver->setInteger(i);
		}
	}
	// 
	{
		OsiSolverInterface* solver = model_.solver();
		OsiClpSolverInterface* si = dynamic_cast<OsiClpSolverInterface*>(solver);
		si->getModelPtr()->scaling(4);
		si->setHintParam(OsiDoReducePrint, true, OsiHintTry);
		si->setSpecialOptions(0x40000000);
	}
	double start = CoinCpuTime();
	/***********求解LP松弛************/
	{
		double time0 = CoinCpuTime();
		OsiSolverInterface* solver = model_.solver();
		OsiClpSolverInterface* si = dynamic_cast<OsiClpSolverInterface*>(solver);
		if (si)
			si->setSpecialOptions(si->specialOptions() | 1024);
		model_.initialSolve();
		ClpSimplex* clpSolver = si->getModelPtr();
		int iStatus = clpSolver->status();
		int iStatus2 = clpSolver->secondaryStatus();
		if (iStatus == 0)
			iStatus2 = 0;
		else if (iStatus == 1)
		{
			iStatus = 0;
			iStatus2 = 1; // say infeasible
		}
		else if (iStatus == 2)
		{
			iStatus = 0;
			iStatus2 = 7; // say unbounded
		}
		else if (iStatus == 3)
		{
			iStatus = 1;
			if (iStatus2 == 9)
				iStatus2 = 4;
			else
				iStatus2 = 3; // Use nodes - as closer than solutions
		}
		else if (iStatus == 4)
		{
			iStatus = 2; // difficulties
			iStatus2 = 0;
		}
		model_.setProblemStatus(iStatus);
		model_.setSecondaryStatus(iStatus2);
		si->setWarmStart(NULL);
		clpSolver->setSpecialOptions(clpSolver->specialOptions() | IN_BRANCH_AND_BOUND);
		if (iStatus > 0)
		{
			const char* msg[] = { "infeasible", "unbounded", "stopped","difficulties", "other" };
			printf("Problem is %s - %.2f seconds", msg[iStatus - 1], CoinCpuTime() - time0);
			return -9999;
		}
		else
		{
			printf("Continuous objective value is %g - %.2f seconds\n",
				solver->getObjValue(), CoinCpuTime() - time0);
		}
		// 设置dualBound
		if (clpSolver->dualBound() == 1.0e10)
		{
			ClpSimplex temp = *clpSolver;
			temp.setLogLevel(0);
			temp.dual(0, 7);
			// user did not set - so modify
			// get largest scaled away from bound
			double largest = 1.0e-12;
			double largestScaled = 1.0e-12;
			int numberRows = temp.numberRows();
			const double* rowPrimal = temp.primalRowSolution();
			const double* rowLower = temp.rowLower();
			const double* rowUpper = temp.rowUpper();
			const double* rowScale = temp.rowScale();
			int iRow;
			for (iRow = 0; iRow < numberRows; iRow++)
			{
				double value = rowPrimal[iRow];
				double above = value - rowLower[iRow];
				double below = rowUpper[iRow] - value;
				if (above < 1.0e12)
					largest = CoinMax(largest, above);
				if (below < 1.0e12)
					largest = CoinMax(largest, below);
				if (rowScale)
				{
					double multiplier = rowScale[iRow];
					above *= multiplier;
					below *= multiplier;
				}
				if (above < 1.0e12)
					largestScaled = CoinMax(largestScaled, above);
				if (below < 1.0e12)
					largestScaled = CoinMax(largestScaled, below);
			}
			int numberColumns = temp.numberColumns();
			const double* columnPrimal = temp.primalColumnSolution();
			const double* columnLower = temp.columnLower();
			const double* columnUpper = temp.columnUpper();
			const double* columnScale = temp.columnScale();
			int iColumn;
			for (iColumn = 0; iColumn < numberColumns; iColumn++)
			{
				double value = columnPrimal[iColumn];
				double above = value - columnLower[iColumn];
				double below = columnUpper[iColumn] - value;
				if (above < 1.0e12)
					largest = CoinMax(largest, above);
				if (below < 1.0e12)
					largest = CoinMax(largest, below);
				if (columnScale)
				{
					double multiplier = 1.0 / columnScale[iColumn];
					above *= multiplier;
					below *= multiplier;
				}
				if (above < 1.0e12)
					largestScaled = CoinMax(largestScaled, above);
				if (below < 1.0e12)
					largestScaled = CoinMax(largestScaled, below);
			}
			clpSolver->setDualBound(CoinMax(1.0001e8, CoinMin(100.0 * largest, 1.00001e10)));
		}
		si->resolve();  // clean up
	}

	/*************预处理**************/
	OsiSolverInterface* saveSolver = NULL;
	CglPreProcess process;
	babModel_ = new CbcModel(model_);
	OsiSolverInterface* solver3 = clpSolver->clone();
	babModel_->assignSolver(solver3);
	OsiClpSolverInterface* clpSolver2 = dynamic_cast<OsiClpSolverInterface*> (babModel_->solver());
	lpSolver = clpSolver2->getModelPtr();

	// 设置因式分解频率，对求解速度影响较大
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

	bool preProcess = true;  // 标志是否进行预处理，对于我们这个算例需要进行
	{
		// see whether to switch off preprocessing
		// only allow SOS and integer
		OsiObject** objects = babModel_->objects();
		int numberObjects = babModel_->numberObjects();
		for (int iObj = 0; iObj < numberObjects; iObj++)
		{
			CbcSOS* objSOS = dynamic_cast <CbcSOS*>(objects[iObj]);
			CbcSimpleInteger* objSimpleInteger = dynamic_cast <CbcSimpleInteger*>(objects[iObj]);
			if (!objSimpleInteger && !objSOS)
			{
				// find all integers anyway
				babModel_->findIntegers(true);
				preProcess = false;
				break;
			}
		}
	}
	{
		double limit;
		clpSolver->getDblParam(OsiDualObjectiveLimit, limit);
		if (clpSolver->getObjValue() * clpSolver->getObjSense() >=
			limit * clpSolver->getObjSense())
			preProcess = false;
	}
	if (preProcess)
	{
		saveSolver = babModel_->solver()->clone();
		OsiSolverInterface* solver2;
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
			if ((babModel_->specialOptions() & 65536) != 0)
				process.setOptions(1);
			// Add in generators
			if ((model_.moreSpecialOptions() & 65536) == 0)
				process.addCutGenerator(&generator1);
			process.passInMessageHandler(babModel_->messageHandler());
			if (!model_.numberObjects())
			{
				/* model may not have created objects
				   If none then create
				*/
				model_.findIntegers(true);
			}
			if (model_.numberObjects())
			{
				OsiObject** oldObjects = babModel_->objects();
				int numberOldObjects = babModel_->numberObjects();
				if (!numberOldObjects)
				{
					oldObjects = model_.objects();
					numberOldObjects = model_.numberObjects();
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
			int tunePreProcess = 6;
			{
				OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (saveSolver);
				osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
				int savePerturbation = osiclp->getModelPtr()->perturbation();
				if ((model_.moreSpecialOptions() & 65536) != 0)
					process.setOptions(2 + 4 + 8); // no cuts
				int saveOptions = osiclp->getModelPtr()->moreSpecialOptions();
				if ((model_.specialOptions() & 16777216) != 0 &&
					model_.getCutoff() > 1.0e30) {
					osiclp->getModelPtr()->setMoreSpecialOptions(saveOptions | 262144);
				}
				solver2 = process.preProcessNonDefault(*saveSolver, 2, numberPasses, tunePreProcess);
				osiclp->getModelPtr()->setPerturbation(savePerturbation);
				osiclp->getModelPtr()->setMoreSpecialOptions(saveOptions);
			}
			// Tell solver we are not in Branch and Cut
			saveSolver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
			if (solver2)
				solver2->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
		}
		if (!solver2)
		{
			delete saveSolver;
			saveSolver = NULL;
			model_.setProblemStatus(0);
			model_.setSecondaryStatus(1);
			babModel_->setProblemStatus(0);
			babModel_->setSecondaryStatus(1);
			printf("Pre-processing says infeasible or unbounded");
			return -9999;
		}
		else
		{
			model_.setProblemStatus(-1);
			babModel_->setProblemStatus(-1);
		}
		if (model_.bestSolution())
		{
			// need to redo - in case no better found in BAB
			// just get integer part right
			const int* originalColumns = process.originalColumns();
			int numberColumns = solver2->getNumCols();
			double* bestSolution = babModel_->bestSolution();
			const double* oldBestSolution = model_.bestSolution();
			for (int i = 0; i < numberColumns; i++)
			{
				int jColumn = originalColumns[i];
				bestSolution[i] = oldBestSolution[jColumn];
			}
		}
		// 预处理后新的整数
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
		babModel_->assignSolver(solver2);
		babModel_->setOriginalColumns(process.originalColumns(), COIN_INT_MAX);
		babModel_->initialSolve();
	}
	printf("Total time is %.3f", CoinCpuTime() - start);
	// now tighten bounds
	{
		OsiClpSolverInterface* si = dynamic_cast<OsiClpSolverInterface*>(babModel_->solver());
		// get clp itself
		ClpSimplex* modelC = si->getModelPtr();
		if (modelC->tightenPrimalBounds() != 0)
		{
			printf("Problem is infeasible!");
			model_.setProblemStatus(0);
			model_.setSecondaryStatus(1);
			delete saveSolver;
			saveSolver = NULL;
			// and in babModel_ if exists
			if (babModel_) {
				babModel_->setProblemStatus(0);
				babModel_->setSecondaryStatus(1);
			}
			return -999;
		}
		si->resolve();
	}

	/*****************启发式********************/
	{
		// Feasibility Pump 
		CbcHeuristicFPump heuristicFPump(*babModel_);
		heuristicFPump.setFractionSmall(0.5);
		heuristicFPump.setMaximumPasses(30);
		/*
			>=10000000 for using obj
			>=1000000 use as accumulate switch
			>=1000 use index+1 as number of large loops
			>=100 use dextra1 as cutoff
			%100 == 10,20 etc for experimentation
			1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds
			4 and static continuous, 5 as 3 but no internal integers
			6 as 3 but all slack basis!
		*/
		// double value = babModel_->solver()->getObjSense() * babModel_->solver()->getObjValue();
		int pumpTune = 1005043;
		int w = pumpTune / 10;
		int i = w % 10;
		w /= 10;
		w /= 10;
		int r = w;
		int accumulate = r / 1000;
		r -= 1000 * accumulate;
		heuristicFPump.setAccumulate(accumulate);
		heuristicFPump.setMaximumRetries(r + 1);
		heuristicFPump.setFeasibilityPumpOptions(i * 10);
		pumpTune = pumpTune % 100;
		if (pumpTune == 6)
			pumpTune = 13;
		heuristicFPump.setWhen((pumpTune % 10) + 10);
		heuristicFPump.setHeuristicName("feasibility pump");
		int whereFromFPump = heuristicFPump.whereFrom();  // 60909
		int whenFPump = heuristicFPump.when(); // 13
		babModel_->addHeuristic(&heuristicFPump);
		// Rounding
		CbcRounding heuristicRounding(*babModel_);
		heuristicRounding.setHeuristicName("rounding");
		int whereFromRounding = heuristicRounding.whereFrom(); // 60909
		int whenRounding = heuristicRounding.when(); // 2
		babModel_->addHeuristic(&heuristicRounding);
		// Greedy
		CbcHeuristicGreedyCover heuristicGreedyCover(*babModel_);
		heuristicGreedyCover.setHeuristicName("greedy cover");
		CbcHeuristicGreedyEquality heuristicGreedyEquality(*babModel_);
		heuristicGreedyEquality.setHeuristicName("greedy equality");
		int whereFromGreedyCover = heuristicGreedyCover.whereFrom(); // 1
		int whereFromGreedyEquality = heuristicGreedyEquality.whereFrom(); // 1
		int whenGreedyCover = heuristicGreedyCover.when(); // 2
		int whenGreedyEquality = heuristicGreedyEquality.when(); // 2
		babModel_->addHeuristic(&heuristicGreedyCover);
		babModel_->addHeuristic(&heuristicGreedyEquality);
		// DiveCoefficient
		CbcHeuristicDiveCoefficient heuristicDC(*babModel_);
		heuristicDC.setHeuristicName("DiveCoefficient");
		heuristicDC.setWhen(2);
		heuristicDC.setMaxIterations(100);
		int whereFromDC = heuristicDC.whereFrom(); // 4605
		int whenDC = heuristicDC.when(); // 2
		babModel_->addHeuristic(&heuristicDC);
		//// RINS
		CbcHeuristicRINS heuristicRINS(*babModel_);
		heuristicRINS.setHeuristicName("RINS");
		heuristicRINS.setFractionSmall(0.5);
		heuristicRINS.setDecayFactor(5.0);
		int whereFromRINS = heuristicRINS.whereFrom(); // 65289
		int whenRINS = heuristicRINS.when(); // 2
		babModel_->addHeuristic(&heuristicRINS);
	}
	/*******************割平面*****************/
	int switches[30];
	int accuracyFlag[30];
	char doAtEnd[30];
	memset(doAtEnd, 0, 30);
	int numberGenerators = 0;
	int numberColumns = babModel_->solver()->getNumCols();
	// Probing
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
	probingGen.setMaxProbeRoot(CoinMin(2000, numberColumns));
	probingGen.setMaxProbeRoot(123);
	probingGen.setMaxProbe(123);
	probingGen.setMaxLookRoot(20);
	babModel_->addCutGenerator(&probingGen, -98, "Probing");
	accuracyFlag[numberGenerators] = 5;
	switches[numberGenerators++] = 0;
	// Gomory
	CglGomory gomoryGen;
	// try larger limit
	gomoryGen.setLimitAtRoot(1000);
	gomoryGen.setLimit(50);
	gomoryGen.setAwayAtRoot(0.005);
	if (numberColumns > 5000)
		gomoryGen.setLimitAtRoot(2000);
	babModel_->addCutGenerator(&gomoryGen, -98, "Gomory");
	accuracyFlag[numberGenerators] = 3;
	switches[numberGenerators++] = 0;
	// Knapsack
	CglKnapsackCover knapsackGen;
	babModel_->addCutGenerator(&knapsackGen, -98, "Knapsack");
	accuracyFlag[numberGenerators] = 1;
	switches[numberGenerators++] = -2;
	// Clique
	CglFakeClique cliqueGen(NULL, false);
	//CglClique cliqueGen(false,true);
	cliqueGen.setStarCliqueReport(false);
	cliqueGen.setRowCliqueReport(false);
	cliqueGen.setMinViolation(0.1);
	babModel_->addCutGenerator(&cliqueGen, -98, "Clique");
	accuracyFlag[numberGenerators] = 0;
	switches[numberGenerators++] = 0;
	// MixedIntegerRounding2
	CglMixedIntegerRounding2 mixedGen(1, true, 1);
	mixedGen.setDoPreproc(1); // safer (and better)
	babModel_->addCutGenerator(&mixedGen, -98, "MixedIntegerRounding2");
	accuracyFlag[numberGenerators] = 2;
	switches[numberGenerators++] = 0;
	// FlowCover
	CglFlowCover flowGen;
	babModel_->addCutGenerator(&flowGen, -98, "FlowCover");
	accuracyFlag[numberGenerators] = 2;
	switches[numberGenerators++] = 0;
	// TwoMirCuts
	CglTwomir twomirGen;
	twomirGen.setMaxElements(250);
	twomirGen.setAwayAtRoot(0.005);
	twomirGen.setAway(0.01);
	babModel_->addCutGenerator(&twomirGen, -98, "TwoMirCuts");
	accuracyFlag[numberGenerators] = 4;
	switches[numberGenerators++] = 1;
	// Say we want timings
	numberGenerators = babModel_->numberCutGenerators();
	int iGenerator;
	for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++)
	{
		CbcCutGenerator* generator = babModel_->cutGenerator(iGenerator);
		int howOften = generator->howOften();
		if (howOften == -98 || howOften == -99 || generator->maximumTries() > 0)
			generator->setSwitchOffIfLessThan(switches[iGenerator]);
		// Use if any at root as more likely later and fairly cheap
		//if (switches[iGenerator]==-2)
		//generator->setWhetherToUse(true);
		generator->setInaccuracy(accuracyFlag[iGenerator]);
		if (doAtEnd[iGenerator])
		{
			generator->setWhetherCallAtEnd(true);
			//generator->setMustCallAgain(true);
		}
		generator->setTiming(true);
	}
	// Could tune more
	double minimumDrop = fabs(babModel_->solver()->getObjValue()) * 1.0e-5 + 1.0e-5;
	babModel_->setMinimumDrop(CoinMin(5.0e-2, minimumDrop));
	if (babModel_->getNumCols() < 500)
		babModel_->setMaximumCutPassesAtRoot(-100); // always do 100 if possible
	else if (babModel_->getNumCols() < 5000)
		babModel_->setMaximumCutPassesAtRoot(100); // use minimum drop
	else
		babModel_->setMaximumCutPassesAtRoot(20);
	babModel_->setMaximumCutPasses(4);
	babModel_->solver()->setIntParam(OsiMaxNumIterationHotStart, 100);
	OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (babModel_->solver());
	// go faster stripes
	if ((osiclp->getNumRows() < 300 && osiclp->getNumCols() < 500))
	{
		osiclp->setupForRepeatedUse(2, 1);
		ClpSimplex* lp = osiclp->getModelPtr();
		int specialOptions = lp->specialOptions();
		lp->setSpecialOptions(specialOptions | (2048 + 4096));
	}
	else
		osiclp->setupForRepeatedUse(0, 1);
	double increment = babModel_->getCutoffIncrement();
	babModel_->setCutoffIncrement(CoinMax(babModel_->getCutoffIncrement(), increment));
	int mipOptions = 1057;
	osiclp->setSpecialOptions(mipOptions);
	babModel_->setSpecialOptions(babModel_->specialOptions() | 2);
	babModel_->setSpecialOptions(babModel_->specialOptions() | 512);
	babModel_->setWhenCuts(999998);
	{
		if (!babModel_->numberStrong() && babModel_->numberBeforeTrust() > 0)
			babModel_->setNumberBeforeTrust(0);
		{
			int numberColumns = babModel_->getNumCols();
			if (numberColumns <= 50)
				babModel_->setNumberBeforeTrust(1000);
			else if (numberColumns <= 100)
				babModel_->setNumberBeforeTrust(100);
			else if (numberColumns <= 300)
				babModel_->setNumberBeforeTrust(50);
		}
		osiclp = dynamic_cast<OsiClpSolverInterface*> (babModel_->solver());
		lpSolver = osiclp->getModelPtr();
		osiclp->setIntParam(OsiMaxNumIterationHotStart, 100);
		if (babModel_->solver()->getNumCols() +
			babModel_->solver()->getNumRows() < 500)
			babModel_->setFastNodeDepth(-12);
		int denseCode = -1;
		int smallCode = -1;
		if (denseCode < 0)
			denseCode = 40;
		if (smallCode < 0 && !lpSolver->factorization()->isDenseOrSmall())
			smallCode = 40;
		if (denseCode > 0)
		{
			lpSolver->factorization()->setGoDenseThreshold(denseCode);
			assert(osiclp == babModel_->solver());
			osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
		}
		if (smallCode > 0 && smallCode > denseCode)
			lpSolver->factorization()->setGoSmallThreshold(smallCode);
		if (lpSolver->factorization()->goOslThreshold() > 1000)
		{
			// use osl in gomory (may not if CglGomory decides not to)
			int numberGenerators = babModel_->numberCutGenerators();
			int nGomory = 0;
			for (int iGenerator = 0; iGenerator < numberGenerators; iGenerator++)
			{
				CbcCutGenerator* generator = babModel_->cutGenerator(iGenerator);
				CglGomory* gomory = dynamic_cast<CglGomory*>
					(generator->generator());
				if (gomory)
				{
					if (nGomory < 2)
						gomory->useAlternativeFactorization();
					else if (gomory->originalSolver())
					{
						OsiClpSolverInterface* clpSolver = dynamic_cast<OsiClpSolverInterface*>(gomory->originalSolver());
						if (clpSolver)
						{
							ClpSimplex* simplex = clpSolver->getModelPtr();
							simplex->factorization()->setGoOslThreshold(0);
						}
					}
					nGomory++;
				}
			}
		}
		babModel_->solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
		babModel_->setStrongStrategy(0);
		babModel_->setMultipleRootTries(0);
		babModel_->setMoreSpecialOptions2(0);

		babModel_->setAllowablePercentageGap(0.0004);

		// 分支定界
		babModel_->branchAndBound(0);
		int truncateColumns = COIN_INT_MAX;
		int truncateRows = -1;
		if (truncateColumns < babModel_->solver()->getNumCols())
		{
			OsiSolverInterface* solverX = babModel_->solver();
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
		}
		int numberSolutions = babModel_->numberSavedSolutions();
		if (numberSolutions > 1)
		{
			for (int iSolution = numberSolutions - 1; iSolution >= 0; iSolution--)
			{
				model_.setBestSolution(babModel_->savedSolution(iSolution),
					model_.solver()->getNumCols(),
					babModel_->savedSolutionObjective(iSolution));
			}
		}
	}
	//printf("result is %.2f\n", babModel_->getObjValue());
	//printf("time is %.2f", CoinCpuTime() - timeStart);
	/*****************结果*******************/
	osiclp = dynamic_cast<OsiClpSolverInterface*> (babModel_->solver());
	printf("Cuts at root node changed objective from %g to %g\n",
		babModel_->getContinuousObjective(), babModel_->rootObjectiveAfterCuts());
	numberGenerators = babModel_->numberCutGenerators();
	int* statistics_number_cuts = new int[numberGenerators];;
	int statistics_number_generators = numberGenerators;
	double statistics_cut_time = 0.0;
	const char** statistics_name_generators = new const char* [numberGenerators];
	char timing[30];
	for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
		CbcCutGenerator* generator = babModel_->cutGenerator(iGenerator);
		statistics_name_generators[iGenerator] = generator->cutGeneratorName();
		statistics_number_cuts[iGenerator] = generator->numberCutsInTotal();
		printf("%s was tried %d times and created %d cuts of which %d were active after adding rounds of cuts",
			generator->cutGeneratorName(),
			generator->numberTimesEntered(),
			generator->numberCutsInTotal() +
			generator->numberColumnCuts(),
			generator->numberCutsActive());
		if (generator->timing())
		{
			printf(" (%.3f seconds)\n", generator->timeInCutGenerator());
			statistics_cut_time += generator->timeInCutGenerator();
		}
		CglStored* stored = dynamic_cast<CglStored*>(generator->generator());
		if (stored && !generator->numberCutsInTotal())
			continue;
		CglImplication* implication = dynamic_cast<CglImplication*>(generator->generator());
		if (implication && !generator->numberCutsInTotal())
			continue;
	}
	double time2 = CoinCpuTime() + CoinCpuTimeJustChildren();
	totalTime += time2 - timeStart;
	// For best solution
	double* bestSolution = NULL;
	if (babModel_->getMinimizationObjValue() < 1.0e50)
	{
		// post process
		int n;
		if (preProcess)
		{
			n = saveSolver->getNumCols();
			bestSolution = new double[n];
			OsiClpSolverInterface* clpSolver = dynamic_cast<OsiClpSolverInterface*> (babModel_->solver());
			// Save bounds on processed model
			const int* originalColumns = process.originalColumns();
			int numberColumns2 = clpSolver->getNumCols();
			double* solution2 = new double[n];
			double* lower2 = new double[n];
			double* upper2 = new double[n];
			for (int i = 0; i < n; i++)
			{
				solution2[i] = COIN_DBL_MAX;
				lower2[i] = COIN_DBL_MAX;
				upper2[i] = -COIN_DBL_MAX;
			}
			const double* columnLower = clpSolver->getColLower();
			const double* columnUpper = clpSolver->getColUpper();
			const double* solution = babModel_->bestSolution();
			for (int i = 0; i < numberColumns2; i++)
			{
				int jColumn = originalColumns[i];
				if (jColumn < n)
				{
					solution2[jColumn] = solution[i];
					lower2[jColumn] = columnLower[i];
					upper2[jColumn] = columnUpper[i];
				}
			}
			ClpSimplex* lpSolver = clpSolver->getModelPtr();
			lpSolver->setSpecialOptions(lpSolver->specialOptions() | IN_BRANCH_AND_BOUND); // say is Cbc (and in branch and bound)
			// put back any saved solutions
			putBackOtherSolutions(babModel_, &model_, &process);
			process.postProcess(*babModel_->solver());
			bool tightenB = false;
			{
				int n = babModel_->numberObjects();
				for (int i = 0; i < n; i++)
				{
					const OsiObject* obj = babModel_->object(i);
					if (!dynamic_cast<const CbcSimpleInteger*>(obj))
					{
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
			for (int i = 0; i < n; i++)
			{
				if (!saveSolver->isInteger(i) && !tightenB)
					continue;
				if (lower2[i] != COIN_DBL_MAX)
				{
					if (lower2[i] != columnLower[i] ||
						upper2[i] != columnUpper[i])
					{
						numberChanged++;
						saveSolver->setColLower(i, lower2[i]);
						saveSolver->setColUpper(i, upper2[i]);
					}
				}
			}
			delete[] solution2;
			delete[] lower2;
			delete[] upper2;
			if (numberChanged)
				printf("%d bounds tightened after postprocessing\n", numberChanged);
			saveSolver->resolve();
			if (!saveSolver->isProvenOptimal())
			{
				// try all slack
				CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel_->solver()->getEmptyWarmStart());
				saveSolver->setWarmStart(basis);
				delete basis;
				saveSolver->initialSolve();
				OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (saveSolver);
				if (osiclp)
					osiclp->getModelPtr()->checkUnscaledSolution();
			}
			// and original solver
			originalSolver->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
			n = originalSolver->getNumCols();
			originalSolver->setColLower(saveSolver->getColLower());
			originalSolver->setColUpper(saveSolver->getColUpper());
			// basis
			CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel_->solver()->getWarmStart());
			originalSolver->setBasis(*basis);
			delete basis;
			originalSolver->resolve();
			if (!originalSolver->isProvenOptimal())
			{
				// try all slack
				CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel_->solver()->getEmptyWarmStart());
				originalSolver->setBasis(*basis);
				delete basis;
				originalSolver->initialSolve();
				OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (originalSolver);
				if (osiclp)
					osiclp->getModelPtr()->checkUnscaledSolution();
			}
			babModel_->assignSolver(saveSolver);
			memcpy(bestSolution, babModel_->solver()->getColSolution(), n * sizeof(double));
		}
		babModel_->deleteSolutions();
		babModel_->setBestSolution(bestSolution, n, babModel_->getMinimizationObjValue());
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
			CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel_->solver()->getWarmStart());
			originalSolver->setBasis(*basis);
			delete basis;
			originalSolver->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
			originalSolver->resolve();
			if (!originalSolver->isProvenOptimal())
			{
				// try all slack
				CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel_->solver()->getEmptyWarmStart());
				originalSolver->setBasis(*basis);
				delete basis;
				originalSolver->initialSolve();
				OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*> (originalSolver);
				if (osiclp)
					osiclp->getModelPtr()->checkUnscaledSolution();
			}

		}
	}
	else if (model_.bestSolution() && model_.getMinimizationObjValue() < 1.0e50 && preProcess)
	{
		printf("Restoring heuristic best solution of %g", model_.getMinimizationObjValue());
		int n = saveSolver->getNumCols();
		bestSolution = new double[n];
		// Put solution now back in saveSolver
		saveSolver->setColSolution(model_.bestSolution());
		babModel_->assignSolver(saveSolver);
		saveSolver = NULL;
		babModel_->setMinimizationObjValue(model_.getMinimizationObjValue());
		memcpy(bestSolution, babModel_->solver()->getColSolution(), n * sizeof(double));
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
			CoinWarmStartBasis* basis = dynamic_cast<CoinWarmStartBasis*> (babModel_->solver()->getWarmStart());
			originalSolver->setBasis(*basis);
			delete basis;
		}
	}
	lpSolver = clpSolver->getModelPtr();
	{
		//move best solution (should be there -- but ..)
		int n = lpSolver->getNumCols();
		if (bestSolution) {
			memcpy(lpSolver->primalColumnSolution(), bestSolution, n * sizeof(double));
			// now see what that does to row solution
			int numberRows = lpSolver->numberRows();
			double* rowSolution = lpSolver->primalRowSolution();
			memset(rowSolution, 0, numberRows * sizeof(double));
			lpSolver->clpMatrix()->times(1.0, bestSolution, rowSolution);
			lpSolver->setObjectiveValue(babModel_->getObjValue());
		}
	}
	delete saveSolver;
	delete[] bestSolution;
	std::string statusName[] = { "", "Stopped on ", "Run abandoned", "", "", "User ctrl-c" };
	std::string minor[] = { "Optimal solution found", "Linear relaxation infeasible", "Optimal solution found (within gap tolerance)", "node limit", "time limit", "user ctrl-c", "solution limit", "Linear relaxation unbounded", "Problem proven infeasible" };
	int iStat = babModel_->status();
	int iStat2 = babModel_->secondaryStatus();
	if (!iStat && !iStat2 && !bestSolution)
		iStat2 = 8;
	if (!iStat && iStat2 == 1 && bestSolution)
		iStat2 = 0; // solution and search completed
	double statistics_seconds = time2 - timeStart;
	double statistics_sys_seconds = CoinSysTime();
	double statistics_elapsed_seconds = CoinWallclockTime();
	double statistics_obj = babModel_->getObjValue();
	double statistics_continuous = babModel_->getContinuousObjective();
	double statistics_tighter = babModel_->rootObjectiveAfterCuts();
	double statistics_nodes = babModel_->getNodeCount();
	double statistics_iterations = babModel_->getIterationCount();
	std::string tatistics_result = statusName[iStat];
	printf("\nResult - %s%s\n", statusName[iStat].c_str(), minor[iStat2].c_str());
	if (babModel_->bestSolution())
		printf("Objective value:                %.8f\n", babModel_->getObjValue());
	else
		printf("No feasible solution found\n");
	if (iStat2 >= 2 && iStat2 <= 6)
	{
		printf("Lower bound:                    %.3f\n",
			babModel_->getBestPossibleObjValue());
		if (babModel_->bestSolution())
		{
			printf("Gap:                            %.2f\n",
				(babModel_->getObjValue() - babModel_->getBestPossibleObjValue()) /
				fabs(babModel_->getBestPossibleObjValue()));
		}
	}
	printf("Enumerated nodes:               %d\n", babModel_->getNodeCount());
	printf("Total iterations:               %d\n", babModel_->getIterationCount());
	printf("Time (CPU seconds):             %.2f\n", CoinCpuTime() - timeStart);
}

