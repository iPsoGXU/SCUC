// This is a variant of RINS which is more appropriate to the problem
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "VRINS.h"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
VRINS::VRINS()
	: CbcHeuristic()
{
	numberSolutions_ = 0;
	numberSuccesses_ = 0;
	numberTries_ = 0;
	stateOfFixing_ = 0;
	shallowDepth_ = 0;
	lastNode_ = -999999;
	howOften_ = 100;
	decayFactor_ = 0.5;
	used_ = NULL;
	whereFrom_ = 1 + 8 + 16 + 255 * 256;
	whereFrom_ = 1 + 8 + 255 * 256;
}

// Constructor with model - assumed before cuts

VRINS::VRINS(CbcModel& model)
	: CbcHeuristic(model)
{
	numberSolutions_ = 0;
	numberSuccesses_ = 0;
	numberTries_ = 0;
	stateOfFixing_ = 0;
	shallowDepth_ = 0;
	lastNode_ = -999999;
	howOften_ = 100;
	decayFactor_ = 0.5;
	assert(model.solver());
	int numberColumns = model.solver()->getNumCols();
	used_ = new char[numberColumns];
	memset(used_, 0, numberColumns);
	whereFrom_ = 1 + 8 + 16 + 255 * 256;
	whereFrom_ = 1 + 8 + 255 * 256;
}

// Destructor
VRINS::~VRINS()
{
	delete[] used_;
}

// Clone
CbcHeuristic*
VRINS::clone() const
{
	return new VRINS(*this);
}

// Assignment operator
VRINS&
VRINS::operator=(const VRINS& rhs)
{
	if (this != &rhs) {
		CbcHeuristic::operator=(rhs);
		numberSolutions_ = rhs.numberSolutions_;
		howOften_ = rhs.howOften_;
		numberSuccesses_ = rhs.numberSuccesses_;
		numberTries_ = rhs.numberTries_;
		stateOfFixing_ = rhs.stateOfFixing_;
		lastNode_ = rhs.lastNode_;
		delete[] used_;
		if (model_ && rhs.used_) {
			int numberColumns = model_->solver()->getNumCols();
			used_ = new char[numberColumns];
			memcpy(used_, rhs.used_, numberColumns);
		}
		else {
			used_ = NULL;
		}
	}
	return *this;
}

// Create C++ lines to get to current state
void
VRINS::generateCpp(FILE* fp)
{
	VRINS other;
	fprintf(fp, "0#include \"CbcHeuristicVRINS.hpp\"\n");
	fprintf(fp, "3  CbcHeuristicVRINS heuristicVRINS(*cbcModel);\n");
	CbcHeuristic::generateCpp(fp, "heuristicVRINS");
	if (howOften_ != other.howOften_)
		fprintf(fp, "3  heuristicVRINS.setHowOften(%d);\n", howOften_);
	else
		fprintf(fp, "4  heuristicVRINS.setHowOften(%d);\n", howOften_);
	fprintf(fp, "3  cbcModel->addHeuristic(&heuristicVRINS);\n");
}

// Copy constructor
VRINS::VRINS(const VRINS& rhs)
	:
	CbcHeuristic(rhs),
	numberSolutions_(rhs.numberSolutions_),
	howOften_(rhs.howOften_),
	numberSuccesses_(rhs.numberSuccesses_),
	numberTries_(rhs.numberTries_),
	stateOfFixing_(rhs.stateOfFixing_),
	lastNode_(rhs.lastNode_)
{
	if (model_ && rhs.used_) {
		int numberColumns = model_->solver()->getNumCols();
		used_ = new char[numberColumns];
		memcpy(used_, rhs.used_, numberColumns);
	}
	else {
		used_ = NULL;
	}
}
// Resets stuff if model changes
void
VRINS::resetModel(CbcModel* /*model*/)
{
	//CbcHeuristic::resetModel(model);
	delete[] used_;
	stateOfFixing_ = 0;
	if (model_ && used_) {
		int numberColumns = model_->solver()->getNumCols();
		used_ = new char[numberColumns];
		memset(used_, 0, numberColumns);
	}
	else {
		used_ = NULL;
	}
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
VRINS::solution(double& solutionValue, double* betterSolution)
{
#pragma region first
	numCouldRun_++;
	int returnCode = 0;
	const double* bestSolution = model_->bestSolution(); // incumbent
	if (!bestSolution)  // no incumbent, return
		return 0;

	OsiSolverInterface* solver = model_->solver();
	const double* currentSolution = solver->getColSolution(); // continuous relaxation
	const int* integerVarIndex = model_->integerVariable();  // integer variable¡®s index
	const int numberInteger = model_->numberIntegers();   // number of integers
	OsiSolverInterface* newSolver = model_->continuousSolver()->clone();
	int numberColumns = newSolver->getNumCols();
	int nFixed = 0;  // record the number of fixed integers
	for (size_t i = 0; i < numberInteger; i++)
	{
		int iColumn = integerVarIndex[i];
		double valueOfInteger = bestSolution[iColumn];
		double valueOfContinuous = currentSolution[iColumn];
		if (fabs(valueOfInteger - valueOfContinuous) < 0.3)   // if the value of var in incumbent == the value in continous relaxation
		{
			// fix the var
			newSolver->setColLower(iColumn, valueOfInteger);
			newSolver->setColUpper(iColumn, valueOfInteger);
			nFixed++;
		}
	}
	if (5 * nFixed > numberInteger) // if more than 1/5 integers are fixed, solve sub MIP
	{
		if (solutionValue == -COIN_DBL_MAX) {
			// return fixings in betterSolution
			const double* colLower = newSolver->getColLower();
			const double* colUpper = newSolver->getColUpper();
			for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
				if (colLower[iColumn] == colUpper[iColumn])
					betterSolution[iColumn] = colLower[iColumn];
				else
					betterSolution[iColumn] = COIN_DBL_MAX;
			}
			delete newSolver;
			return 0;
		}
		returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue, model_->getCutoff(), "VRINS");
	}
	numberTries_++;
	if ((numberTries_ % 10) == 0 && numberSuccesses_ * 3 < numberTries_)
		howOften_ += static_cast<int> (howOften_ * decayFactor_);
	delete newSolver;

	return returnCode;
#pragma endregion
}
// update model
void VRINS::setModel(CbcModel* model)
{
	model_ = model;
	// Get a copy of original matrix
	assert(model_->solver());
	delete[] used_;
	int numberColumns = model->solver()->getNumCols();
	used_ = new char[numberColumns];
	memset(used_, 0, numberColumns);
}



