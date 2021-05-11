#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <algorithm>

#include "ShiftAndPropagate.hpp"
#include "CbcModel.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
ShiftAndPropagate::ShiftAndPropagate()
	:CbcHeuristic()
{

}

// Constructor with model
ShiftAndPropagate::ShiftAndPropagate(CbcModel& model)
	:CbcHeuristic(model)
{

}

// Copy constructor
ShiftAndPropagate::ShiftAndPropagate(const ShiftAndPropagate& rhs)
	:CbcHeuristic(rhs)
{

}

// Destructor
ShiftAndPropagate::~ShiftAndPropagate()
{

}

// Clone
CbcHeuristic*
ShiftAndPropagate::clone() const
{
	return new ShiftAndPropagate(*this);
}

// Assignment operator
ShiftAndPropagate&
ShiftAndPropagate::operator=(const ShiftAndPropagate& rhs)
{
	if (this != &rhs)
		CbcHeuristic::operator=(rhs);
	return *this;
}

void 
ShiftAndPropagate::resetModel(CbcModel* model)
{ }

void 
ShiftAndPropagate::setModel(CbcModel* model)
{ }

// transform integer variable
void
ShiftAndPropagate::transformIntVar(CoinPackedMatrix* matrixByRow, double* colUpper, double* colLower, double* rowUpper, double* rowLower, int iCol)
{
	double shiftValue;
	bool negateCoeff;

	if (fabs(colLower[iCol]) < fabs(colUpper[iCol]))
	{
		shiftValue = colLower[iCol];
		negateCoeff = false;
		if (colUpper[iCol] != COIN_DBL_MAX)
			colUpper[iCol] -= colLower[iCol];
	}
	else
	{
		shiftValue = colUpper[iCol];
		negateCoeff = true;
		if (colLower[iCol] != -COIN_DBL_MAX)
			colUpper[iCol] = -colLower[iCol] + colUpper[iCol];
		else
			colUpper[iCol] = COIN_DBL_MAX;

	}
	colLower[iCol] = 0;

	int rowNum = matrixByRow->getNumRows();
	for (int i = 0; i < rowNum; i++)
	{
		double aij = matrixByRow->getCoefficient(i, iCol);
		if (aij == 0)
			continue;
		if (rowLower[i] != -COIN_DBL_MAX)
			rowLower[i] -= aij * shiftValue;
		if (rowUpper[i] != COIN_DBL_MAX)
			rowUpper[i] -= aij * shiftValue;
		if (negateCoeff)
			matrixByRow->modifyCoefficient(i, iCol, -aij);
	}
}

// relax continuous variable
void
ShiftAndPropagate::relaxContVar(CoinPackedMatrix* matrix, const double* colUpper, const double* colLower, double* rowUpper, double* rowLower, int iCol)
{
	int rowNum = matrix->getNumRows();
	for (int i = 0; i < rowNum; i++)
	{
		double aij = matrix->getCoefficient(i, iCol);
		if (aij == 0)
			continue;
		else if (aij < 0)
		{
			if (rowUpper[i] != COIN_DBL_MAX)
			{
				rowUpper[i] -= aij * colUpper[iCol];
				if (isinf(rowUpper[i]))
					rowUpper[i] = COIN_DBL_MAX;
			}
			if (rowLower[i] != -COIN_DBL_MAX)
			{
				rowLower[i] -= aij * colLower[iCol];
				if (isinf(rowLower[i]))
					rowLower[i] = -COIN_DBL_MAX;
			}
		}
		else
		{
			if (rowUpper[i] != COIN_DBL_MAX)
			{
				rowUpper[i] -= aij * colLower[iCol];
				if (isinf(rowUpper[i]))
					rowUpper[i] = COIN_DBL_MAX;
			}
			if (rowLower[i] != -COIN_DBL_MAX)
			{
				rowLower[i] -= aij * colUpper[iCol];
				if (isinf(rowLower[i]))
					rowLower[i] = -COIN_DBL_MAX;
			}
		}
		matrix->modifyCoefficient(i, iCol, 0);
	}
}

// do problem transformation
OsiSolverInterface*
ShiftAndPropagate::transformProblem(CbcModel* model_)
{
	int i, j;
	int rowNum = model_->getNumRows();
	int colNum = model_->getNumCols();
	const CoinPackedMatrix* matrix = model_->getMatrixByRow();
	const double* rowLower = model_->getRowLower();
	const double* rowUpper = model_->getRowUpper();
	double* transedRowLower = new double[rowNum]();
	memcpy(transedRowLower, rowLower, sizeof(double) * rowNum);
	double* transedRowUpper = new double[rowNum]();
	memcpy(transedRowUpper, rowUpper, sizeof(double) * rowNum);
	CoinPackedMatrix transedMatrix(*matrix);

	/*modify the lhs and rhs, i.e divide the lhs and rhs by the maximum absolute value of the row*/
	for (i = 0; i < rowNum; i++)
	{
		CoinPackedVector ithVec = matrix->getVector(i); 

		double maxCoeff = ithVec.infNorm();
		if (transedRowLower[i] != -COIN_DBL_MAX)
			transedRowLower[i] /= maxCoeff;
		if (transedRowUpper[i] != COIN_DBL_MAX)
			transedRowUpper[i] /= maxCoeff;

		const int* indices = ithVec.getIndices();
		const double* elements = ithVec.getElements();
		int eleNum = ithVec.getNumElements();

		for (j = 0; j < eleNum; j++)
		{
			double newElement = elements[j] / maxCoeff;
			int iCol = indices[j];
			transedMatrix.modifyCoefficient(i, iCol, newElement);
		}
	}

	const double* colLower = model_->getColLower();
	const double* colUpper = model_->getColUpper();
	double* transedColLower = new double[colNum]();
	memcpy(transedColLower, colLower, sizeof(double) * colNum);
	double* transedColUpper = new double[colNum]();
	memcpy(transedColUpper, colUpper, sizeof(double) * colNum);
	for (j = 0; j < colNum; j++)
	{
		if (model_->isInteger(j))
		{
			if (colLower[j] == -COIN_DBL_MAX && colUpper[j] == COIN_DBL_MAX)  // 不考虑自由整数变量
				return NULL;
			transformIntVar(&transedMatrix, transedColUpper, transedColLower, transedRowUpper, transedRowLower, j);
		}
		else
		{
			relaxContVar(&transedMatrix, transedColUpper, transedColLower, transedRowUpper, transedRowLower, j);
		}
	}

	// const int* integerIndices = model_->integerVariable();
	int integerNum = model_->numberIntegers();
	int continuousNum = colNum - integerNum;
	const double* coefficients = model_->getObjCoefficients();
	double* transedCoeffs = new double[integerNum];
	int* continousIndices = new int[continuousNum];
	double* finalColLower = new double[integerNum];
	double* finalColUpper = new double[integerNum];
	int coeffCount = 0;
	int continousCount = 0;
	int colBoundCount = 0;
	for (j = 0; j < colNum; j++)
	{
		if (model_->isInteger(j))
		{
			transedCoeffs[coeffCount++] = coefficients[j];
			finalColLower[colBoundCount] = transedColLower[j];
			finalColUpper[colBoundCount] = transedColUpper[j];
			colBoundCount++;
		}
		else
			continousIndices[continousCount++] = j;
	}
	transedMatrix.deleteCols(continousCount, continousIndices);

	OsiSolverInterface* returnSolver = model_->solver()->clone();

	returnSolver->loadProblem(transedMatrix, finalColLower, finalColUpper, transedCoeffs, transedRowLower, transedRowUpper);

	for (int i = 0; i < returnSolver->getNumCols(); i++)
		returnSolver->setInteger(i);
	delete[] transedColLower;
	delete[] transedColUpper;
	delete[] transedRowLower;
	delete[] transedRowUpper;
	delete[] transedCoeffs;
	delete[] continousIndices;
	delete[] finalColLower;
	delete[] finalColUpper;

	return returnSolver;
}

/*compare function for variables' importance*/
bool 
ShiftAndPropagate::cmpForVarImportance(std::pair<int, double>a, std::pair<int, double>b)
{
	return a.second > b.second;
}

int 
ShiftAndPropagate::decideBestShift(OsiSolverInterface* transProblem, int colIndex)
{
	const CoinPackedMatrix* matrix = transProblem->getMatrixByRow();
	const double* rowLower = transProblem->getRowLower();
	const double* rowUpper = transProblem->getRowUpper();
	const double* colLower = transProblem->getColLower();
	const double* colUpper = transProblem->getColUpper();
	int rowNumber = transProblem->getNumRows();
	int t;
	std::set<std::pair<int, int>>Q;
	for (int i = 0; i < rowNumber; i++)
	{
		double aij = matrix->getCoefficient(i, colIndex);
		if (aij == 0) continue;

		if (rowUpper[i] != COIN_DBL_MAX)
		{
			if (rowUpper[i] < 0 && aij < 0)
			{
				t = std::ceil(rowUpper[i] / aij);
				if (colLower[colIndex] <= t && t <= colUpper[colIndex])
					Q.insert(std::pair<int, int>(t, -1));
			}
			else if (rowUpper[i] >= 0 && aij > 0)
			{
				t = std::ceil(rowUpper[i] / aij + 1e-5);   // 必须要使得rowUpper[i]严格小于0才行，否则若rowUpper[i] / aij刚好为整数，那么此时不加1e-5就会使得rowUpper[i]=0
				if (colLower[colIndex] <= t && t <= colUpper[colIndex])
					Q.insert(std::pair<int, int>(t, 1));
			}
		}

		if (rowLower[i] != -COIN_DBL_MAX)
		{
			if (rowLower[i] > 0 && aij > 0)
			{
				t = std::ceil(rowLower[i] / aij);
				if (colLower[colIndex] <= t && t <= colUpper[colIndex])
					Q.insert(std::pair<int, int>(t, -1));
			}
			else if (rowLower[i] <= 0 && aij < 0)
			{
				t = std::ceil(rowLower[i] / aij + 1e-5);
				if (colLower[colIndex] <= t && t <= colUpper[colIndex])
					Q.insert(std::pair<int, int>(t, 1));
			}
		}
	}

	if (Q.empty())
		return 0;

	bool exsitZero = false;
	bool exsitOne = false;
	for (auto it = Q.begin(); it != Q.end(); it++)
	{
		if ((*it).first == 0)
			exsitZero = true;
		else
			exsitOne = true;
	}
	if (!exsitOne)
		return 0;
	if (!exsitZero)
		return 1;

	int numberOne = 0;
	int numberZero = 0;
	for (auto it = Q.begin(); it != Q.end(); it++)
	{
		if ((*it).first == 0)
			numberZero += (*it).second;
		else
			numberOne += (*it).second;
	}
	if (numberZero <= numberZero)
		return 1;
	else
		return 0;
}

// MAIN PROCEDURE 
int
ShiftAndPropagate::solution(double& solutionValue, double* betterSolution)
{
	int integerNumber = model_->numberIntegers();
	if (integerNumber == 0)
		return 0;  // there are no integer variables, return

	int depth = model_->currentDepth();
	if (depth > 0)
		return 0; // we only do shift and propagate at root node

	int feasSolNumber = model_->getSolutionCount();
	if (feasSolNumber != 0)
		return 0;  // if there is already a primarily feasible solution, return

	int rowNumber = model_->getNumRows();
	int columnNumber = model_->getNumCols();
	if (columnNumber == 0)
		return 0; // it is not a valid model, return

	const double* originalColLower = model_->getColLower();
	const double* originalColUpper = model_->getColUpper();
	const CoinPackedMatrix* originalMatrix = model_->getMatrixByRow();
	const int* originalIntegerIndices = model_->integerVariable();

	// transform problem
	OsiSolverInterface* transedProblem = transformProblem(model_);

	const CoinPackedMatrix* transedMatrixByCol = transedProblem->getMatrixByCol(); 
	int originalColNumber = model_->getNumCols();
	int transedRowNumber = transedMatrixByCol->getNumRows();
	int transedColNumber = transedMatrixByCol->getNumCols();

	const double* transedColLower = transedProblem->getColLower();
	const double* transedColUpper = transedProblem->getColUpper();
	const double* transedRowLower = transedProblem->getRowLower();
	const double* transedRowUpper = transedProblem->getRowUpper();

	std::vector<std::pair<int, double>>varImportance;
	for (int j = 0; j < integerNumber; j++)
	{
		CoinPackedVector ithColVector = transedMatrixByCol->getVector(j);
		int nonZeroCoefNum = 0;
		double absCoefSum = 0.0;
		for (int i = 0; i < transedRowNumber; i++)
		{
			if (ithColVector[i] != 0)
			{
				nonZeroCoefNum++;
				absCoefSum += std::fabs(ithColVector[i]);
			}
		}
		double importance = nonZeroCoefNum + absCoefSum;
		varImportance.push_back(std::pair<int, double>(j, importance));
	}
	//std::sort(varImportance.begin(), varImportance.end(), cmpForVarImportance);

	/*start shift and propagate*/
	double* constant = new double[transedRowNumber]();    // 用于在每固定一个变量之后，计算对于每一行所产生的的常量值
	double* bestshiftForEveryVar = new double[transedColNumber]();  // 用于记录每个变量的bestshift值
	bool isSatisfy = true;
	OsiSolverInterface* solver = model_->solver()->clone();
	for (auto it = varImportance.begin(); it != varImportance.end(); it++)
	{
		int iCol = (*it).first;
		if (transedColLower[iCol] == transedColUpper[iCol])  // 如果相等，肯定都为0，因为问题转换已经将下界变为0，所以不需要更改constant的值
		{
			solver->setColLower(iCol, transedColLower[iCol]);
			solver->setColUpper(iCol, transedColLower[iCol]);
			continue;
		}
		for (int i = 0; i < transedRowNumber; i++)  // 判断0解是否满足
		{
			if (!(transedRowLower[i] - constant[i] <= 0 && transedRowUpper[i] - constant[i] >= 0))
			{
				isSatisfy = false;
				break;
			}
			isSatisfy = true;
		}
		if (isSatisfy) break;  // 0解满足，则停止

		int bestshift = decideBestShift(transedProblem, iCol);

		const double* currentColLower = transedProblem->getColLower();
		const double* currentColUpper = transedProblem->getColUpper();
		double* backtrackColLower = new double[transedColNumber]();
		double* backtrackColUpper = new double[transedColNumber]();
		memcpy(backtrackColLower, currentColLower, sizeof(double) * transedColNumber);
		memcpy(backtrackColUpper, currentColUpper, sizeof(double) * transedColNumber);

		transedProblem->setColLower(iCol, bestshift);
		transedProblem->setColUpper(iCol, bestshift);
		CglPreProcess process;
		int infeasibility = process.tightenPrimalBounds(*transedProblem);
		if (infeasibility)
		{
			transedProblem->setColLower(backtrackColLower);
			transedProblem->setColUpper(backtrackColUpper);
			//showMatrix(transedProblem);
		}
		else
		{
			transedProblem->setColLower(backtrackColLower);
			transedProblem->setColUpper(backtrackColUpper);
			transedProblem->setColLower(iCol, bestshift);
			transedProblem->setColUpper(iCol, bestshift);
			bestshiftForEveryVar[iCol] = bestshift;
			//showMatrix(transedProblem);

			if (bestshift != 0)
			{
				for (int i = 0; i < transedRowNumber; i++)
				{
					double temp = (transedMatrixByCol->getCoefficient(i, iCol) * bestshift);
					constant[i] += temp;
				}
			}
		}
		delete[] backtrackColLower;
		delete[] backtrackColUpper;
	}
	int returnCode = isSatisfy;
	if (isSatisfy)
	{
		// 转换回对应的原问题中的值
		for (int j = 0; j < integerNumber; j++)
		{
			int iCol = originalIntegerIndices[j];
			if (std::fabs(originalColLower[iCol]) <= std::fabs(originalColUpper[iCol]))
			{
				double originalValue = bestshiftForEveryVar[j] + originalColLower[iCol];
				solver->setColLower(iCol, originalValue);
				solver->setColUpper(iCol, originalValue);
			}
			else
			{
				double originalValue = originalColUpper[iCol] - bestshiftForEveryVar[j];
				solver->setColLower(iCol, originalValue);
				solver->setColUpper(iCol, originalValue);
			}
		}
		//showMatrix(&solver);

		solver->initialSolve();
		returnCode = solver->isProvenOptimal();
		if (solver->isProvenOptimal())
		{
			const double* solverSol = solver->getColSolution();
			memcpy(betterSolution, solverSol, sizeof(double)* columnNumber);
			double solverObjVal = solver->getObjValue();
			solutionValue = solverObjVal;
		}
	}
	return returnCode;
}