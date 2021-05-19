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
#include "OsiClpSolverInterface.hpp"

// Default Constructor
ShiftAndPropagate::ShiftAndPropagate()
	:CbcHeuristic()
{
	maxBacktrackLimit_ = -1;
}

// Constructor with model
ShiftAndPropagate::ShiftAndPropagate(CbcModel& model)
	:CbcHeuristic(model)
{
	maxBacktrackLimit_ = -1;
}

// Copy constructor
ShiftAndPropagate::ShiftAndPropagate(const ShiftAndPropagate& rhs)
	:CbcHeuristic(rhs)
{
	maxBacktrackLimit_ = -1;
}

// Destructor
ShiftAndPropagate::~ShiftAndPropagate()
{
	maxBacktrackLimit_ = -1;
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


// relax continuous variable
void
ShiftAndPropagate::relaxContinuousVar(CoinPackedMatrix* matrixByCol, double colUpper, double colLower, double* rowUpper, double* rowLower, int iCol)
{
	const int* indicesByCol = matrixByCol->getIndices(); // 每列非0元素索引
	const int* startsByCol = matrixByCol->getVectorStarts(); // 每列第一个非0元素为所有非0元素中的第几个
	const int* lengthsByCol = matrixByCol->getVectorLengths(); // 每列非0元素个数
	int elenum = matrixByCol->getNumElements();
	int involvedRowNum = lengthsByCol[iCol]; // 包含该连续变量的行数
	int start = startsByCol[iCol];
	int i = 0;
	while (i < involvedRowNum)
	{
		int iRow = indicesByCol[start];
		double aij = matrixByCol->getCoefficient(iRow, iCol);
		if (aij < 0)
		{
			if (rowUpper[iRow] != COIN_DBL_MAX)
			{
				rowUpper[iRow] -= aij * colUpper;
				if (isinf(rowUpper[iRow]))
					rowUpper[iRow] = COIN_DBL_MAX;
			}
			if (rowLower[iRow] != -COIN_DBL_MAX)
			{
				rowLower[iRow] -= aij * colLower;
				if (isinf(rowLower[iRow]))
					rowLower[iRow] = -COIN_DBL_MAX;
			}
		}
		else
		{
			if (rowUpper[iRow] != COIN_DBL_MAX)
			{
				rowUpper[iRow] -= aij * colLower;
				if (isinf(rowUpper[iRow]))
					rowUpper[iRow] = COIN_DBL_MAX;
			}
			if (rowLower[iRow] != -COIN_DBL_MAX)
			{
				rowLower[iRow] -= aij * colUpper;
				if (isinf(rowLower[iRow]))
					rowLower[iRow] = -COIN_DBL_MAX;
			}
		}
		matrixByCol->modifyCoefficient(iRow, iCol, 0, true);
		start++;
		i++;
	}
}

// do problem transformation
OsiSolverInterface*
ShiftAndPropagate::transformProblem(CbcModel* model_)
{
	std::cout << "start transform problem method!" << std::endl;

	int i, j;
	int originalRowNum = model_->getNumRows();
	int originalColNum = model_->getNumCols();
	const CoinPackedMatrix* originalMatrixByRow = model_->getMatrixByRow();
	const CoinPackedMatrix* originalMatrixByCol = model_->getMatrixByCol();
	CoinPackedMatrix transedMatrixByCol(*originalMatrixByCol);

	const double* originalRowLower = model_->getRowLower();
	const double* originalRowUpper = model_->getRowUpper();
	double* transedRowLower = new double[originalRowNum];
	memcpy(transedRowLower, originalRowLower, sizeof(double) * originalRowNum);
	double* transedRowUpper = new double[originalRowNum];
	memcpy(transedRowUpper, originalRowUpper, sizeof(double) * originalRowNum);

	/*	矩阵系数标准化：
		1. 矩阵所有系数变为1
		2. 修改对应的行左、右端项
	*/
	for (i = 0; i < originalRowNum; i++)
	{
		CoinPackedVector ithVec = originalMatrixByRow->getVector(i);
		double maxCoeff = ithVec.infNorm();

		if (transedRowLower[i] != -COIN_DBL_MAX)
			transedRowLower[i] /= maxCoeff;
		if (transedRowUpper[i] != COIN_DBL_MAX)
			transedRowUpper[i] /= maxCoeff;

		const int* indices = ithVec.getIndices();
		const double* elements = ithVec.getElements();  // 只包含非0元素
		int eleNum = ithVec.getNumElements();

		for (j = 0; j < eleNum; j++)
		{
			double newElement = elements[j] / maxCoeff;
			int iCol = indices[j];
			transedMatrixByCol.modifyCoefficient(i, iCol, newElement);
		}
	}

	int integerNum = model_->numberIntegers();
	int continuousNum = originalColNum - integerNum;
	const double* originalColLower = model_->getColLower();
	const double* originalColUpper = model_->getColUpper();
	const double* coefficients = model_->getObjCoefficients();
	double* transedCoeffs = new double[integerNum];
	int* continousIndices = new int[continuousNum];
	double* transedColLower = new double[integerNum];
	double* transedColUpper = new double[integerNum];
	int coeffCount = 0;
	int continousCount = 0;
	int colBoundCount = 0;
	/*	1. 松弛连续变量，因为整数变量均为0-1变量，所以不需要处理
		2. 得到转换后目标函数系数、变量上下界（此处由于均为0-1，其实可以直接上界全设置为1，下界全设置为0）
		3. 统计连续变量所在索引，最后需要在矩阵中将这些变量删除掉
	*/
	for (j = 0; j < originalColNum; j++)
	{
		if (model_->isInteger(j))
		{
			transedCoeffs[coeffCount++] = coefficients[j];
			transedColLower[colBoundCount] = originalColLower[j];
			transedColUpper[colBoundCount] = originalColUpper[j];
			colBoundCount++;
		}
		else
		{
			relaxContinuousVar(&transedMatrixByCol, originalColUpper[j], originalColLower[j], transedRowUpper, transedRowLower, j);
			continousIndices[continousCount++] = j;
		}
	}
	transedMatrixByCol.deleteCols(continousCount, continousIndices);
	OsiSolverInterface* transedSolver = model_->solver()->clone();
	transedSolver->loadProblem(transedMatrixByCol, transedColLower, transedColUpper, transedCoeffs, transedRowLower, transedRowUpper);

	for (i = 0; i < transedSolver->getNumCols(); i++)
		transedSolver->setInteger(i);
	delete[] transedColLower;
	delete[] transedColUpper;
	delete[] transedRowLower;
	delete[] transedRowUpper;
	delete[] transedCoeffs;
	delete[] continousIndices;

	return transedSolver;
}

/*compare function for variables' importance*/
bool 
ShiftAndPropagate::cmpForVarImportance(std::pair<int, double>a, std::pair<int, double>b)
{
	if (a.second == b.second)
		return a.first < b.first;
	return a.second > b.second;
}

int 
ShiftAndPropagate::bestShift(OsiSolverInterface* transedProblem, int iCol, double* currentRowLower, double* currentRowUpper)
{
	const CoinPackedMatrix* matrixByCol = transedProblem->getMatrixByCol();
	double ithColUpper = transedProblem->getColUpper()[iCol];
	double ithColLower = transedProblem->getColLower()[iCol];
	const CoinPackedVector ithColVector = matrixByCol->getVector(iCol);
	const double* ithColElements = ithColVector.getElements();
	const int* ithColIndices = ithColVector.getIndices();
	int ithColElementsNum = ithColVector.getNumElements();
	std::set<std::pair<int, int>>Q;
	for (int i = 0; i < ithColElementsNum; i++)
	{
		int iRow = ithColIndices[i];
		double rlb = currentRowLower[iRow];
		double rub = currentRowUpper[iRow];
		double aij = ithColElements[i];
		if (rub != COIN_DBL_MAX)
		{
			if (rub < 0 && aij < 0)
			{
				int t = std::ceil(rub / aij);
				if (ithColLower <= t && t <= ithColUpper)
					Q.insert(std::pair<int, int>(t, -1));
			}
			if (rub >= 0 && aij > 0)
			{
				int t = std::ceil(rub / aij + 1e-5);   // 为了避免rub / aij 刚好为整数的情况
				if (ithColLower <= t && t <= ithColUpper)
					Q.insert(std::pair<int, int>(t, 1));
			}
		}

		if (rlb != -COIN_DBL_MAX)
		{
			if (rlb > 0 && aij > 0)
			{
				int t = std::ceil(rlb / aij);
				if (ithColLower <= t && t <= ithColUpper)
					Q.insert(std::pair<int, int>(t, -1));
			}
			if (rlb <= 0 && aij < 0)
			{
				int t = std::ceil(rlb / aij + 1e-5);  // 为了避免rlb / aij刚好为整数情况
				if (rlb <= t && t <= rlb)
					Q.insert(std::pair<int, int>(t, 1));
			}
		}
	}

	if (Q.empty())  // 这个地方文章中说明这样更改
		return std::floor(ithColUpper);

	int sigma = 0;
	int bestVal = 0;
	int tBefore = 0;
	int minRowViolationSum = 0;
	for (auto it = Q.begin(); it != Q.end(); it++)
	{
		if ((*it).first == tBefore)
			sigma += (*it).second;
		else
		{
			if (sigma < minRowViolationSum)
			{
				minRowViolationSum = sigma;
				bestVal = tBefore;
			}
			tBefore = (*it).first;
			sigma = 0;
			sigma += (*it).second;
		}
		auto iter = Q.end(); // 判断是否为最后一个元素
		iter--;
		if (it == iter && sigma < minRowViolationSum)
		{
			minRowViolationSum = sigma;
			bestVal = tBefore;
		}
	}
	return bestVal;
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
	OsiSolverInterface* transedProblem0 = transformProblem(model_);
	OsiClpSolverInterface* transedProblem = dynamic_cast<OsiClpSolverInterface*>(transedProblem0);

	// Order of variable importance
	int transedRowNum = transedProblem->getNumRows();
	int transedColNum = transedProblem->getNumCols();
	const CoinPackedMatrix* transedMatrixByCol = transedProblem->getMatrixByCol();
	std::vector<std::pair<int, double>>varImportance;
	for (int j = 0; j < transedColNum; j++)
	{
		CoinPackedVector ithColVector = transedMatrixByCol->getVector(j);
		double* ithColElements = ithColVector.getElements();
		int ithColElementsNum = ithColVector.getNumElements();

		int nonZeroCoefNum = ithColElementsNum;
		double absCoefSum = 0.0;
		for (int i = 0; i < ithColElementsNum; i++)
			absCoefSum += std::fabs(ithColElements[i]);
		double importance = nonZeroCoefNum + absCoefSum;
		varImportance.push_back(std::pair<int, double>(j, importance));
	}
	std::sort(varImportance.begin(), varImportance.end(), cmpForVarImportance);

	/***********main loop *************/
	double* bestShiftForCols = new double[transedColNum](); // 用于存储每个整数变量的best shift，对于循环中没有确定的，皆取best shift为下界即可
	// double* constant = new double[transedRowNum](); // 用于记录每次固定变量后所产生的常量值，用于判断zero-solution的可行性。  
	bool isSatisfy = true;
	int backtrackNum = 0;   // 记录回溯的次数
	const double* transedRowLower = transedProblem->getRowLower();
	const double* transedRowUpper = transedProblem->getRowUpper();
	double* currentRowLower = new double[transedRowNum];   // 用于存储实时行上下界，因为在某一变量固定后，计算新的上下界用于判断0解的可行性
	double* currentRowUpper = new double[transedRowNum];
	std::memcpy(currentRowLower, transedRowLower, sizeof(double) * transedRowNum);
	std::memcpy(currentRowUpper, transedRowUpper, sizeof(double) * transedRowNum);

	for (auto it = varImportance.begin(); it != varImportance.end(); it++)
	{
		if (backtrackNum > maxBacktrackLimit_)
			break;
		int iCol = (*it).first;
		bool isLastEqual = false;  // 表示是否为最后一列，且上下界已经相等
		// 获取当前列上下界，因为propagate可能会改变列的上下界
		const double* currentColLower = transedProblem->getColLower();
		const double* currentColUpper = transedProblem->getColUpper();
		if (currentColLower[iCol] == currentColUpper[iCol])
		{
			bestShiftForCols[iCol] = currentColLower[iCol];
			CoinPackedVector ithColVector = transedMatrixByCol->getVector(iCol);
			double* ithColElements = ithColVector.getElements();
			int* ithColIndices = ithColVector.getIndices();
			int ithElementsNum = ithColVector.getNumElements();
			for (int i = 0; i < ithElementsNum; i++)
			{
				int iRow = ithColIndices[i];
				double element = ithColElements[i];
				double constant = element * currentColLower[iCol];
				if (currentRowLower[iRow] != -COIN_DBL_MAX)
					currentRowLower[iRow] -= constant;
				if (currentRowUpper[iRow] != COIN_DBL_MAX)
					currentRowUpper[iRow] -= constant;
			}

			auto iter = varImportance.end();
			iter--;
			if (it != iter)
				continue;
			isLastEqual = true;
		}
		// 判断0解是否满足
		for (int i = 0; i < transedRowNum; i++)
		{
			if (currentRowLower[i] > 0 || currentRowUpper[i] < 0)  // 如果存在某行下界大于0或者某行上界小于0，不满足
			{
				isSatisfy = false;
				break;
			}
			isSatisfy = true;
		}
		if (isSatisfy || isLastEqual)
			break;

		int bestVal = bestShift(transedProblem, iCol, currentRowLower, currentRowUpper);

		double* backtrackColLower = new double[transedColNum];
		double* backtrackColUpper = new double[transedColNum];
		memcpy(backtrackColLower, currentColLower, sizeof(double) * transedColNum);
		memcpy(backtrackColUpper, currentColUpper, sizeof(double) * transedColNum);

		transedProblem->setColLower(iCol, bestVal);
		transedProblem->setColUpper(iCol, bestVal);
		int tightenNum = transedProblem->tightenBounds();
		if (tightenNum == -1)
		{
			transedProblem->setColLower(backtrackColLower);
			transedProblem->setColUpper(backtrackColUpper);
			backtrackNum++;
		}
		else
		{
			bestShiftForCols[iCol] = bestVal;
			CoinPackedVector ithColVector = transedMatrixByCol->getVector(iCol);
			double* ithColElements = ithColVector.getElements();
			int* ithColIndices = ithColVector.getIndices();
			int ithElementsNum = ithColVector.getNumElements();
			for (int i = 0; i < ithElementsNum; i++)
			{
				int iRow = ithColIndices[i];
				double element = ithColElements[i];
				double constant = element * bestVal;
				if (currentRowLower[iRow] != -COIN_DBL_MAX)
					currentRowLower[iRow] -= constant;
				if (currentRowUpper[iRow] != COIN_DBL_MAX)
					currentRowUpper[iRow] -= constant;
			}
		}
		delete[] backtrackColLower;
		delete[] backtrackColUpper;
	}
	int returnCode = isSatisfy;
	if (isSatisfy)
	{
		const int* integerIndices = model_->integerVariable();
		OsiSolverInterface* solver = model_->solver()->clone();
		for (int j = 0; j < transedColNum; j++)
		{
			int iCol = integerIndices[j];
			solver->setColLower(iCol, bestShiftForCols[j]);
			solver->setColUpper(iCol, bestShiftForCols[j]);
		}
		solver->initialSolve();
		if (solver->isProvenOptimal())
		{
			const double* solverSol = solver->getColSolution();
			memcpy(betterSolution, solverSol, sizeof(double) * columnNumber);
			double solverObjVal = solver->getObjValue();
			solutionValue = solverObjVal;
		}
	}
	return returnCode;
}
