#ifndef VRINS_H
#define VRINS_H

#include "CbcHeuristic.hpp"
// for backward compatibility include 3 other headers
#include "CbcHeuristicRENS.hpp"
#include "CbcHeuristicDINS.hpp"
#include "CbcHeuristicVND.hpp"
/** LocalSearch class
 */

class VRINS : public CbcHeuristic
{
public:

	// Default Constructor
	VRINS();

	/* Constructor with model - assumed before cuts
	   Initial version does not do Lps
	*/
	VRINS(CbcModel& model);

	// Copy constructor
	VRINS(const VRINS&);

	// Destructor
	~VRINS();

	/// Clone
	virtual CbcHeuristic* clone() const;


	/// Assignment operator
	VRINS& operator=(const VRINS& rhs);

	/// Create C++ lines to get to current state
	virtual void generateCpp(FILE* fp);

	/// Resets stuff if model changes
	virtual void resetModel(CbcModel* model);

	/// update model (This is needed if cliques update matrix etc)
	virtual void setModel(CbcModel* model);

	using CbcHeuristic::solution;
	/** returns 0 if no solution, 1 if valid solution.
		Sets solution values if good, sets objective value (only if good)
		This does Relaxation Induced Neighborhood Search
	*/
	virtual int solution(double& objectiveValue,
		double* newSolution);
	/// This version fixes stuff and does IP
	int solutionFix(double& objectiveValue,
		double* newSolution,
		const int* keep);

	/// Sets how often to do it
	inline void setHowOften(int value)
	{
		howOften_ = value;
	}
	/// Used array so we can set
	inline char* used() const
	{
		return used_;
	}
	/// Resets lastNode
	inline void setLastNode(int value)
	{
		lastNode_ = value;
	}
	/// Resets number of solutions
	inline void setSolutionCount(int value)
	{
		numberSolutions_ = value;
	}

protected:
	// Data

	/// Number of solutions so we can do something at solution
	int numberSolutions_;
	/// How often to do (code can change)
	int howOften_;
	/// Number of successes
	int numberSuccesses_;
	/// Number of tries
	int numberTries_;
	/** State of fixing continuous variables -
		0 - not tried
		+n - this divisor makes small enough
		-n - this divisor still not small enough
	*/
	int stateOfFixing_;
	/// Node when last done
	int lastNode_;
	/// Whether a variable has been in a solution
	char* used_;
};
#endif


