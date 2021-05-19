// This is the shift and propagate primal heuristic
#ifndef ShifAndPropagate_H
#define ShifAndPropagate_H

#include "CbcHeuristic.hpp"

class ShiftAndPropagate : public CbcHeuristic
{
public:

    /// Default Constructor
    ShiftAndPropagate();

    /// Constructor with model 
    ShiftAndPropagate(CbcModel& model);

    /// Copy constructor
    ShiftAndPropagate(const ShiftAndPropagate&);

    /// Destructor
    ~ShiftAndPropagate();

    /// Clone
    virtual CbcHeuristic* clone() const;

    /// Assignment operator
    ShiftAndPropagate& operator=(const ShiftAndPropagate& rhs);

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel* model);

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel* model);

    /// relax continuous variable
    void relaxContinuousVar(CoinPackedMatrix* matrixByCol, double colUpper, double colLower, double* rowUpper, double* rowLower, int iCol);

    /// do problem transformation
   virtual OsiSolverInterface* transformProblem(CbcModel* model_);

   /// decide best shift
   int bestShift(OsiSolverInterface* transedProblem, int iCol, double* currentRowLower, double* currentRowUpper);

   /*compare function for variables' importance*/
   bool cmpForVarImportance(std::pair<int, double>a, std::pair<int, double>b);

    using CbcHeuristic::solution;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        This does Relaxation Induced Neighborhood Search
    */
    virtual int solution(double& objectiveValue,
        double* newSolution);

    /// Sets how often to do it
    inline void setMaxBacktrack(int value)
    {
        maxBacktrackLimit_ = value;
    }
    /// Gets how often to do it
    inline int getHowOften()
    {
        return maxBacktrackLimit_;
    }


protected:
    /// max backtrack limit
    int maxBacktrackLimit_;

};
#endif


