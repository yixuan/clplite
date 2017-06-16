OBJS = ClpCholeskyBase.o \
	ClpCholeskyDense.o \
	ClpConstraint.o \
	ClpConstraintLinear.o \
	ClpConstraintQuadratic.o \
	Clp_C_Interface.o \
	ClpDualRowDantzig.o \
	ClpDualRowPivot.o \
	ClpDualRowSteepest.o \
	ClpDummyMatrix.o \
	ClpDynamicExampleMatrix.o \
	ClpDynamicMatrix.o \
	ClpEventHandler.o \
	ClpFactorization.o \
	ClpGubDynamicMatrix.o \
	ClpGubMatrix.o \
	ClpHelperFunctions.o \
	ClpInterior.o \
	ClpLinearObjective.o \
	ClpMatrixBase.o \
	ClpMessage.o \
	ClpModel.o \
	ClpNetworkBasis.o \
	ClpNetworkMatrix.o \
	ClpNonLinearCost.o \
	ClpNode.o \
	ClpObjective.o \
	ClpPackedMatrix.o \
	ClpPlusMinusOneMatrix.o \
	ClpPredictorCorrector.o \
	ClpPdco.o \
	ClpPdcoBase.o \
	ClpLsqr.o \
	ClpPresolve.o \
	ClpPrimalColumnDantzig.o \
	ClpPrimalColumnPivot.o \
	ClpPrimalColumnSteepest.o \
	ClpQuadraticObjective.o \
	ClpSimplex.o \
	ClpSimplexDual.o \
	ClpSimplexNonlinear.o \
	ClpSimplexOther.o \
	ClpSimplexPrimal.o \
	ClpSolve.o \
	Idiot.o \
	IdiSolve.o

all: $(OBJS)
	$(AR) cru libClp.a $(OBJS)

clean:
	$(RM) *.o *.a

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -DHAVE_CONFIG_H -DCLP_BUILD -I../CoinUtils -c $< -o $@
