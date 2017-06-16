OBJS = CoinAlloc.o \
	CoinBuild.o \
	CoinDenseFactorization.o \
	CoinDenseVector.o \
	CoinError.o \
	CoinFactorization1.o \
	CoinFactorization2.o \
	CoinFactorization3.o \
	CoinFactorization4.o \
	CoinFileIO.o \
	CoinFinite.o \
	CoinIndexedVector.o \
	CoinLpIO.o \
	CoinMessage.o \
	CoinMessageHandler.o \
	CoinModel.o \
	CoinModelUseful.o \
	CoinModelUseful2.o \
	CoinMpsIO.o \
	CoinOslFactorization.o \
	CoinOslFactorization2.o \
	CoinOslFactorization3.o \
	CoinPackedMatrix.o \
	CoinPackedVector.o \
	CoinPackedVectorBase.o \
	CoinParam.o \
	CoinParamUtils.o \
	CoinPostsolveMatrix.o \
	CoinPrePostsolveMatrix.o \
	CoinPresolveDoubleton.o \
	CoinPresolveDual.o \
	CoinPresolveDupcol.o \
	CoinPresolveEmpty.o \
	CoinPresolveFixed.o \
	CoinPresolveForcing.o \
	CoinPresolveHelperFunctions.o \
	CoinPresolveImpliedFree.o \
	CoinPresolveIsolated.o \
	CoinPresolveMatrix.o \
	CoinPresolveMonitor.o \
	CoinPresolvePsdebug.o \
	CoinPresolveSingleton.o \
	CoinPresolveSubst.o \
	CoinPresolveTighten.o \
	CoinPresolveTripleton.o \
	CoinPresolveUseless.o \
	CoinPresolveZeros.o \
	CoinRational.o \
	CoinSearchTree.o \
	CoinShallowPackedVector.o \
	CoinSimpFactorization.o \
	CoinSnapshot.o \
	CoinStructuredModel.o \
	CoinWarmStartBasis.o \
	CoinWarmStartDual.o \
	CoinWarmStartPrimalDual.o \
	CoinWarmStartVector.o

all: $(OBJS)
	$(AR) cru libCoinUtils.a $(OBJS)

clean:
	$(RM) *.o *.a

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -DHAVE_CONFIG_H -DCOINUTILS_BUILD -c $< -o $@
