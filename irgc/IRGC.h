/* 
	This code implements the IRGC and IRGC+expansion algorithms
	described in the paper 

	"Iteratively Reweighted Graph Cut for Multi-label MRFs with Non-convex Priors", 
	Thalaiyasingam Ajanthan, Richard Hartley, Mathieu Salzmann and Hongdong Li,
	IEEE Conference on Computer Vision and Pattern Recognition,
	June 2015.

	*Code Assumptions:
		1. The MRF energy has the following form
			E(x) = \sum \theta_{i}(x_i) + \sum \theta_{ij} (x_i, x_j), 
			where \theta_{ij} (x_i, x_j) = \gamma_{ij} \theta(|x_i - x_j|).
		2. The code currently supports MRF with 4-connected grid structure only
			nodes labelled from 0 --> width * height - 1, in a grid structure.k
			E.g. width = 3, height = 2
			0 -- 1 -- 2
			|	 |	  |
			3 -- 4 -- 5

	If you use this code, please consider citing the aforementioned paper 
	in any resulting publication.
	
	*Contact: thalaiyasingam.ajanthan@nicta.com.au
	
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>

#include "maxflow/graph_irgc.h"
#include "alpha_exapnsion_west/GCoptimization.h"

#define DEBUG 1
#define INFTY 100000000	
#define INIT_W 0.5

/* possible data types
	nodeid  - short, int, long (signed and unsigned)
	labelid - short, int, long (signed only) 
	captype - int, float, double (signed only) 
	==> possible instances are in irgc_instances.inc */
template <typename nodeid, typename labelid, typename captype> 
class IRGC
{
public:
	
	typedef GraphIrgc<captype, captype, captype> IshikawaGraph;
	typedef GCoptimizationGridGraph ExpansionGraph;
	
	/* Constructor 
		width - image width, height - image height, labels - number of labels, 
		hybrid - run the IRGC+expansion, logf - log file name */
	IRGC(nodeid width, nodeid height, labelid labels, bool hybrid = false, const char* logf = NULL);

	/* Destructor */
	~IRGC();

	/* Set the MRF energy		
		unaryPotential(i, li)  - returns the unary potential \theta_{i}(li)
		binaryWeights(i, j)    - returns the constant binary weights \gamma_{ij}
		binaryFunction(li, lj) - returns the smoothness cost \theta(li, lj) 
		binaryPotential(i, j, li, lj) - returns binaryWeights(i, j) * binaryFunction(li, lj) 
		binary - [TL, TQ, CA] 
		lambda > 0 */
	void setEnergy(captype(*unaryPotential)(nodeid, labelid), captype(*binaryWeights)(nodeid, nodeid), 
		captype(*binaryFunction)(labelid, labelid), captype(*binaryPotential)(nodeid, nodeid, labelid, labelid), 
		char* binary, labelid lambda);

	/* Run algorithm
		maxiter - maximum number of iterations 
		returns true if algorithm is converged */
	bool optimize(size_t maxiter = 100);

	/* Get the label
		returns the label of the given node */
	labelid getLabel(nodeid ii)
	{
		return labelling[ii];
	}

private:
	nodeid width, height;
	labelid labels, lambda;
	std::ofstream fout;
	bool hybrid;
	nodeid imagexy, numnodes, numedges;
	size_t iter;
	labelid* labelling;

	captype(*unaryPotential)(nodeid, labelid);
	captype(*binaryWeights)(nodeid, nodeid);
	captype(*binaryFunction)(labelid, labelid);
	captype(*binaryPotential)(nodeid, nodeid, labelid, labelid);

	enum BFunction {TL, TQ, CA} bf;

	IshikawaGraph *IG;
	ExpansionGraph *EG;

	void errFunction(const char* msg)
	{
		std::cout << endl << "Error!\n " << msg << std::endl;
		fout << endl << "Error!\n " << msg << std::endl;
		exit(1);
	}

	// cauchy loss function
	float cauchy(float x)
	{
		return std::pow(lambda, 2) * log(1 + std::pow((x / lambda), 2)) / 2;
	}

	void init();
	captype computeEnergy();
	void minimumCut();
	captype binaryEdgeCapacities(nodeid, nodeid, labelid, labelid);
	void addCrossEdges(nodeid, nodeid);
	void construct();
	double irgcIteration();
	void decompile();

	double alphaIteration();
};

