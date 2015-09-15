
#include <ctime>

#include "IRGC.h"

template <typename nodeid, typename labelid, typename captype>
IRGC<nodeid, labelid, captype>::IRGC(nodeid width, nodeid height, labelid labels,
	bool hybrid, const char* logf)
	: width(width), height(height), labels(labels), hybrid(hybrid),
	lambda(0), numnodes(0), numedges(0), iter(0), IG(NULL)
{
	if (!logf) {
		if (hybrid) fout.open("irgc_exp.out");
		else fout.open("irgc.out");
	}
	else fout.open(logf);
	imagexy = width * height;

	if (hybrid) EG = new ExpansionGraph(width, height, labels);
	labelling = new labelid[imagexy]();
}

template <typename nodeid, typename labelid, typename captype>
IRGC<nodeid, labelid, captype>::~IRGC()
{
	delete[] labelling;
	if (IG) delete IG;
	if (EG) delete EG;
	fout.close();
}

template <typename nodeid, typename labelid, typename captype>
void IRGC<nodeid, labelid, captype>::setEnergy(captype(*unaryPotential)(nodeid, labelid),
	captype(*binaryWeights)(nodeid, nodeid), captype(*binaryFunction)(labelid, labelid),
	captype(*binaryPotential)(nodeid, nodeid, labelid, labelid), char* binary, labelid lambda)
{
	this->unaryPotential = unaryPotential;
	this->binaryWeights = binaryWeights;
	this->binaryFunction = binaryFunction;
	this->binaryPotential = binaryPotential;

	if (!binary) errFunction("Invalid binary function type!");
	else if (!strcmp(binary, (char*)"TL")) bf = TL;	// truncated linear
	else if (!strcmp(binary, (char*)"TQ")) bf = TQ;	// truncated quadratic
	else if (!strcmp(binary, (char*)"CA")) bf = CA;	// cauchy
	else errFunction("Invalid binary function type!");
	if (lambda <= 0) errFunction("Invalid lambda!");
	this->lambda = lambda;

	if (hybrid) {
		// set up the needed data to pass to function for the data costs
		EG->setDataCost(unaryPotential);
		// smoothness comes from function pointer
		EG->setSmoothCost(binaryPotential);
	}
}

template <typename nodeid, typename labelid, typename captype>
void IRGC<nodeid, labelid, captype>::init()
{
	if (lambda <= 0) errFunction("Invalid energy!");

	numnodes = imagexy * (labels - 1);
	numedges = imagexy * (labels - 2);	// vertical edges
	nodeid numedges_p = (width - 1) * height + (height - 1) * width;
	switch (bf) {
	case TL:
		numedges += numedges_p * (labels - 1);	// horizontal
		break;
	case TQ:
	case CA:
		numedges += numedges_p * lambda * (2 * labels - 3 - lambda);	// cross-edges - horizontal
		numedges += numedges_p * (labels - 1);	// horizontal
		break;
	}

	IG = new IshikawaGraph(numnodes, numedges);
}

template <typename nodeid, typename labelid, typename captype>
captype IRGC<nodeid, labelid, captype>::computeEnergy()
{
	captype energy = 0;

	for (nodeid i = 0; i < height; ++i) {
		for (nodeid j = 0; j < width; ++j) {
			nodeid ii = i * width + j;
			energy += unaryPotential(ii, labelling[ii]);	//unary potential
			if (j != width - 1) {	// right
				nodeid jj = ii + 1;
				energy += binaryPotential(ii, jj, labelling[ii], labelling[jj]);
			}
			if (i != height - 1) {	// bottom
				nodeid jj = ii + width;
				energy += binaryPotential(ii, jj, labelling[ii], labelling[jj]);
			}
		}
	}
	return energy;
}

template <typename nodeid, typename labelid, typename captype>
void IRGC<nodeid, labelid, captype>::minimumCut()
{
	for (nodeid i = 0; i < imagexy; ++i) {	// image grid
		nodeid ii = i;
		bool marked = false;
		for (labelid k = 1; k < labels; ++k) {	// label dimension
			if (!marked && IG->what_segment(ii)) {
				labelling[i] = k - 1; // no of zeros
				marked = true;
			}
			assert((marked && IG->what_segment(ii)) || (!marked && !IG->what_segment(ii)));	// only one vertical edge can be cut
			ii += imagexy;
		}
		if (!marked) labelling[i] = labels - 1;
	}
}

template <typename nodeid, typename labelid, typename captype>
captype IRGC<nodeid, labelid, captype>::binaryEdgeCapacities(nodeid ii, nodeid jj, labelid li, labelid lj)
{
	captype w = binaryWeights(ii, jj);
	labelid delta = std::abs(labelling[ii] - labelling[jj]);
	labelid labelDiff = std::abs(li - lj);
	float derivative = INIT_W;
	switch (bf) {
	case TL:{	// truncated linear
		if (iter) {
			if (delta == 0) derivative = 1;
			else if (delta < lambda) derivative = 1;
			else if (delta == lambda) derivative = 0.5;
			else derivative = 0;
		}
		if (li == lj) return (captype)w * derivative;
		errFunction("Invalid parameters for binaryEdgeCapacities for truncated-linear!");
		break;
	}
	case TQ:{	// truncated quadratic
		if (iter) {
			if (delta < lambda) derivative = 1;
			else if (delta == lambda) derivative = 0.5;
			else derivative = 0;
		}
		if (labelDiff == 0) return (captype)w * derivative;
		else if (labelDiff < lambda) return (captype)w * derivative * 2;
		else if (labelDiff == lambda) return (captype)w * derivative;
		errFunction("Invalid parameters for binaryEdgeCapacities for truncated-quadratic!");
		break;
	}
	case CA:{	// cauchy function
		if (iter) {
			if (delta <= lambda) derivative = 1;
			//else if(delta == lambda) derivative = 1;
			else derivative = 2 * lambda * delta / (std::pow(lambda, 2) + std::pow(delta, 2));//(cauchy(delta)-cauchy(delta-1))/(cauchy(lambda+1)-cauchy(lambda));//
		}
		//float sd = pow(lambda, 2)*(pow(lambda, 2)-pow(labelDiff, 2))/pow((pow(lambda, 2)+pow(labelDiff, 2)),2);
		//captype sd = binaryFunction(li + 1, lj) + binaryFunction(li, lj + 1) - binaryFunction(li, lj) - binaryFunction(li + 1, lj + 1);
		float sd = cauchy(labelDiff + 1) + cauchy(labelDiff - 1) - 2 * cauchy(labelDiff);
		if (sd * 100000000 < 0) sd = 0;
		assert(sd >= 0);
		if (labelDiff == 0) return (captype)w * derivative * sd / 2.0;
		else if (labelDiff < lambda) return (captype)w * derivative * sd;
		else if (labelDiff == lambda) return (captype)w * derivative * sd;
		errFunction("Invalid parameters for binaryEdgeCapacities for cauchy!");
		break;
	}
	}
	return 0;
}

template <typename nodeid, typename labelid, typename captype>
void IRGC<nodeid, labelid, captype>::addCrossEdges(nodeid ii, nodeid jj)
{
	nodeid iiOriginal = ii;
	nodeid jjOriginal = jj;
	for (labelid k = 1; k < labels; ++k) {
		nodeid ii_o = ii;
		nodeid jj_o = jj;
		switch (bf) {
		case TL:{
			captype w = binaryEdgeCapacities(iiOriginal, jjOriginal, k, k);
			if (w < 0) w = 0;
			IG->add_edge(ii_o, jj_o, w, w);	// bidirectional horizontal edges
			break;
		}
		case TQ:
		case CA:{
			// left to right (upward)
			for (labelid l = k + 1; (l <= k + lambda && l < labels); ++l) {	// l > k
				ii += imagexy;
				captype w = binaryEdgeCapacities(iiOriginal, jjOriginal, l, k);
				if (w < 0) w = 0;
				IG->add_edge(ii, jj_o, w, 0);	// other cross edges (ii-->jj)
			}

			// right to left (upward)
			for (labelid l = k + 1; (l <= k + lambda && l < labels); ++l) {	// l > k
				jj += imagexy;
				captype w = binaryEdgeCapacities(iiOriginal, jjOriginal, l, k);
				if (w < 0) w = 0;
				IG->add_edge(jj, ii_o, w, 0);	// other cross edges (jj-->ii)
			}

			captype w = binaryEdgeCapacities(iiOriginal, jjOriginal, k, k);
			if (w < 0) w = 0;
			IG->add_edge(ii_o, jj_o, w, w);	// bidirectional horizontal edges
			break;
		}
		}
		ii = ii_o + imagexy;
		jj = jj_o + imagexy;
	}
}

template <typename nodeid, typename labelid, typename captype>
void IRGC<nodeid, labelid, captype>::construct()
{
	IG->add_node(numnodes);	// add all nodes at once 

	// unary potentials
	for (nodeid i = 0; i < imagexy; ++i) {	// image grid
		nodeid ii = i;
		for (labelid k = 1; k < labels; ++k) {	// label dimension
			if (k == labels - 1) {
				IG->add_tweights(ii, 0, unaryPotential(i, k));	// sink
				break;
			}
			else if (k == 1) {
				IG->add_tweights(ii, unaryPotential(i, 0), 0);	// source
			}

			nodeid jj = ii + imagexy;
			IG->add_edge(ii, jj, unaryPotential(i, k), INFTY);	// vertical edges
			ii += imagexy;
		}
	}

	// binary potentials
	for (nodeid i = 0; i < height; ++i) {
		for (nodeid j = 0; j < width; ++j) {
			nodeid ii = i * width + j;
			if (j != width - 1) {	// right
				nodeid jj = ii + 1;
				addCrossEdges(ii, jj);
			}
			if (i != height - 1) {	// bottom
				nodeid jj = ii + width;
				addCrossEdges(ii, jj);
			}
		}
	}
}

template <typename nodeid, typename labelid, typename captype>
double IRGC<nodeid, labelid, captype>::irgcIteration()
{
	if (!IG) errFunction("\nIshikawa graph not created!!");

	// construct graph
	construct();
#if DEBUG >= 2
	std::cout << "\nGet Nodes: " << IG->get_node_num() << ", Get Edges: " << IG->get_arc_num() << std::endl;
	fout << "\nGet Nodes: " << IG->get_node_num() << ", Get Edges: " << IG->get_arc_num() << std::endl;
#endif

	// run maxflow
	clock_t startTime = clock();
	IG->maxflow();
	clock_t endTime = clock();
	double timeTaken = ((double)endTime - (double)startTime) / CLOCKS_PER_SEC;

	// compute min-cut
	minimumCut();
	// reset graph
	IG->reset();

	return timeTaken;
}

template <typename nodeid, typename labelid, typename captype>
double IRGC<nodeid, labelid, captype>::alphaIteration()
{
	if (!EG) errFunction("Alpha expansion graph not created!");

	for (nodeid i = 0; i < imagexy; ++i) {	// initialize
		EG->setLabel(i, labelling[i]);
	}

	// set random label order
	//EG->setLabelOrder(true);

	// run expansion
	clock_t startTime = clock();
	EG->expansion(1);
	clock_t endTime = clock();
	double timeTaken = ((double)endTime - (double)startTime) / CLOCKS_PER_SEC;

	// compute min-cut
	for (nodeid i = 0; i < imagexy; ++i) {
		labelling[i] = EG->whatLabel(i);
	}

	return timeTaken;
}

template <typename nodeid, typename labelid, typename captype>
void IRGC<nodeid, labelid, captype>::decompile()
{
	if (!IG) return;

	std::cout << "Graph: " << height << "x" << width << "x" << labels << std::endl;
	fout << "Graph: " << height << "x" << width << "x" << labels << std::endl;
	std::cout << "Binary Function: ";
	fout << "Binary Function: ";
	switch (bf) {
	case TL:  	   std::cout << "Truncated Linear (lambda=" << lambda << ")\n"; fout << "Truncated Linear (lambda=" << lambda << ")\n"; break;
	case TQ:  	   std::cout << "Truncated Quadratic (lambda=" << lambda << ")\n"; fout << "Truncated Quadratic (lambda=" << lambda << ")\n"; break;
	case CA:  	   std::cout << "Cauchy Function (lambda=" << lambda << ")\n"; fout << "Cauchy Function (lambda=" << lambda << ")\n"; break;
	}

	std::cout << "Nodes: " << numnodes << ", Edges: " << numedges << std::endl;
	fout << "Nodes: " << numnodes << ", Edges: " << numedges << std::endl;

	fout << "\nMinimum cut\n";
	for (nodeid i = 0; i < imagexy; ++i) {	// image grid
		fout << i << ": " << labelling[i] << std::endl;
	}

#if DEBUG >= 3
	std::cout << "\nSource-Sink capacities\n";
	fout << "\nSource-Sink capacities\n";
	for (nodeid i = 0; i < imagexy; ++i) {	// image grid
		nodeid ii = i;
		for (labelid k = 1; k < labels; ++k) {	// label dimension
			if (k == labels - 1) {
				std::cout << ", " << ii << ": " << g->get_trcap(ii);
				fout << ", " << ii << ": " << g->get_trcap(ii);
				break;
			}
			else if (k == 1) {
				std::cout << ", " << ii << ": " << g->get_trcap(ii);
				fout << ", " << ii << ": " << g->get_trcap(ii);
				ii += imagexy;
				continue;
			}
			std::cout << ", " << ii << ": " << g->get_trcap(ii);
			fout << ", " << ii << ": " << g->get_trcap(ii);
			ii += imagexy;
		}
		std::cout << std::endl;
		fout << std::endl;
	}

	std::cout << "\nArc capacities\n";
	fout << "\nArc capacities\n";
	GraphType::arc_id a = g->get_first_arc();
	for (nodeid e = 0; e < numedges; ++e) {
		nodeid i = 0, j = 0;
		capType rc = 0, revRc = 0;
		g->get_arc_ends(a, i, j);
		g->get_rcap(a, rc, revRc);
		std::cout << ", " << i << "->" << j << ": (" << rc << ", " << revRc << ")";
		fout << ", " << i << "->" << j << ": (" << rc << ", " << revRc << ")";
		//a = g->get_next_arc(a);
		//std::cout << ", " << j << "->" << i << ": " << g->get_rcap(a);
		//fout << ", " << j << "->" << i << ": " << g->get_rcap(a);
		a = g->get_next_arc(a);
		std::cout << "\n";
		fout << "\n";
	}
#endif
}

template <typename nodeid, typename labelid, typename captype>
bool IRGC<nodeid, labelid, captype>::optimize(size_t maxiter)
{
	if (hybrid) {
		std::cout << "\n## IRGC + expansion ##\n";
		fout << "\n## IRGC + expansion ##\n";
	}
	else {
		std::cout << "\n## IRGC ##\n";
		fout << "\n## IRGC ##\n";
	}
	init();

	captype energy = computeEnergy();
	std::cout << "Init-Energy = " << energy << "\n";
	fout << "Init-Energy = " << energy << "\n";

	captype oldEnergy = energy;
	double timeTaken = 0;
	size_t count = 0;

	for (size_t i = 1; i <= maxiter; ++i) {
		double t = 0;
#if DEBUG >= 1
		std::cout << "\nIteration: " << i;
		fout << "\nIteration: " << i;
#endif
		if (hybrid) {
			if (i % 2) t = irgcIteration();	// first iteration is here
			else t = alphaIteration();	// even numbers
		}
		else {
			t = irgcIteration();
		}
		timeTaken += t;
		++iter;

		energy = computeEnergy();
#if DEBUG >= 1
		std::cout << " (energy = " << energy << ", time = [" << t << ", " << timeTaken << "])";
		fout << " (energy = " << energy << ", time = [" << t << ", " << timeTaken << "])";
#endif
		assert(round(oldEnergy) >= round(energy) || i == 1);
		if (abs(oldEnergy - energy) < 0.00001) ++count;
		else count = 0;
		if (count >= 1) break;
		oldEnergy = energy;
	}
	std::cout << "\nNo of iterations: " << iter << std::endl;
	fout << "\nNo of iterations: " << iter << std::endl;
	decompile();
	std::cout << ":::Min - Energy = " << energy << "\n:::Time taken = " << timeTaken << " seconds\n";
	fout << ":::Min - Energy = " << energy << "\n:::Time taken = " << timeTaken << " seconds\n";
	fout.flush();

	return (iter < maxiter);
}

template <typename nodeid, typename labelid, typename captype>
typename labelid IRGC<nodeid, labelid, captype>::getLabel(nodeid ii)
{
	return labelling[ii];
}

#include "irgc_instances.inc"
