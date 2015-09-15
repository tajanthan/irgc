
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "IRGC.h"

using namespace std;

typedef float captype;
typedef long nodeid;
typedef int labelid;

static nodeid width;
static nodeid height;
static labelid labels;
static nodeid numNodes;
static labelid lambda;
static captype *unaryPotentialArray;
static captype *binaryWeightArray;
static captype *binaryFunctionArray;

void errorFunction(char* msg)
{
	cout << endl << "Error!\n " << msg << endl;
	exit(1);
}

void readUnaryPotential(char* unaryF)
{
	ifstream fs(unaryF);
	assert(fs);

	char c;
	nodeid ii;
	labelid l;
	captype val;
	nodeid count = 0;
	while (fs >> ii >> c >> l >> c >> val) {	// read comma separated values		
		unaryPotentialArray[ii * labels + l] = val;
		++count;
	}
	assert(count == numNodes * labels);
	fs.close();
}

void readBinaryWeights(char* binaryWF)
{
	ifstream fs(binaryWF);
	assert(fs);

	nodeid numEdges = width * (height - 1) + height * (width - 1);

	char c;
	nodeid i, j;
	captype val;
	nodeid count = 0;
	while (fs >> i >> c >> j >> c >> val) {	// read comma separated values
		if (j == i + 1) binaryWeightArray[i] = val;	// right neighbor
		if (j == i + width) binaryWeightArray[i + numNodes] = val;	// down neighbor
		++count;
	}
	assert(count == numEdges);
	fs.close();
}

void cauchy()
{
	for (labelid li = 0; li < labels; ++li) {
		for (labelid lj = 0; lj < labels; ++lj) {
			labelid x = std::abs(li - lj);
			binaryFunctionArray[li * labels + lj] = (captype)std::pow(lambda, 2) * log(1 + std::pow(((float)x / lambda), 2)) / 2;
		}
	}
}

void truncatedLinear()
{
	for (labelid li = 0; li < labels; ++li) {
		for (labelid lj = 0; lj < labels; ++lj) {
			labelid x = std::abs(li - lj);
			binaryFunctionArray[li * labels + lj] = (captype)std::min(x, lambda);
		}
	}
}

void truncatedQuadratic()
{
	for (labelid li = 0; li < labels; ++li) {
		for (labelid lj = 0; lj < labels; ++lj) {
			labelid x = std::abs(li - lj);
			binaryFunctionArray[li * labels + lj] = (captype)std::pow(std::min(x, lambda), 2);
		}
	}
}

void populateBinaryFunction(char* binary) {
	if (!binary) errorFunction("Invalid binary function type!");
	else if (!strcmp(binary, (char*)"TL")) truncatedLinear();	// truncated linear
	else if (!strcmp(binary, (char*)"TQ")) truncatedQuadratic();	// truncated quadratic
	else if (!strcmp(binary, (char*)"CA")) cauchy();	// cauchy
	else errorFunction("Invalid binary function type!");
}

captype unaryPotential(nodeid i, labelid li)
{
	return unaryPotentialArray[i * labels + li];
}

captype binaryWeights(nodeid i, nodeid j)
{
	if (i + 1 == j) return binaryWeightArray[i];	// right
	else if (i + width == j) return binaryWeightArray[i + numNodes];	// down
	else if (j + 1 == i) return binaryWeightArray[j];	// left
	else if (j + width == i) return binaryWeightArray[j + numNodes];	// up
	else assert(0);
	return 0;
}

captype binaryFunction(labelid li, labelid lj)
{
	return binaryFunctionArray[li * labels + lj];
}

captype binaryPotential(nodeid i, nodeid j, labelid li, labelid lj)
{
	return binaryWeights(i, j) * binaryFunction(li, lj);
}
 
void createLocalArrays(nodeid w, nodeid h, labelid l, char* unaryF, char* binaryWF, char* binary, labelid lam)
{
	width = w;
	height = h;
	labels = l;
	lambda = lam;
	numNodes = width * height;
	unaryPotentialArray = new captype[numNodes*labels]();
	binaryWeightArray = new captype[numNodes * 2]();	// right, down neighbor
	binaryFunctionArray = new captype[labels * labels]();

	readUnaryPotential(unaryF);
	readBinaryWeights(binaryWF);
	populateBinaryFunction(binary);	
}

void deleteLocalArrays()
{
	delete[] unaryPotentialArray;
	delete[] binaryWeightArray;
	delete[] binaryFunctionArray;
}

int main(int argc, char **argv)
{
	char* usage = " [width] [height] [labels] [unary potential file (V*L)] [binary weights file (E)] \
		[binary potential function (TL, TQ, CA)] [lambda (inflection point)] [hybrid] [max-iter]";
	if (argc < 8) {
		cout << "Error!" << endl;
		cout << "Usage: " << argv[0] << usage << endl;
		return 1;
	}
	createLocalArrays((nodeid)atoi(argv[1]), (nodeid)atoi(argv[2]), (labelid)atoi(argv[3]), argv[4], argv[5], argv[6], (labelid)atoi(argv[7]));

	bool hybrid = false;
	if (argc > 8) if (atoi(argv[8])) hybrid = true;
	size_t maxiter = 100;
	if (argc > 9) {
		int m = atoi(argv[9]);
		if (m > 0) maxiter = (size_t)m;
	}

	IRGC<nodeid, labelid, captype> irgc(width, height, labels, hybrid);
	irgc.setEnergy(unaryPotential, binaryWeights, binaryFunction, binaryPotential, argv[6], (labelid)atoi(argv[7]));
	irgc.optimize(maxiter);

	deleteLocalArrays();

	return 0;
}




