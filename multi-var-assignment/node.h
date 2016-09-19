#pragma once
#ifndef NODE_H
#define NODE_H

#define MULTIPLE (1.5)

#include <vector>
#include <stack>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include "hash_extend.h"

using namespace std;

typedef vector<string> AssMat; // assignment matrix
typedef vector<int> MintermVector; // cube vector [0, 1, 2, 1]
typedef vector<int> AssVec; // assignment vector [0, 1, 5, 4] ([000 001 101 100])
typedef vector<string> StrAssVec; // string assignment vector ["000", "001", "101", "100"]

// !!! IMPORTANT! the first int is the log2 of the line
// i.e., if a CubeDecomposition has 16 lines,
// then the int should be 4!
typedef pair<int, vector<MintermVector>> CubeDecomposition; 

// use pair instead
//struct cubeDecompostion
//{
//    cubeDecompostion(int L, vector<MintermVector> CV) : line(L), cubeVectors(CV) {}
//
//    int line;
//    vector<MintermVector> cubeVectors;
//};

struct Node
{
    // TODO: constructor to be written
    Node(AssMat newAssMat, MintermVector newProblemVector, int newLevel, unordered_multiset<CubeDecomposition> newAssignedCubeDecompositions, CubeDecomposition lastAssignedCubeDecomposition, int literalCount);

    // assigned "assignment matrix"
    AssMat _assignedAssMat;

    // the remaining problem vector
    MintermVector _remainingProblemVector;

    // the level in the tree; root is _level 0
    int _level;

    // the set of cubes that are already assigned
    unordered_multiset<CubeDecomposition> _assignedCubeDecompositions;

    // the last cube vector assigned
    CubeDecomposition _lastAssignedCubeDecomposition;

    // the literal count so far
    int _literalCountSoFar;
};

class SolutionTree
{
public:

    // methods

    // constructor
    // TODO: complete the constructor
    SolutionTree(vector<int> initialProblemVector, vector<int> degrees, int accuracy, int caseNumber);

    // process the tree
    void ProcessTree();

    // process a single Node
    vector<Node> ProcessNode(Node currentNode);

    // process the node vector, which has 
    vector<Node> ProcessNodeVector(vector<Node> nodeVecToBeProcessed);

    // returns all of the possible cube decompositions
    // INPUT: the log2 of the cube size
    // OUTPUT: the vector of the possible CubeDecomposition
    // NOTE: this method will use _accuracy and the _degrees
    vector<CubeDecomposition> PossibleCubeDecompositions(int log2CubeSize);
    vector<CubeDecomposition> PossibleCubeDecompositionsHelper(int remainingLog2CubeSize, vector<CubeDecomposition> partialDecompositions, vector<int> remainingDegrees);
    vector<MintermVector> PossibleLineCubeVectors(int log2CubeSize, int degree) const;

    // originally "containCube" method
    //static bool CapacityConstraintSatisfied(vector<int> problemVector, MintermVector cubeVector);

    AssMat AssignMatrixByEspresso(AssMat originalAssMat, CubeDecomposition cubeDecompositionToBeAssigned);

    // subtract the cube from the problem vector
    //static MintermVector SubtractCube(MintermVector problemVector, MintermVector cube);




    // data
    int _caseNumber;
    int _lengthOfTotalCube; // (d1+1)(d2+1)...(dn+1)
    vector<int> _degrees;
    int _accuracy;

    vector<Node> _nodeVector;
    int _minLiteralCount;   // current minimum literal count
    Node _optimalNode;
    unordered_map<int, int> _minLiteralCountOfLevel;

    vector<int> _rowGrayCode;
    vector<int> _colGrayCode;

    unordered_map<unordered_multiset<CubeDecomposition>, bool> _processedCubeDecompositions; // true if the set is processed before

    int _updateTime;
    int _nodeNumber;
    int _maxLevel;
    unordered_map<int, int> _sizeOfCubeInLevel;
};

// multiply number and a vector
MintermVector multiply(int line, MintermVector cubeVec);

// multiply two cube vectors
MintermVector multiply(MintermVector cubeVec1, MintermVector cubeVec2);

// multiply n cube vectors, the size of cubeVecs should be greater than 1
MintermVector multiply(vector<MintermVector> cubeVecs);

string IntToBin(int num, int highestDegree);

bool CapacityConstraintSatisfied(vector<int> problemVector, MintermVector cubeVector);

MintermVector SubtractCube(MintermVector problemVector, MintermVector cube);

bool IsZeroMintermVector(MintermVector vec);

long long int choose(int n, int k);

#endif