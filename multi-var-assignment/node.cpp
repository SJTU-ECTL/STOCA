#include "node.h"

SolutionTree::SolutionTree(vector<int> initialProblemVector, vector<int> degrees, int accuracy, int caseNumber)
{
    _lengthOfTotalCube = static_cast<int>(initialProblemVector.size());
    _degrees = degrees;
    _accuracy = accuracy;
    _caseNumber = caseNumber;
    _minLiteralCount = INT_MAX;

    // initialize root node
    // NOTE: _lengthOfTotalCube equals to (d1+1)(d2+1)...(dn+1)
    int powDegree = static_cast<int>(pow(2, _lengthOfTotalCube));
    int powAccuracy = static_cast<int>(pow(2, _accuracy));
    AssMat assignmentMatrix = AssMat(powAccuracy, string(powDegree, '0'));

    Node root = Node(); // TODO: fill this after complete the constructor of Node

    _nodeVector.push_back(root);

    // initialize Gray code vector TODO

    // initialize other data
    _updateTime = 0;
    _nodeNumber = 1;
    _maxLevel = 0;
}

void SolutionTree::ProcessTree()
{
    int processedNodeNumber = _nodeVector.size();

    while (!_nodeVector.empty())
    {
        // process all of the nodes in the vector
        vector<Node> nodeVectorOfNextLevel = ProcessNodeVector(_nodeVector);

        // update the node vector
        _nodeVector = nodeVectorOfNextLevel;

        // add to the count
        processedNodeNumber += _nodeVector.size();
    }

    // TODO finish the display part
}

vector<Node> SolutionTree::ProcessNode(Node currentNode)
{
    // the vector of nodes to be returned
    // they should have the same literal count
    // and the same size, the size must be the largest size
    vector<Node> ret;

    // currentNode is originally "cNode"

    // at first, we sum up the number of minterms
    //// mintermCount is originally numMinTerm
    int mintermCount = 0;
    for (auto i : currentNode._remainingProblemVector)
    {
        mintermCount += i;
    }

    // then, we find the maximal possible cubes

    bool found = false;

    // log2MintermCount is originally MTNinCube
    // It is the log2 of minterm count in the cube
    int log2MintermCount = static_cast<int>(floor(log2(mintermCount)));

    // the main loop;
    while (!found && (log2MintermCount >= 0))
    {
        // originally, the function here is "possibleCubes"
        // but for multi-var, it should return possible cube 
        // decompositions -- CubeDecomposition, 
        // i.e., vector<int, vector<MintermVector>>
        vector<CubeDecomposition> possibleCubeDecompositionVector = PossibleCubeDecompositions(log2MintermCount);

        // traverse all of the possible cube decompositions

        for (auto cubeDecomposition : possibleCubeDecompositionVector)
        {
            // first we need to get the vector form
            // i.e., we multiply the vectors to get
            // the total vector:
            // (v0, v1, v2)X(u0, u1)=(w00, w10, w20, w01, w11, w21)
            MintermVector totalCubeVector = multiply(static_cast<int>(pow(2, cubeDecomposition.first)), multiply(cubeDecomposition.second));

            // checking the capacity constraint
            if (!CapacityConstraintSatisfied(currentNode._remainingProblemVector, totalCubeVector))
            {
                // if "currentNode._remainingProblemVector" does not contain the
                // cube with the vector of "totalCubeVector"
                // skip this "cubeVector" and process the next one
                continue;
            }

            // "newAssignedCubeDecompositions" is the so called "cube set"
            // it is "new" because the new set will replace the old one
            // Here we declare a new set instead of inserting
            // a vector into the old one because in some cases
            // we will prune this node and keep the old set unchanged
            auto newAssignedCubeDecompositions = currentNode._assignedCubeDecompositions;

            // insert the new cube decomposition
            newAssignedCubeDecompositions.insert(cubeDecomposition);

            // checking for the dupe cube set
            if (_processedCubeDecompositions[newAssignedCubeDecompositions] == true)
            {
                // if this cube set is processed before
                // prune it and display the deletion
                // comment out
                cout << "** BEGIN:Pruning DUPE cube set **" << endl;
                cout << "Following decomposition is duped: " << endl;
                // TODO: display the decomposition
            }

            // set the processed cubeset
            _processedCubeDecompositions[newAssignedCubeDecompositions] = true;

            // assign the cube
            AssMat newAssMat = AssignMatrixByEspresso(currentNode._assignedAssMat, cubeDecomposition);

            // if it's unassignable, go to next cube decomposition
            if (newAssMat[0][0] == 'x')
            {
                continue;
            }

            // evaluate the assigned matrix

            // create a temporary .pla file
            ofstream ofs("temp.pla");
            ofs << ".i " << _lengthOfTotalCube + _accuracy << endl;
            ofs << ".o 1" << endl;
            for (auto line = 0; line < newAssMat.size(); ++line)
            {
                for (auto col = 0; col < newAssMat[line].size(); ++col)
                {
                    if (newAssMat[line][col] == '0') continue;
                    ofs << IntToBin(line, _accuracy - 1) << IntToBin(col, _lengthOfTotalCube - 1) << " 1" << endl;
                }
            }
            ofs << ".e" << endl;

            // invoke the espresso through MVSIS
            system("mvsis50720.exe -c \"read_pla temp.pla; espresso; print_stats -s;\" > tempres.txt");

            ifstream ifs("tempres.txt");

            string tempString;
            int literalCount = 0;
            while (ifs >> tempString)
            {
                if (tempString == "lits")
                {
                    ifs >> tempString >> literalCount;
                    break;
                }
            }

            // prune if the current level's literal count exceeds
            if ((_minLiteralCountOfLevel[currentNode._level] != 0) && (_minLiteralCountOfLevel[currentNode._level] * MULTIPLE <= literalCount))
            {
                // display the pruning message
                // continue;
            }

            // update the current level's literal count
            if ((_minLiteralCountOfLevel[currentNode._level] > literalCount) || (_minLiteralCountOfLevel[currentNode._level] == 0))
            {
                _minLiteralCountOfLevel[currentNode._level] = literalCount;
            }

            // prune if the literal count exceed the overall minimum literal count
            if (literalCount > _minLiteralCount)
            {
                // display the pruning message
                continue;
            }

            // sub node
            MintermVector remainingProblemVector = SubtractCube(currentNode._remainingProblemVector, totalCubeVector);
            Node newNode(newAssMat, remainingProblemVector, currentNode._level + 1, newAssignedCubeDecompositions, cubeDecomposition, literalCount);

            // update the max level
            if (_maxLevel <= currentNode._level)
            {
                _maxLevel = currentNode._level;
            }

            // the current level's cubes' size must be 0 (un-initialized) or <= 2^log2MintermCount
            assert((_sizeOfCubeInLevel[currentNode._level] == 0) || (_sizeOfCubeInLevel[currentNode._level] <= static_cast<int>(pow(2, log2MintermCount))));

            // if it is less than 2^log2MintermCount strictly: reset the following levels
            if (_sizeOfCubeInLevel[currentNode._level] < static_cast<int>(pow(2, log2MintermCount)))
            {
                _sizeOfCubeInLevel[currentNode._level] = static_cast<int>(pow(2, log2MintermCount));
                for (auto i = currentNode._level + 1; i <= _maxLevel; ++i)
                {
                    _sizeOfCubeInLevel[i] = 0;
                }
            }

            found = true;

            // if this is leaf node, update the solution
            if (IsZeroMintermVector(remainingProblemVector))
            {
                // update the solution
                _minLiteralCount = literalCount;
                _optimalNode = newNode;

                // display the update message TODO

                ++_updateTime;

                // continue to skip insertion
                continue;
            }

            // insert the node to the vector
            ret.push_back(newNode);
        } // end of for (auto cubeDecomposition : possibleCubeDecompositionVector)

        --log2MintermCount;
        if ((_sizeOfCubeInLevel[currentNode._level] != 0) && (_sizeOfCubeInLevel[currentNode._level] > static_cast<int>(pow(2, log2MintermCount)))) break;
    } // end of while (!found && (log2MintermCount >= 0))

    return ret;
}

vector<Node> SolutionTree::ProcessNodeVector(vector<Node> nodeVecToBeProcessed)
{
    vector<Node> resultSubNodeVector;

    // first process each node in the vector and get the optimal nodes
    // use for loop to traverse all of the nodes

    for (auto traversedNode : nodeVecToBeProcessed)
    {
        // we need a helper method to process the node
        // INPUT: Node
        // OUTPUT: a vector of sub nodes that have the same literal count and size
        vector<Node> tempNodeVector = ProcessNode(traversedNode);

        // insert the nodes into the vector
        resultSubNodeVector.insert(resultSubNodeVector.end(), tempNodeVector.begin(), tempNodeVector.end());
    }

    // if the level is the leaf _level, then return
    if (resultSubNodeVector.empty()) return resultSubNodeVector;

    auto countSize = [](vector<int> v)
    {
        int s = 0;
        for (auto i : v)
        {
            s += i;
        }
        return s;
    };

    // for each sub node we need to sort them by the size of extracted cube and the literal count
    sort(resultSubNodeVector.begin(), resultSubNodeVector.end(), [countSize](Node n1, Node n2)
    {
        // if the literal count not equal
        if (n1._literalCountSoFar < n2._literalCountSoFar) return true;
        if (n1._literalCountSoFar > n2._literalCountSoFar) return false;

        auto v1 = n1._lastAssignedCubeDecomposition;
        auto v2 = n2._lastAssignedCubeDecomposition;

        int size1 = countSize(multiply(v1.second));
        int size2 = countSize(multiply(v2.second));

        return (size1 < size2);
    });

    int smallestLiteralCount = resultSubNodeVector[0]._literalCountSoFar;
    int smallestCubeSize = countSize(multiply(resultSubNodeVector[0]._lastAssignedCubeDecomposition.second));

    auto delIt = resultSubNodeVector.end();

    for (auto it = resultSubNodeVector.begin(); it != resultSubNodeVector.end(); ++it)
    {
        int size = countSize(multiply((it->_lastAssignedCubeDecomposition).second));
        if ((it->_literalCountSoFar != smallestLiteralCount) || (size != smallestCubeSize))
        {
            delIt = it;
            break;
        }
    }

    // delete all nodes greater than the smallest ones
    resultSubNodeVector.erase(delIt, resultSubNodeVector.end());

    return resultSubNodeVector;
}

vector<CubeDecomposition> SolutionTree::PossibleCubeDecompositions(int log2CubeSize)
{
    return PossibleCubeDecompositionsHelper(log2CubeSize, vector<CubeDecomposition>{make_pair(0, vector<MintermVector>())}, _degrees);
}

vector<CubeDecomposition> SolutionTree::PossibleCubeDecompositionsHelper(int remainingLog2CubeSize, vector<CubeDecomposition> partialDecompositions, vector<int> remainingDegrees)
{
    vector<CubeDecomposition> ret;
    // here remainingDegrees keeps d1, d2, ..., d(k-1)
    if (remainingDegrees.size() == 1)
    {
        // if it is the last variable, use loop to find all decompositions
        // l + s_1 = N'
        for (int s = 0; s < *(remainingDegrees.rbegin()); ++s)
        {
            int l = remainingLog2CubeSize - s;
            if (l > _accuracy)
            {
                // the line number exceeds the hight of the matrix
                continue;
            }
            // if l satisfies the constraint, then find all possible vectors
            // for s_1. We should notice that there are 2^s_1 minterms and the vector
            // is a line cube vector.
            vector<MintermVector> possibleCubes = PossibleLineCubeVectors(s, *remainingDegrees.rbegin());
            for (auto cube : possibleCubes)
            {
                for (auto cubeDecomp : partialDecompositions) // TODO: if partialDecompositions is empty
                {
                    auto newCubeDecomp = cubeDecomp;
                    // assume newCubeDecomp = (0, ([1, 0])), here first = 0 means unset.
                    // assign the log2 line number l, assume it's 2
                    newCubeDecomp.first = l;
                    // then newCubeDecomp = (2, ([1, 0]))
                    // assume cube = [1, 1, 0],
                    newCubeDecomp.second.insert(newCubeDecomp.second.begin(), cube);
                    // now newCubeDecomp = (2, ([1, 1, 0], [1, 0]))
                    // it means 2^2 X [1, 1, 0] X [1, 0] = [4, 4, 0, 0, 0, 0]
                    ret.push_back(newCubeDecomp);
                }
            }
        }

        return ret;
    }

    vector<CubeDecomposition> partialDecompositionsToPass;

    // if the remainingDegrees has more than one term
    // i.e., d1, d2, ..., dk
    // first get the partial decompositions
    for (int s = 0; s < *(remainingDegrees.rbegin()); ++s)
    {
        // if this is the first level of recursion, we assume partialDecompositions.second
        // is empty, i.e., it is (0, ())
        vector<MintermVector> possibleCubes = PossibleLineCubeVectors(s, *remainingDegrees.rbegin());

        if (partialDecompositions.empty())
        {
            for (auto cube : possibleCubes)
            {
                partialDecompositionsToPass.push_back(make_pair(0, vector<MintermVector>{cube}));
            }
        }
        else // partialDecomposition is not empty
        {
            for (auto cube : possibleCubes)
            {
                for (auto cubeDecomp : partialDecompositions)
                {
                    auto newCubeDecomp = cubeDecomp;
                    newCubeDecomp.second.insert(newCubeDecomp.second.begin(), cube);
                    partialDecompositionsToPass.push_back(newCubeDecomp);
                }
            }
        }

        auto remainingDegreesToPass = remainingDegrees;
        remainingDegreesToPass.pop_back();

        auto results = PossibleCubeDecompositionsHelper(remainingLog2CubeSize - s, partialDecompositionsToPass, remainingDegreesToPass);
        ret.insert(ret.end(), results.begin(), results.end());
    }

    return ret;
}

vector<MintermVector> SolutionTree::PossibleLineCubeVectors(int log2CubeSize, int degree) const
{
    vector<MintermVector> ret;

    auto powAcc = static_cast<int>(pow(2, _accuracy));

    for (auto zeroBefore = 0; zeroBefore <= degree - log2CubeSize; ++zeroBefore)
    {
        MintermVector line;

        for (auto i = 0; i < zeroBefore; ++i)
        {
            line.push_back(0);
        }

        for (auto i = 0; i <= log2CubeSize; ++i)
        {
            long long int mintermCountAtPositionI = choose(log2CubeSize, i);
            line.push_back(mintermCountAtPositionI);
        }

        for (auto i = 0; i < degree - log2CubeSize - zeroBefore; ++i)
        {
            line.push_back(0);
        }

        ret.push_back(line);
    }
    return ret;
}

AssMat SolutionTree::AssignMatrixByEspresso(AssMat originalAssMat, CubeDecomposition cubeDecompositionToBeAssigned)
{
    // TODO: finish this
    auto powAcc = static_cast<int>(pow(2, _accuracy));

    // input example 2^1 [1,1,0] [1,1]
    // we could use std::set to impliment assignment vector
    // or we should call it assignment column set
    // first, for this example, we need to get
    // the assignment set for [1,1,0]
    // which is {00,01} and {00,10}
    // and the assignment set for [1,1]
    // which is {0,1}

    // vector<MintermVector> vectors
    // for each MintermVector in the cubeDecompositionToBeAssigned
    // we find its assignment set
    
    auto assignmentSetsVector = vector<vector<set<string>>>();

    // find all assignment sets for each cube vector
    for (auto mintermVec : cubeDecompositionToBeAssigned.second)
    {
        // for cube minterm vector [1, 1, 0], we find the assignment sets as following
        // {00, 01} and {00, 10}.
        // vector<set<int>> assignmentSets
        auto assignmentSets = FindAssignmentSetsOfStringForMintermVector(mintermVec);
        assignmentSetsVector.push_back(assignmentSets);
    }

    // it contains all of the possible assignment sets
    // for this cube decomposition in string
    auto assignmentSetsOfStringForCubeDecomposition = FindAssignmentSetsOfStringForCubeDecomposition(assignmentSetsVector);

    // transform the string version to integer version
    auto assignmentSetsOfIntForCubeDecomposition = vector<set<int>>();
    for (auto setOfString : assignmentSetsOfStringForCubeDecomposition)
    {
        auto setOfInt = set<int>();
        for (auto str : setOfString)
        {
            setOfInt.insert(BinToInt(str));
        }
        assignmentSetsOfIntForCubeDecomposition.push_back(setOfInt);
    }

    return AssMat();
}

vector<set<string>> SolutionTree::FindAssignmentSetsOfStringForMintermVector(MintermVector lineCubeVector) const
{
    auto countOfOne = 0;

    // find the first non-zero term in lineCubeVector
    // because the number of zeros in the beginning
    // is the number of un-complemented X-variables
    for (; countOfOne < lineCubeVector.size(); ++countOfOne)
    {
        if (lineCubeVector[countOfOne] == 0) continue;
        break;
    }

    auto countOfMinterm = 0;

    // find the count of minterms
    for (auto i = countOfOne; i < lineCubeVector.size(); ++i)
    {
        countOfMinterm += lineCubeVector[i];
    }

    // the number of zeros in the end is the 
    // number of complemented X-variables
    // it's d - t - r, d = size - 1, t = countOfOne, r = log2(countOfMinterm)
    // since countOfMinterm = choose(r, 0) + choose(r, 1) + ... + choose(r, r)
    auto countOfZero = lineCubeVector.size() - 1 - countOfOne - static_cast<int>(log2(countOfMinterm));

    // then we have all of the information we need

    // next we need to find a basic assignment set according to the 
    // count of minterms. For example, [1, 1, 0] has two minterms,
    // the basic assignment set will be {0, 1}, then insert 0/1
    // according to the countOfOne and countOfZero we have got.
    auto basicAssignmentSet = BuildBasicAssignmentSet(countOfMinterm);

    // then we will find all of the assignment sets
    auto assignmentSetsOfString = BuildAssignmentSet(basicAssignmentSet, countOfZero, countOfOne);
    return assignmentSetsOfString;
}

set<string> SolutionTree::BuildBasicAssignmentSet(int mintermCount) const
{
    auto ret = set<string>();

    // if the cube has only one minterm, like [0, 1, 0]
    // its basic assignment set is defined as {""}
    if (mintermCount == 1)
    {
        ret.insert("");
        return ret;
    }

    auto log2MintermCount = static_cast<int>(log2(mintermCount));
    auto grayCode = ConstructGrayCode(log2MintermCount);

    for (auto i : grayCode)
    {
        ret.insert(IntToBin(i, log2MintermCount - 1));
    }

    return ret;
}

vector<set<string>> SolutionTree::BuildAssignmentSet(set<string> basicAssignmentSet, int countOfZero, int countOfOne)
{
    auto ret = vector<set<string>>();
    auto zeroOneTwoPermutation = BuildZeroOneTwoPermutation(countOfZero, basicAssignmentSet.begin()->size(), countOfOne);

    for (auto pattern : zeroOneTwoPermutation)
    {
        auto tempAssignmentSet = set<string>();

        for (auto str : basicAssignmentSet)
        {
            auto pos = 0;
            auto res = string();

            for (auto i = 0; i < pattern.size(); ++i)
            {
                if (pattern[i] == '0')
                {
                    res += "0";
                }
                else if (pattern[i] == '2')
                {
                    res += "1";
                }
                else
                {
                    res += str[pos++];
                }
            }
            tempAssignmentSet.insert(res);
        }
        ret.push_back(tempAssignmentSet);
    }
    return ret;
}

vector<set<string>> SolutionTree::FindAssignmentSetsOfStringForCubeDecomposition(vector<vector<set<string>>> assignmentSetsVector) const
{
    return FindAssignmentSetsOfStringForCubeDecompositionHelper(assignmentSetsVector, vector<set<string>>{set<string>{""}});
}

vector<set<string>> SolutionTree::FindAssignmentSetsOfStringForCubeDecompositionHelper(vector<vector<set<string>>> remainingAssignmentSetsVector, vector<set<string>> currentAssignmentSets) const
{
    assert(!remainingAssignmentSetsVector.empty());

    auto setVec = vector<set<string>>();
    for (auto set1 : currentAssignmentSets)
    {
        for (auto set2 : remainingAssignmentSetsVector[0])
        {
            setVec.push_back(MultiplyAssignmentSets(set1, set2));
        }
    }

    if (remainingAssignmentSetsVector.size() == 1)
    {
        return setVec;
    }

    // remainingAssignmentSetsVector.size() >= 2
    remainingAssignmentSetsVector.erase(remainingAssignmentSetsVector.begin());
    return FindAssignmentSetsOfStringForCubeDecompositionHelper(remainingAssignmentSetsVector, setVec);
}

MintermVector multiply(int line, MintermVector cubeVec)
{
    for (auto i : cubeVec)
    {
        i *= line;
    }
    return cubeVec;
}

MintermVector multiply(MintermVector cubeVec1, MintermVector cubeVec2)
{
    MintermVector product;
    for (auto u : cubeVec2)
    {
        for (auto v : cubeVec1)
        {
            product.push_back(v * u);
        }
    }
    return product;
}

MintermVector multiply(vector<MintermVector> cubeVecs)
{
    assert(!cubeVecs.empty());
    if (cubeVecs.size() == 1)
    {
        return cubeVecs[0];
    }

    auto cv = multiply(cubeVecs[0], cubeVecs[1]);

    if (cubeVecs.size() == 2)
    {
        return cv;
    }

    // size >= 3
    cubeVecs.erase(cubeVecs.begin());
    cubeVecs[0] = cv;

    return multiply(cubeVecs);
}

string IntToBin(int num, int highestDegree)
{
    assert(num < static_cast<int>(pow(2, highestDegree + 1)));

    string ret;

    while (num)
    {
        char digit = static_cast<char>((num & 1) + '0');
        ret.insert(ret.begin(), 1, digit);

        num >>= 1;
    }

    ret.insert(ret.begin(), highestDegree + 1 - ret.size(), '0');

    return ret;
}

int BinToInt(string str)
{
    auto ret = 0;
    for (auto ch : str)
    {
        ret <<= 1;
        ret |= static_cast<int>(ch - '0');
    }
    return ret;
}

bool CapacityConstraintSatisfied(vector<int> problemVector, MintermVector cubeVector)
{
    assert(problemVector.size() == cubeVector.size());

    for (auto i = 0; i < problemVector.size(); ++i)
    {
        if (problemVector[i] < cubeVector[i]) return false;
    }
    return true;
}

MintermVector SubtractCube(MintermVector problemVector, MintermVector cube)
{
    assert(CapacityConstraintSatisfied(problemVector, cube));
    for (auto i = 0; i < problemVector.size(); ++i)
    {
        problemVector[i] -= cube[i];
    }
    return problemVector;
}

bool IsZeroMintermVector(MintermVector vec)
{
    for (auto i : vec)
    {
        if (i != 0) return false;
    }
    return true;
}

long long int choose(int n, int k)
{
    if ((n < k) || (k < 0)) return 0;

    long long int ret = 1;

    for (auto i = 1; i <= k; ++i)
    {
        ret *= n--;
        ret /= i;
    }

    return ret;
}

vector<int> ConstructGrayCode(int size)
{
    auto grayCode = vector<int>();
    grayCode.push_back(0);
    grayCode.push_back(1);

    for (auto i = 1; i < size; ++i)
    {
        grayCode = ConstructGrayCodeHelper(grayCode);
    }

    return grayCode;
}

vector<int> ConstructGrayCodeHelper(vector<int> grayCode)
{
    auto ret = vector<int>();
    for (auto i = 0; i < grayCode.size(); ++i)
    {
        if (i % 2 == 0)
        {
            ret.push_back(grayCode[i] << 1);
            ret.push_back((grayCode[i] << 1) | 1);
        }
        else
        {
            ret.push_back((grayCode[i] << 1) | 1);
            ret.push_back(grayCode[i] << 1);
        }
    }
    return ret;
}

vector<string> BuildZeroOneTwoPermutation(int countOfZero, int countOfOne, int countOfTwo)
{
    auto vic = string(countOfZero, '0') + string(countOfOne, '1') + string(countOfTwo, '2');
    auto ret = vector<string>();

    do
    {
        ret.push_back(vic);
    } while (next_permutation(vic.begin(), vic.end()));

    return ret;
}

set<string> MultiplyAssignmentSets(set<string> set1, set<string> set2)
{
    auto ret = set<string>();
    for (auto str1 : set1)
    {
        for (auto str2 : set2)
        {
            ret.insert(str1 + str2);
        }
    }
    return ret;
}

set<string> MultiplyAssignmentSets(vector<set<string>> sets)
{
    assert(!sets.empty());
    if (sets.size() == 1)
    {
        return sets[0];
    }

    auto prod = MultiplyAssignmentSets(sets[0], sets[1]);
    if (sets.size() == 2)
    {
        return prod;
    }

    // if the size >= 3
    sets.erase(sets.begin());
    sets[0] = prod;

    return MultiplyAssignmentSets(sets);
}

Node::Node(AssMat newAssMat, MintermVector newProblemVector, int newLevel, unordered_multiset<CubeDecomposition> newAssignedCubeDecompositions, CubeDecomposition lastAssignedCubeDecomposition, int literalCount)
{
    _assignedAssMat = newAssMat;
    _remainingProblemVector = newProblemVector;
    _level = newLevel;
    _assignedCubeDecompositions = newAssignedCubeDecompositions;
    _lastAssignedCubeDecomposition = lastAssignedCubeDecomposition;
    _literalCountSoFar = literalCount;
}
