#include "node.h"

using namespace std;

int main(int argc, char* argv[])
{
    // first line: number of variables, n and accuracy m
    // second line: degrees
    // following lines: one problem vector per line

    int variableNumber;
    int accuracy;
    vector<int> degrees;
    int temp;

    ifstream ifs("vector.txt");
    auto line = string();
    auto caseCount = 0;
    getline(ifs, line); // first line
    stringstream ss;
    ss << line;
    ss >> variableNumber >> accuracy;
    ss.str(string());
    ss.clear();
    getline(ifs, line); // second line
    ss << line;

    while (ss >> temp)
    {
        degrees.push_back(temp);
    }

    assert(degrees.size() == variableNumber);
    

    cout << "Variable number: " << variableNumber << endl;
    cout << "Accuracy: " << accuracy << endl;
    cout << "Degrees: ";
    for (auto d : degrees)
    {
        cout << d << " ";
    }
    cout << endl << endl;

    while (getline(ifs, line))
    {
        cout << "Case No. " << caseCount << endl;
        cout << "======" << endl;

        ss.str(string());
        ss.clear();
        ss << line;
        MintermVector problemVector;
        while (ss >> temp)
        {
            problemVector.push_back(temp);
        }

        auto solutionTree = SolutionTree(problemVector, degrees, accuracy, caseCount);
        solutionTree.ProcessTree();
        ++caseCount;
    }

     return 0;
}
