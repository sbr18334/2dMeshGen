#include <iostream>
#include <stdio.h>
#include <vector>

#include <string>
#include <sstream>
#include <fstream>

#include <unordered_set>

using namespace std;

int main() {
    vector<vector<float> > points;
    std::ifstream infile("../input/A.1.node");
    std::string line;
    std::getline(infile, line);
    std::istringstream iss(line);
    int a_c;
    if (!(iss >> a_c)) {}
    while (a_c > 0)
    {
        std::getline(infile, line);
        std::istringstream issT(line);
        float a2, b, c, d;
        if (!(issT >> a2 >> b >> c >> d)) { break; }
        vector<float> v;
        v.push_back(b);
        v.push_back(c);
        points.push_back(v);
        a_c--;
    }

    vector<vector<int> > elements;
    std::ifstream infile2("../input/A.1.ele");
    std::string line2;
    std::getline(infile2, line2);
    std::istringstream iss2(line2);
    int a;
    if (!(iss2 >> a)) {}
    cout << a << endl;
    while (a > 0)
    {
        std::getline(infile2, line2);
        std::istringstream iss2T(line2);
        int a2, b, c, d;
        if (!(iss2T >> a2 >> b >> c >> d)) { break; }
        cout << a2 << b << c << d << endl;
        vector<int> v1;
        v1.push_back(b);
        v1.push_back(c);
        v1.push_back(d);
        elements.push_back(v1);
        a--;
    }

    std::ofstream outfile ("../output/finalInput.txt", ios_base::app);

    for(int i=0;i<elements.size();i++) {
        outfile << points[elements[i][0]-1][0] << "," <<points[elements[i][0]-1][1] << endl;
        outfile << points[elements[i][1]-1][0] << "," <<points[elements[i][1]-1][1] << endl;
        outfile << points[elements[i][2]-1][0] << "," <<points[elements[i][2]-1][1] << endl;
    }
    return 0;
}