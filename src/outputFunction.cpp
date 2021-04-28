#include <iostream>
#include <stdio.h>
#include <vector>

#include <string>
#include <sstream>
#include <fstream>

#include <unordered_set>

using namespace std;

void writeTo(vector<int> &partElements, vector<int> &eind, vector<vector<float> > &points){
    std::ofstream outfile ("../output/finalOutput.txt", ios_base::app);
    cout << "from here" << endl;
    for(int n=0;n<partElements.size();n++) {
        outfile << points[eind[partElements[n]*3]][0] << "," <<points[eind[partElements[n]*3]][1] << endl;
        outfile << points[eind[partElements[n]*3+1]][0] << "," <<points[eind[partElements[n]*3+1]][1] << endl;
        outfile << points[eind[partElements[n]*3+2]][0] << "," <<points[eind[partElements[n]*3+2]][1] << endl;
    }
}

int main() {
    int nProcesses = 4;
    for(int i=0;i<nProcesses;i++){
        std::ifstream infile("../output/test"+std::to_string(i) +".txt");
        std::string line;  
        std::getline(infile, line);
        std::istringstream iss(line);

        int a, m ,n;
        if (!(iss >> a)) {}
        vector<int> partElements;
        while (a > 0)
        {
            std::getline(infile, line);
            std::istringstream issT(line);
            float b;
            if (!(issT >> b)) { break; }
            partElements.push_back(b);
            a--;
        }
        std::getline(infile, line);
        std::istringstream iss2(line);
        if (!(iss2 >> m)) {}
        vector<int> eind;
        while (m > 0)
        {
            std::getline(infile, line);
            std::istringstream iss2T(line);
            float b;
            if (!(iss2T >> b)) { break; }
            eind.push_back(b);
            m--;
        }
        std::getline(infile, line);
        std::istringstream iss3(line);
        if (!(iss3 >> n)) {}
        cout << n;
        vector<vector<float> > points;
        while (n > 0)
        {
            std::getline(infile, line);
            std::istringstream iss3T(line);
            float b,c;
            if (!(iss3T >> b >> c)) { break; }
            vector<float> v;
            v.push_back(b);
            v.push_back(c);
            points.push_back(v);
            n--;
        }
        writeTo(partElements, eind, points);
    }
    return 0;
}