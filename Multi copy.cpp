
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <string>

#include <iterator>
#include <algorithm>
#include <array>
#include <random>
#include <chrono>

#include "include/vector2.h"
#include "include/triangle.h"
#include "include/delaunay.h"

using namespace std;

ifstream inFile ("input.txt");

vector<vector<double>> vect;
vector<double> col;

int main(int argc, char** argv)
{
    int x;
    inFile >> x;
    cout << x << endl;
    for(int i=0;i<x;i++) {
        int x1,x2;
        inFile >> x1 >> x2;
        col.push_back(x1);
        col.push_back(x2);
        vect.push_back(col);
        col.clear();
    }

	std::vector<dt::Vector2<double>> points;
	for(int i = 0; i < x; ++i) {
        points.push_back(dt::Vector2<double>{4.33, 2.13});
	}

    dt::Delaunay<double> triangulation;
	const auto start = std::chrono::high_resolution_clock::now();
	const std::vector<dt::Triangle<double>> triangles = triangulation.triangulate(points);
	const auto end = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double> diff = end - start;

    return 0;
}
