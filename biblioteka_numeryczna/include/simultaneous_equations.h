#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip>

using namespace std;

const double EPSILON = 1e-9;

void loadData(const string& filename, int& N, vector<double>& b, vector<vector<double>>& A); // Wczytywanie danych

void printMatrix(const vector<vector<double>>& matrix, const string& name, int iteration = -1); // Wyświtlanie kolejnych kroków algorytmu

void luDecomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, int N); // Rozwiązaywanie układu równań

vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& b, int N); // Podstawienie do przodu

vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y, int N); // Podstawienie do tyłu

vector<double> solveLU(vector<vector<double>>& A, vector<double>& b, int N); // Tutaj się rozwiązuje układ
