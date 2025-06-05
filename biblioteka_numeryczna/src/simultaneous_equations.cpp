#include "../include/simultaneous_equations.h"

void loadData(const string& filename, int& N, vector<double>& b, vector<vector<double>>& A) 
{
    ifstream file(filename);
    if (!file) {
        cerr << "Nie mozna otworzyc pliku!" << endl;
        exit(1);
    }
    string line;
    while (getline(file, line)) {
        if (line.find("N = ") != string::npos) {
            stringstream ss(line.substr(4));
            ss >> N;
            break;
        }
    }
    while (getline(file, line)) {
        if (line.find("b:") != string::npos) {
            getline(file, line);
            stringstream ss(line);
            double val;
            while (ss >> val) {
                b.push_back(val);
            }
            break;
        }
    }
    while (getline(file, line)) {
        if (line.find("A:") != string::npos) {
            for (int i = 0; i < N; ++i) {
                getline(file, line);
                stringstream ss(line);
                vector<double> row;
                double val;
                while (ss >> val) {
                    row.push_back(val);
                }
                A.push_back(row);
            }
            break;
        }
    }
    file.close();
}

void printMatrix(const vector<vector<double>>& matrix, const string& name, int iteration = -1) 
{
    cout << "\n" << name;
    if (iteration != -1) cout << " (po iteracji " << iteration << ")";
    cout << ":\n";
    for (const auto& row : matrix) {
        for (auto val : row)
            cout << setw(8) << setprecision(5) << val;
        cout << endl;
    }
    cout << endl;
}

void luDecomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, int N) 
{
    L.assign(N, vector<double>(N, 0.0));
    U.assign(N, vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            U[i][j] = A[i][j];
            for (int k = 0; k < i; ++k)
                U[i][j] -= L[i][k] * U[k][j];
        }

        for (int j = i; j < N; ++j) {
            if (fabs(U[i][i]) < EPSILON) {
                cerr << "Blad: Macierz osobliwa, nie mozna zapisac macierzy A w postaci\nmacierzy dolnej trojkatnej L i macierzy gornej trojkatnej U." << endl;
                exit(1);
            }
            if (i == j)
                L[i][i] = 1;
            else {
                L[j][i] = A[j][i];
                for (int k = 0; k < i; ++k)
                    L[j][i] -= L[j][k] * U[k][i];
                L[j][i] /= U[i][i];
            }
        }
        printMatrix(L, "Macierz L (po iteracji " + to_string(i) + ")");
        printMatrix(U, "Macierz U (po iteracji " + to_string(i) + ")");
        cout << "\na11 = u11: " << A[0][0] << " = " << U[0][0] << endl;
    }
}

vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& b, int N) 
{
    vector<double> y(N, 0.0);
    for (int i = 0; i < N; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j)
            y[i] -= L[i][j] * y[j];
    }
    cout << "\nWektor z:" << endl;
    for (double val : y) {
        cout << setw(8) << setprecision(4) << val << " ";
    }
    cout << endl;
    return y;
}

vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y, int N) 
{
    vector<double> x(N, 0.0);
    for (int i = N - 1; i >= 0; --i) {
        if (fabs(U[i][i]) < EPSILON) {
            cerr << "Błąd: Macierz osobliwa, dzielenie przez zero." << endl;
            exit(1);
        }
        x[i] = y[i];
        for (int j = i + 1; j < N; ++j)
            x[i] -= U[i][j] * x[j];
        x[i] /= U[i][i];
    }
    cout << "\nWektor x:" << endl;
    for (double val : x) {
        cout << setw(8) << setprecision(4) << val << " ";
    }
    cout << endl;
    return x;
}

vector<double> solveLU(vector<vector<double>>& A, vector<double>& b, int N) 
{
    vector<vector<double>> L, U;
    luDecomposition(A, L, U, N);
    vector<double> y = forwardSubstitution(L, b, N);
    vector<double> x = backwardSubstitution(U, y, N);

    // Sprawdzenie poprawności rozwiązania
    vector<double> Ax(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Ax[i] += A[i][j] * x[j];
        }
    }

    vector<vector<double>> LU(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                LU[i][j] += L[i][k] * U[k][j];
            }
        }
    }
    cout << "\nSprawdzenie poprawnosci rozwiazania (L * U powinno byc rowne A):\n";
    printMatrix(LU, "Iloczyn L * U");
    printMatrix(A, "Macierz A");

    cout << "\nSprawdzenie poprawnosci rozwiazania (Ax powinno byc rowne b):\n";
    for (int i = 0; i < N; ++i) {
        double diff = fabs(Ax[i] - b[i]);
        cout << "Element " << i + 1 << ": Ax = " << setw(8) << setprecision(2) << Ax[i]
            << ", b = " << setw(8) << setprecision(2) << b[i]
            << ", Roznica = " << setw(8) << setprecision(2) << diff << endl;
    }

    return {};
}