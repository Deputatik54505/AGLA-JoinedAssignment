//Ivan Platonov
//i.platonov@innopolis.university
#include <iostream>
#include <vector>
#include<cmath>

using namespace std;
#define tableType double
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

class SquareMatrix;

class ColumnVector;

//Base matrix class
class Matrix {
protected:
    int rows;
    int cols;
    vector<vector<tableType>> table;

//function for creating table of Identity matrix
    void createIdentity(int n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                table[i][j] = i == j ? 1 : 0;
            }
        }
    }

//function for printing matrix table
    static void printTable(const vector<vector<tableType>> &a) {
        cout << Matrix(a);
    }

    static void printTable(const vector<vector<tableType>> &a, const vector<vector<tableType>> &b) {
        vector<vector<tableType>> tab;
        tab.resize(a.size());
        for (int i = 0; i < a.size(); i++) {
            tab[i].resize(a[0].size() + b[0].size());
        }
        for (int i = 0; i < b.size(); i++) {
            for (int j = 0; j < a.size(); j++) {
                tab[i][j] = a[i][j];
            }
            for (int j = 0; j < b.size(); j++) {
                tab[i][a[0].size() + j] = b[i][j];
            }
        }
        cout << Matrix(tab);
    }

public:
    Matrix(int r, int c) {
        rows = r;
        cols = c;
        table.resize(rows);
        for (int i = 0; i < rows; i++) {
            table[i].resize(cols);
        }
    }

    explicit Matrix(const vector<vector<tableType>> &a) {
        rows = (int) a.size();
        cols = (int) a[0].size();
        table = a;
    }

    explicit operator SquareMatrix() const;

    explicit operator ColumnVector() const;

    friend istream &operator>>(istream &is, Matrix &m) {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.cols; j++) {
                is >> m.table[i][j];
            }
        }
        return is;
    }

    friend ostream &operator<<(ostream &os, const Matrix &m) {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.cols; j++) {

                printf("%.4f", m.table[i][j]);

                if (j != m.cols - 1) {
                    os << " ";
                }
            }
            os << endl;
        }
        return os;
    }

    Matrix &operator=(const Matrix &m) {
        rows = m.rows;
        cols = m.cols;
        table.resize(rows);
        for (int i = 0; i < rows; i++) {
            table[i].resize(cols);
            for (int j = 0; j < cols; j++) {
                table[i][j] = m.table[i][j];
            }
        }
        return *this;
    }

    Matrix operator+(const Matrix &m) const {
        if (rows != m.rows || cols != m.cols) {
            throw exception();
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.table[i][j] = table[i][j] + m.table[i][j];
            }
        }
        return result;
    }

    Matrix operator-(const Matrix &m) const {
        if (rows != m.rows || cols != m.cols) {
            throw exception();
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.table[i][j] = table[i][j] - m.table[i][j];
            }
        }
        return result;
    }

    Matrix operator*(const Matrix &m) const {
        if (cols != m.rows) {
            throw exception();
        }
        Matrix result(rows, m.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < m.cols; j++) {
                result.table[i][j] = 0;
                for (int k = 0; k < cols; k++) {
                    result.table[i][j] += table[i][k] * m.table[k][j];
                }
            }
        }
        return result;
    }

    Matrix getDiagonal() {
        vector<vector<tableType>> diagTab;
        diagTab.resize(rows);

        for (int i = 0; i < min(rows, cols); i++) {
            diagTab[i].resize(cols);
            diagTab[i][i] = table[i][i];
        }
        Matrix diag(diagTab);
        return diag;
    }

    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                result.table[i][j] = this->table[j][i];
            }
        }
        return result;
    }

    bool isDiagonalDominant() {
        tableType sum = 0;
        for (int i = 0; i < rows; i++) {
            sum = 0;
            for (int j = 0; j < cols; j++) {
                if (j != i) {
                    sum += abs(table[i][j]);
                }
            }
            if (sum > abs(table[i][i]))
                return false;
        }
        return true;
    }

    Matrix getLowerPart() {
        vector<vector<tableType>> lowTab;
        lowTab.resize(rows);

        for (int i = 0; i < min(rows, cols); i++) {
            for (int j = 0; j <= i; j++) {
                lowTab[i].resize(cols);
                lowTab[i][j] = table[i][j];
            }
        }
        Matrix low(lowTab);
        return Matrix(lowTab);
    }


    int getRowsNum() const {
        return rows;
    }

    int getColsNum() const {
        return cols;
    }

    vector<vector<tableType>> getTable() const {
        return table;
    };
};

class ColumnVector : public Matrix {
public:
    explicit ColumnVector(int n) : Matrix(n, 1) {};

    explicit ColumnVector(const vector<vector<tableType>> &a) : Matrix(a) {
        if (this->cols != 1) {
            throw exception();
        }
    };
};

Matrix::operator ColumnVector() const {
    return ColumnVector(this->getTable());
}

class SquareMatrix : public Matrix {
private:
    void pivot(vector<vector<tableType>> &a, vector<vector<tableType>> &b, int row) const {
        int max_index = row;
        tableType max_value = abs(a[row][row]);

        for (int i = row + 1; i < n; i += 1) {
            if (abs(a[i][row]) > max_value) {
                max_index = i;
                max_value = abs(a[i][row]);
            }
        }
        if (max_index != row) {
            swap(a[max_index], a[row]);
            swap(b[max_index], b[row]);
            printTable(a, b);
        }
    }

    void solveLin(vector<vector<tableType>> &a, vector<vector<tableType>> &b) const {
        int i = 0;
        //eliminateColumn matrix to upper triangular form
        for (; i < n; i++) {
            eliminateColumn(a, b, i, 1);
        }
        i--;
        //eliminateColumn matrix to getDiagonal form
        for (; i >= 0; i--) {
            eliminateColumn(a, b, i, -1);
        }
        i = 0;
        for (; i < n; i++) {
            for (int j = 0; j < b[0].size(); j++) {
                b[i][j] = b[i][j] / a[i][i];
            }
            a[i][i] = 1;
        }
    }

    void eliminateColumn(vector<vector<tableType>> &a, vector<vector<tableType>> &b, int column, int step) const {
        for (int i = column + step; i < n && i >= 0; i += step) {
            bool flag = true;

            tableType factor = a[i][column] / a[column][column];
            flag = factor != 0;

            for (int j = 0; flag && j < n && j >= 0; j++) {
                a[i][j] -= factor * a[column][j];
                b[i][j] -= factor * b[column][j];
            }
        }
    }

protected:
    int n;

public:
    explicit SquareMatrix(const vector<vector<tableType>> &a) : Matrix(a) {
        if (a.size() != a[0].size()) {
            throw exception();
        }
        this->n = this->rows;
    }


    explicit SquareMatrix(int n) : Matrix(n, n) { this->n = n; };

    SquareMatrix inverse() {
        vector<vector<tableType>> a = this->getTable();
        SquareMatrix identity(this->n);
        identity.createIdentity(this->n);
        vector<vector<tableType>> b = identity.getTable();
        solveLin(a, b);
        return SquareMatrix(b);
    }

    tableType JacobiIteration(const ColumnVector &b, double error) {
        if (!isDiagonalDominant()) {
            cout << "The method is not applicable!" << endl;
            return -1;
        }
        SquareMatrix diag = SquareMatrix(getDiagonal());
        diag = diag.inverse();
        SquareMatrix identity(n);
        identity.createIdentity(n);
        SquareMatrix alpha = SquareMatrix(identity - (diag * *this));
        ColumnVector betta = ColumnVector(diag * b);

        cout << "alpha:" << endl;
        cout << alpha;

        cout << "beta:" << endl;
        cout << betta;

        ColumnVector x = betta;

        double currentError;

        cout << "x(" << 0 << "):" << endl;
        cout << x;

        for (int iteration = 1;; iteration++) {


            ColumnVector adjustment = ColumnVector(diag * (b - *this * x));
            x = ColumnVector(x + adjustment);
            currentError = 0;
            for (int i = 0; i < x.getRowsNum(); i++) {
                currentError += pow(adjustment.getTable()[i][0], 2);
            }
            currentError = pow(currentError, 0.5);
            printf("e: %.4f \n", currentError);
            cout << "x(" << iteration << "):" << endl;
            cout << x;
            if (currentError <= error) {
                break;
            }
        }

        return currentError;

    }

    tableType Seidel(const ColumnVector &b, double error) {
        if (!isDiagonalDominant()) {
            cout << "The method is not applicable!" << endl;
            return -1;
        }
        SquareMatrix diag = SquareMatrix(getDiagonal());
        diag = diag.inverse();
        SquareMatrix identity(n);
        identity.createIdentity(n);
        SquareMatrix alpha = SquareMatrix(identity - (diag * *this));
        ColumnVector betta = ColumnVector(diag * b);
        SquareMatrix B = SquareMatrix(alpha.getLowerPart());
        SquareMatrix C = SquareMatrix(alpha - B);


        cout << "beta:" << endl;
        cout << betta;

        cout << "alpha:" << endl;
        cout << alpha;

        cout << "B:" << endl;
        cout << B;

        cout << "C:" << endl;
        cout << C;

        B = SquareMatrix(identity - B);
        cout << "I-B:" << endl;
        cout << B;

        B = B.inverse();
        cout << "(I-B)_-1:" << endl;
        cout << B;

        ColumnVector x = betta;
        double currentError;

        cout << "x(" << 0 << "):" << endl;
        cout << x;

        B = SquareMatrix(this->getLowerPart());
        C = SquareMatrix(*this - B);
        B = B.inverse();
        for (int iteration = 1;; iteration++) {

            ColumnVector prevX = x;
            x = ColumnVector(B * (b - (C * prevX)));
            currentError = 0;
            for (int i = 0; i < x.getRowsNum(); i++) {
                currentError += pow(x.getTable()[i][0] - prevX.getTable()[i][0], 2);
            }
            currentError = pow(currentError, 0.5);
            printf("e: %.4f \n", currentError);
            cout << "x(" << iteration << "):" << endl;
            cout << x;
            if (currentError <= error) {
                break;
            }
        }

        return currentError;


    }

};

Matrix::operator SquareMatrix() const {
    return SquareMatrix(this->getTable());
}


class IdentityMatrix : public SquareMatrix {
public:
    explicit IdentityMatrix(int n) : SquareMatrix(n) {
        createIdentity(n);
    };
};


class PermutationMatrix : public SquareMatrix {
public:
    explicit PermutationMatrix(int n) : SquareMatrix(n) {
        createIdentity(n);
    };

    void changeRows(int first, int second) {
        table[first][first] = 0;
        table[second][second] = 0;
        table[first][second] = 1;
        table[second][first] = 1;
    }
};


int main() {
    FILE *pipe = _popen(GNUPLOT_NAME, "w");

    if (pipe == nullptr) {
        cout << "Could not open pipe" << endl;
        exit(0);
    }
    int experimentsCount;
    cin >> experimentsCount;
    fprintf(pipe, "%s\n", "plot [-3:3][0:5] '-' with points, '-' with lines");

    vector<double> xValues;
    vector<vector<double>> yValues;
    xValues.resize(experimentsCount);
    yValues.resize(experimentsCount);

    for (int i = 0; i < experimentsCount; i++) {
        cin >> xValues[i];
        yValues[i].resize(1);
        cin >> yValues[i][0];

        fprintf(pipe, "%f\t%f\n", xValues[i], yValues[i][0]);
    }
    fprintf(pipe, "%s\n", "e");

    int polynomialDegree;
    cin >> polynomialDegree;

    vector<vector<double>> aTable;
    aTable.resize(experimentsCount);

    for (int i = 0; i < experimentsCount; i++) {
        aTable[i].resize(polynomialDegree + 1);
        for (int j = 0; j <= polynomialDegree; j++) {
            aTable[i][j] = pow(xValues[i], j);
        }
    }


    Matrix a(aTable);
    cout << "A:" << endl;
    cout << a;
    SquareMatrix ata(a.transpose() * a);
    cout << "A_T*A:" << endl;
    cout << ata;
    SquareMatrix inv = ata.inverse();
    cout << "(A_T*A)^-1:" << endl;
    cout << inv;
    ColumnVector b(yValues);


    Matrix rightSight = a.transpose() * b;
    cout << "A_T*b:" << endl;
    cout << rightSight;
    Matrix res = inv * rightSight;
    cout << "x~:" << endl;
    cout << res;

    fprintf(pipe, "%f\t%f\n", 0.0, res.getTable()[0][0]);
    tableType yValue = 0;
    for (auto a: res.getTable()) {
        yValue += a[0];
    }
    fprintf(pipe, "%f\t%f\n", 1.0, yValue);

    fprintf(pipe, "%s\n", "e");

    fflush(pipe);

}



