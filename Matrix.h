#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <ostream>

using std::vector;

class Matrix
{
  public:
    // part 1
    Matrix(int r, int c, vector<int> nums);
    Matrix(int ver, vector<vector<int>> nums);
    void printMatrix();
    vector<vector<int>> getRows() const;
    vector<vector<int>> getCols() const;
    vector<int> getElems() const;
    unsigned int numRows() const;
    unsigned int numCols() const;
    Matrix dotWith(Matrix b);
    Matrix transposeOf();
    void transpose();
    // part 2
    void multiply(int b);
    void multiply(const Matrix& b);
    void add(const Matrix& b);
    Matrix operator+(const Matrix& b);    
    Matrix operator*(const Matrix& b);
    Matrix operator+=(const Matrix& b);
    bool operator==(const Matrix& b) const;
    Matrix operator*(int a);
    // part 3 (aka the real stuff)
    // return true if b is a submatrix of this (this contains b)
    bool contains(const Matrix& b) const;
    // return true if this is a submatrix of b (b contains this)
    bool containedIn(const Matrix& b) const;
    // no repeated columns
    bool isSimple();
    // return true if no configuration of b is in this
    bool avoid(const Matrix& b) const;
    

  public:
    struct Node {
      Node(vector<int> r, vector<int> s) 
        : rem(r), sel(s) {}
      vector<int> rem, sel;
    };
    unsigned int rows, cols;
    vector<int> elems;
    vector<vector<int>> potentialCols(int row, int k, int pos, const Matrix& b) const;
    bool rowsContain(const Matrix& a, const Matrix& b) const;
    bool avoid(const Node& col_node, const Matrix& b) const;
    bool avoid(const vector<int> &pot_cols, const Node& row_node, const Matrix& b) const;
};

Matrix dot(Matrix a, Matrix b);
int dot(vector<int>& a, vector<int>& b);
Matrix operator*(int a, const Matrix& b);
std::ostream& operator<<(std::ostream& output, const vector<int>& a);
std::ostream& operator<<(std::ostream& output, const vector<vector<int>>& a);
std::ostream& operator<<(std::ostream& output, const Matrix& a);
vector<int> makeVec(unsigned int i);
vector<vector<int>> rowColSwap(vector<vector<int>>& a);

#endif