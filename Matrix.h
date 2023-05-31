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
    Matrix complementOf();
    void complement();
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
    // part 4
    // return true if no configuration of b is in this
    bool avoids(const Matrix& b) const;
    // return true if b is a row permuation of this
    bool row_permutation_of(const Matrix& b) const;
    // part 5
    void shiftUp(unsigned int row);
    void shiftDown(unsigned int row);

  private:
    struct Node {
      Node(vector<int> r, vector<int> s) 
        : rem(r), sel(s) {}
      vector<int> rem, sel;
    };
    unsigned int rows, cols;
    vector<int> elems;
    vector<vector<int>> potentialCols(int row, int k, int pos, const Matrix& b) const;
    bool rowsContain(const Matrix& a, const Matrix& b) const;
    bool avoids(const Node& col_node, const Matrix& b) const;
    bool avoids(const vector<int> &pot_cols, const Node& row_node, const Matrix& b) const;
    bool row_permutation_of(const Node& row_node, const Matrix &b) const;
};

Matrix dot(Matrix a, Matrix b);
int dot(vector<int>& a, vector<int>& b);
Matrix operator*(int a, const Matrix& b);
std::ostream& operator<<(std::ostream& output, const vector<int>& a);
std::ostream& operator<<(std::ostream& output, const vector<vector<int>>& a);
std::ostream& operator<<(std::ostream& output, const Matrix& a);
std::ostream& operator<<(std::ostream& output, const vector<Matrix>& a);
vector<int> makeVec(unsigned int i);
vector<vector<int>> rowColSwap(vector<vector<int>>& a);
template<class E>
bool list_contains(E &elem, vector<E> &list);

// part 6
vector<vector<int>> columns_of_col_sum(unsigned int height, unsigned int sum);
vector<vector<int>> columns_of_K(unsigned int height);
Matrix generate_K(unsigned int height);
Matrix generate_T(unsigned int height);

// part 7
Matrix matrixCross(const Matrix& a, const Matrix& b);
Matrix matrixCombine(const Matrix& a, const Matrix& b);

// part 8
vector<Matrix> Avoid(unsigned int m, const Matrix& F);
unsigned int max_col_count(vector<Matrix>& list);
unsigned int forb(unsigned int m, const Matrix& F);
vector<Matrix> ext(unsigned int m, const Matrix& F);
vector<Matrix> match(unsigned int m, const Matrix& F, unsigned int bound);
vector<Matrix> ext_match_helper(unsigned int bound, vector<Matrix>& list);
vector<Matrix> remove_row_perms(vector<Matrix>& list);
// part 9
void what_is_missing(unsigned int k, Matrix &A);
void remove_duplicates(vector<vector<int>> &columns);
// part 10
unsigned int forbDown(unsigned int m, const Matrix &F, unsigned int start);
vector<Matrix> extDown(unsigned int m, const Matrix &F, unsigned int start);

#endif