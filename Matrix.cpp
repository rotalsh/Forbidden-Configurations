#include "Matrix.h"
#include <iostream>
#include <algorithm>
#include <ostream>
#include <vector>

using std::vector;

/**
 * Returns the number of rows in the matrix
 */
unsigned int Matrix::numRows() const {
  return rows;
}

/**
 * Returns the number of columns in the matrix
 */
unsigned int Matrix::numCols() const {
  return cols;
}

/**
 * Constructor for Matrix
 * @param r the number of rows of the matrix
 * @param c the number of columns in the matrix
 * @param nums the elements of the matrix - the order in nums is such that each row is filled first
 * REQUIRE: r*c == nums.size()
 */
Matrix::Matrix(int r, int c, vector<int> nums) {
  // checking size
  if ((unsigned int) r * (unsigned int) c != nums.size()) {
    std::cout << "Error, wrong size" << std::endl;
    return;
  }
  rows = r;
  cols = c;
  elems = nums;
}

/**
 * Constructor for Matrix
 * @param ver 0 means nums is vector of rows, 1 means nums is vector of columns
 * @param nums elements of the matrix, represented either as a vector of rows or vector of columns
 */
Matrix::Matrix(int ver, vector<vector<int>> nums) {
  if (nums.empty()) return;

  if (ver == 0) { // vector of rows
    rows = nums.size();
    cols = nums[0].size();
    for (unsigned int i = 0; i < nums.size(); i++) {
      for (unsigned int j = 0; j < nums[0].size(); j++) {
        elems.push_back(nums[i][j]);
      }
    }
  } else { // vector of columns
    cols = nums.size();
    rows = nums[0].size();
    for (unsigned int i = 0; i < nums[0].size(); i++) {
      for (unsigned int j = 0; j < nums.size(); j++) {
        elems.push_back(nums[j][i]);
      }
    }
  }
}

/**
 * Return the result of matrix multiplication
 * @param b the matrix that comes second in the multiplication
 */
Matrix Matrix::dotWith(Matrix b) {
  return dot(*this, b);
}

/**
 * Print contents of the matrix
 */
void Matrix::printMatrix() {
  std::cout << "_" << std::endl;
  if (elems.empty()) return;
  for (unsigned int i = 0; i < rows; i++) {
    std::cout << "| ";
    for (unsigned int j = 0; j < cols; j++) {
      std::cout << elems[i*cols + j] << " ";
    }
    std::cout << "|" << std::endl;
  }
}

/**
 * Return the elements in the matrix as a vector of rows
 */
vector<vector<int>> Matrix::getRows() const {
  vector<vector<int>> ret;
  for (unsigned int i = 0; i < rows; i++) {
    vector<int> row;
    for (unsigned int j = 0; j < cols; j++) {
      row.push_back(elems[i*cols + j]);
    }
    ret.push_back(row);
  }
  return ret;
}

/**
 * Return the elements in the matrix as a vector of columns
 */
vector<vector<int>> Matrix::getCols() const {
  vector<vector<int>> ret;
  for (unsigned int i = 0; i < cols; i++) {
    vector<int> col;
    for (unsigned int j = 0; j < rows; j++) {
      col.push_back(elems[i + j*cols]);
    }
    ret.push_back(col);
  }
  return ret;
}

/**
 * Return the elements of the matrix as a single vector of elements
 */
vector<int> Matrix::getElems() const
{
  return elems;
}

/**
 * Return the transpose of the matrix
 */
Matrix Matrix::transposeOf() {
  return Matrix(0, getCols());
}

/**
 * Transpose the matrix (rotate along the main diagonal)
 */
void Matrix::transpose() {
  Matrix A(0, getCols());
  rows = A.rows;
  cols = A.cols;
  elems = A.elems;
}

/**
 * Return the 0-1 complement of the matrix
 */
Matrix Matrix::complementOf() {
  vector<int> comp;
  for (unsigned int i = 0; i < elems.size(); i++) {
    comp.push_back(elems[i] ? 0 : 1);
  }
  return Matrix(rows, cols, comp);
}

/**
 * 0-1 complement the matrix (set 0's to 1's and vice versa)
 */
void Matrix::complement() {
  Matrix other = complementOf();
  elems = other.elems;
}

/**
 * Multiply each element in the matrix by an integer a
 * @param a the integer to multiply
 */
void Matrix::multiply(int a) {
  for (unsigned int i = 0; i < elems.size(); i++) {
    elems[i] *= a;
  }
}

/**
 * Perform elementwise multiplication of this matrix and b
 * @param b the matrix to elementwise multiply to
 * REQUIRE: this->rows == b.rows && this->cols == b.cols
 */
void Matrix::multiply(const Matrix &b){
  if (rows != b.rows || cols != b.cols) {
    std::cout << "The matrices must have the same size\n";
    return;
  }
  for (unsigned int i = 0; i < elems.size(); i++) {
    elems[i] *= b.elems[i];
  }
}

/**
 * Perform elementwise addition of this matrix and b
 * @param b the matrix to elementwise add to
 * REQUIRE: this->rows == b.rows && this->cols == b.cols
 */
void Matrix::add(const Matrix& b) {
  if (rows != b.rows || cols != b.cols) {
    std::cout << "The matrices must have the same size\n";
    return;
  }
  for (unsigned int i = 0; i < elems.size(); i++) {
    elems[i] += b.elems[i];
  }
}

/**
 * Return the result of performing elementwise addition on this matrix and b
 * @param b the matrix to elementwise add to
 * REQUIRE: this->rows == b.rows && this->cols == b.cols
 */
Matrix Matrix::operator+(const Matrix &b)
{
  if (rows != b.rows || cols != b.cols) {
    std::cout << "The matrices must have the same size\n";
    return *this;
  }
  vector<int> nums;
  for (unsigned int i = 0; i < elems.size(); i++) {
    nums.push_back(elems[i] + b.elems[i]);
  }
  return Matrix(rows, cols, nums);
}

/**
 * Return the result of performing elementwise multiplication on this matrix and b
 * the matrix to elementwise multiply to
 * REQUIRE: this->rows == b.rows && this->cols == b.cols
 */
Matrix Matrix::operator*(const Matrix &b)
{
  if (rows != b.rows || cols != b.cols) {
    std::cout << "The matrices must have the same size\n";
    return *this;
  }
  vector<int> nums;
  for (unsigned int i = 0; i < elems.size(); i++) {
    nums.push_back(elems[i] * b.elems[i]);
  }
  return Matrix(rows, cols, nums);
}

/**
 * Return the result of multiplying each element in the matrix by an integer a
 * @param a the integer to multiply
 */
Matrix Matrix::operator*(int a)
{
  vector<int> nums;
  for (unsigned int i = 0; i < elems.size(); i++) {
    nums.push_back(elems[i] * a);
  }
  return Matrix(rows, cols, nums);
}

/**
 * Return true if this matrix has a submatrix that is equal to b
 * @param b the matrix that may or may not be a submatrix
 */
bool Matrix::contains(const Matrix &b) const
{
  if (rows < b.rows || cols < b.cols) return false;
  if (*this == b) return true;
  for (unsigned int i = 0; i < rows; i++) {
    vector<vector<int>> potCols = potentialCols(i, 0, 0, b);
    if (potCols.empty()) continue;
    for (unsigned int j = 0; j < potCols.size(); j++) {
      vector<vector<int>> selectedCols;
      vector<vector<int>> thisCols = getCols();
      for (unsigned int k = 0; k < potCols[j].size(); k++) {
        selectedCols.push_back(thisCols[potCols[j][k]]);
      }
      Matrix A(1, selectedCols);
      if (rowsContain(A, b)) 
        return true;
    }
  }
  return false;
}

/**
 * Return true if b has a submatrix that is equal to this matrix
 * @param b the matrix that may or may not contain this matrix as a submatrix
 */
bool Matrix::containedIn(const Matrix &b) const
{
  return b.contains(*this);
}

/**
 * Helper function for contains - special case where both matrices have the same column length, only
 *  need to check row equality in this case
 * @param a larger matrix that b could be in
 * @param b smaller matrix that a could contain
 */
bool Matrix::rowsContain(const Matrix &a, const Matrix &b) const {
  int pos = 0;
  vector<vector<int>> aRows = a.getRows();
  vector<vector<int>> bRows = b.getRows();
  for (unsigned int i = 0; i < a.rows; i++) {
    if (aRows[i] == bRows[pos]) {
      if ((unsigned int) ++pos >= b.rows) 
        return true;
    }
  }
  return false;
}

/**
 * Return true if this matrix does not contain b as a configuration
 * @param b the matrix that may or may not be a configuration
 */
bool Matrix::avoids(const Matrix &b) const
{
  if (rows < b.rows || cols < b.cols) return true;
  if (*this == b) return false;
  return avoids(Node(makeVec(b.cols), vector<int>()), b);
}

/**
 * Private helper function for avoid, checking for different columns first
 * Check different column permutations of b to see if this matrix has that column permutation as a submatrix
 * @param col_node the permutation of columns to be chosen
 * @param b the matrix that many or may not be a configuration
 */
bool Matrix::avoids(const Node& col_node, const Matrix &b) const
{
  if (col_node.rem.empty()) {
    return avoids(col_node.sel, Node(makeVec(b.rows), vector<int>()), b);
  }
  for (unsigned int i = 0; i < col_node.rem.size(); i++) {
    vector<int> curr_sel = col_node.sel;
    vector<int> curr_rem = col_node.rem;
    curr_sel.push_back(curr_rem[i]);
    curr_rem.erase(curr_rem.begin() + i);    
    bool attempt = avoids(Node(curr_rem, curr_sel), b);
    if (!attempt) return false;
  }
  return true;
}

/**
 * Private helper function for avoid, given a fixed column check different rows
 * Check different row permutations of the given column permutations in pot_cols
 *  to see if this matrix has that row and column permutation as a submatrix
 * @param pot_cols the permutation of columns being checked currently
 * @param row_node the permutation of rows to be chosen
 * @param the matrix that may or may not be a configuration
 */
bool Matrix::avoids(const vector<int> &pot_cols, const Node &row_node, const Matrix &b) const
{
  if (row_node.rem.empty()) {
    vector<vector<int>> temp;
    for (unsigned int i = 0; i < pot_cols.size(); i++) {
      temp.push_back(b.getCols()[pot_cols[i]]);
    }
    vector<vector<int>> mat;
    temp = rowColSwap(temp);
    for (unsigned int j = 0; j < row_node.sel.size(); j++) {
      mat.push_back(temp[row_node.sel[j]]);
    }
    return !contains(Matrix(0, mat));
  }
  for (unsigned int i = 0; i < row_node.rem.size(); i++) {
    vector<int> curr_sel = row_node.sel;
    vector<int> curr_rem = row_node.rem;
    curr_sel.push_back(curr_rem[i]);
    curr_rem.erase(curr_rem.begin() + i);    
    bool attempt = avoids(pot_cols, Node(curr_rem, curr_sel), b);
    if (!attempt) return false;
  }
  return true;
}

/**
 * Return a vector of vector of ints, where each vector of ints represents the columns of this
 *  you can select so that on given row, entries in the selected columns match the entries of 
 *  the first row of b for entries pos onwards
 * @param row the row of this matrix to match with
 * @param k the leftmost column to search from
 * @param pos the current position to match with in matrix b
 * @param b the matrix of which to match the first row with 
 */
vector<vector<int>> Matrix::potentialCols(int row, int k, int pos, const Matrix& b) const {
  vector<vector<int>> ret;
  for (unsigned int i = k; i < cols; i++) {
    if (elems[row*cols + i] == b.elems[pos]) {
      if ((unsigned int) pos+1 >= b.cols) {
        ret.push_back({(int) i});
      } else {
        vector<vector<int>> afters = potentialCols(row, i+1, pos+1, b);
        for (unsigned int j = 0; j < afters.size(); j++) {
          afters[j].emplace(afters[j].begin(), i);
          ret.push_back(afters[j]);
        }
      }
    }
  }
  return ret;
}

/**
 * Return true if b is a row permutation of this matrix
 * @param b the matrix that may or may not be a row permutation
 */
bool Matrix::row_permutation_of(const Matrix &b) const
{
  if (cols != b.cols || rows != b.rows) 
    return false;
  return row_permutation_of(Node(makeVec(b.rows), vector<int>()), b);
}

/**
 * Private helper function for row_permutation_of
 * Return true if b is a row permutation of this matrix
 * @param row_node a node keeping track of which rows have been added so far
 * @param b the matrix that may or may not be a row permutation
 */
bool Matrix::row_permutation_of(const Node& row_node, const Matrix &b) const {
  if (row_node.rem.empty()) {
    vector<vector<int>> mat; 
    vector<vector<int>> temp = b.getRows();
    for (unsigned int j = 0; j < row_node.sel.size(); j++) {
      mat.push_back(temp[row_node.sel[j]]);
    }
    return *this == (Matrix(0, mat));
  }
  for (unsigned int i = 0; i < row_node.rem.size(); i++) {
    vector<int> curr_sel = row_node.sel;
    vector<int> curr_rem = row_node.rem;
    curr_sel.push_back(curr_rem[i]);
    curr_rem.erase(curr_rem.begin() + i);    
    bool attempt = row_permutation_of(Node(curr_rem, curr_sel), b);
    if (attempt) return true;
  }
  return false;
}

/**
 * Function for shifting down (subtract 1 from given row for each column unless it makes a copy)
 * 0 is preserved as 0
 * @param row the row to shift (0-indexed)
 * REQUIRES: matrix is simple, row < this->rows
 */
void Matrix::shiftDown(unsigned int row)
{
  if (row >= rows || !this->isSimple()) {
    std::cout << "Shifting undefined for this case" << std::endl;
  }
  vector<vector<int>> toShift = getCols();
  for (unsigned int i = 0; i < toShift.size(); i++) {
    vector<int> temp = toShift[i];
    temp[row] = 0;
    if (!list_contains(temp, toShift)) {
      toShift[i] = temp;
    }
  }
  Matrix other(1, toShift);
  elems = other.elems;
}

/**
 * Function for shifting up (add 1 from given row for each column unless it makes a copy)
 * 1 is preserved as 1
 * @param row the row to shift (0-indexed)
 * REQUIRES: matrix is simple, row < this->rows
 */
void Matrix::shiftUp(unsigned int row)
{
  if (row >= rows || !this->isSimple()) {
    std::cout << "Shifting undefined for this case" << std::endl;
  }
  vector<vector<int>> toShift = getCols();
  for (unsigned int i = 0; i < toShift.size(); i++) {
    vector<int> temp = toShift[i];
    temp[row] = 1;
    if (!list_contains(temp, toShift)) {
      toShift[i] = temp;
    }
  }
  Matrix other(1, toShift);
  elems = other.elems;
}

/**
 * Return true if the matrix has no repeated columns
 */
bool Matrix::isSimple()
{
  vector<vector<int>> c = getCols();
  for (unsigned int i = 0; i < c.size(); i++) {
    for (unsigned int j = i + 1; j < c.size(); j++) {
      if (c[i] == c[j]) return false;
    }
  }
  return true;
}

/**
 * Add the matrix b to this matrix then return the result of the assignment
 * @param b the matrix to add
 */
Matrix Matrix::operator+=(const Matrix &b)
{
  return *this = *this + b;
}

/**
 * Return true if this matrix and b have the same number of rows, columns, and the same elements
 * @param b the matrix to check equality
 */
bool Matrix::operator==(const Matrix &b) const
{
  if (this == &b) return true;
  return rows == b.rows && cols == b.cols && elems == b.elems;
}

/**
 * Return the result of doing matrix multiplication on two matrices
 * @param a the first matrix in the multiplication
 * @param b the second matrix in the multiplication
 * REQUIRE: a.numCols() == b.numRows()
 */
Matrix dot(Matrix a, Matrix b) {
  if (a.numCols() != b.numRows()) {
    std::cout << "The number of columns in the first matrix have to equal the number of rows in the second\n";
    return a;
  }
  int r = a.numRows();
  int c = b.numCols();
  vector<int> nums;
  vector<vector<int>> rows = a.getRows();
  vector<vector<int>> cols = b.getCols();
  for (unsigned int i = 0; i < rows.size(); i++) {
    for (unsigned int j = 0; j < cols.size(); j++) {
      nums.push_back(dot(rows[i], cols[j]));
    }
  }
  Matrix A(r, c, nums);
  return A;
}

/**
 * Return the result of multiplying each element in b by an integer a
 * @param a the integer to multiply
 * @param b the matrix to multiply
 */
Matrix operator*(int a, const Matrix &b)
{
  vector<int> nums;
  for (unsigned int i = 0; i < b.getElems().size(); i++) {
    nums.push_back(b.getElems()[i] * a);
  }
  return Matrix(b.numRows(), b.numCols(), nums);
}

/**
 * Definition of the ostream operator<< for a vector of ints
 * @param output the output to add to
 * @param a the vector of ints
 */
std::ostream &operator<<(std::ostream &output, const vector<int> &a)
{
  for (unsigned int i = 0; i < a.size(); i++) {
    output << a[i] << (i == a.size()-1 ? "" : " ");
  }
  return output;
}

/**
 * Definition of the ostream operator<< for a vector of vector of ints
 * @param output the output to add to
 * @param a the vector of vector of ints
 */
std::ostream& operator<<(std::ostream& output, const vector<vector<int>>& a) {
  for (unsigned int i = 0; i < a.size(); i++) {
    output << a[i] << (i == a.size()-1 ? "" : "\n");
  }
  return output;
}

/**
 * Definition of the ostream operator<< for a matrix
 * @param output the output to add to
 * @param a the matrix
 */
std::ostream &operator<<(std::ostream &output, const Matrix &a)
{
  output << "_" << std::endl;
  unsigned int rows = a.numRows(), cols = a.numCols();
  vector<int> elems = a.getElems();
  if (elems.empty()) {
    return output;
  }
  for (unsigned int i = 0; i < rows; i++) {
    output << "| ";
    for (unsigned int j = 0; j < cols; j++) {
      output << elems[i*cols + j] << " ";
    }
    output << "|" << (i == rows-1 ? "" : "\n");
  }
  return output;
}

/**
 * Definition of the ostream operator<< for a vector of matrices
 * @param output the output to add to
 * @param a the vector of matrices
 */
std::ostream& operator<<(std::ostream& output, const vector<Matrix>& a) {
  for (unsigned int i = 0; i < a.size(); i++) {
    output << a[i] << (i == a.size()-1 ? "" : "\n");
  }
  return output;
}

/**
 * Make a vector of size i containing the elements [0, ..., i - 1] in that order
 * @param i the number of elements the output vector should have
 */
vector<int> makeVec(unsigned int i)
{
  vector<int> ret;
  for (unsigned int k = 0; k < i; k++) {
    ret.push_back(k);
  }
  return ret;
}

/**
 * Assume a represents a vector of rows of a matrix.
 *  Return the vector of columns of this same matrix (or vice versa - assume vector of columns, return rows).
 * @param a the vector of vector of ints
 * REQUIRE: a[i].size() is the same for all i in [0, a.size()-1]
 */
vector<vector<int>> rowColSwap(vector<vector<int>> &a)
{
  vector<vector<int>> ret;
  for (unsigned int i = 0; i < a[0].size(); i++) {
    vector<int> vec;
    for (unsigned int j = 0; j < a.size(); j++) {
      vec.push_back(a[j][i]);
    }
    ret.push_back(vec);
  }
  return ret;
}

/**
 * Return a vector of vector of ints, each whose size is equal to height and have the same number of 1's as sum.
 * @param height the size of each vector
 * @param sum the number of 1's each vector should have
 * REQUIRES: sum <= height
 */
vector<vector<int>> columns_of_col_sum(unsigned int height, unsigned int sum)
{
  vector<vector<int>> ret;
  if (sum == 0) {
    vector<int> zeros;
    for (unsigned int i = 0; i < height; i++) {
      zeros.push_back(0);
    }
    ret.push_back(zeros);
    return ret;
  }
  if (height == sum) {
    vector<int> ones;
    for (unsigned int i = 0; i < height; i++) {
      ones.push_back(1);
    }
    ret.push_back(ones);
    return ret;
  }
  ret = columns_of_col_sum(height - 1, sum);
  for (unsigned int i = 0; i < ret.size(); i++) {
    ret[i].emplace(ret[i].begin(), 0);
  }
  vector<vector<int>> add_ones = columns_of_col_sum(height - 1, sum - 1);
  for (unsigned int j = 0; j < add_ones.size(); j++) {
    add_ones[j].emplace(add_ones[j].begin(), 1);
  }
  add_ones.insert(add_ones.end(), ret.begin(), ret.end());
  return add_ones;
}

/**
 * Return K_height - a height * 2^height matrix of all possible columns of size height
 * @param height the height the generated matrix should have
 */
Matrix generate_K(unsigned int height)
{
  return Matrix(1, columns_of_K(height));
}

/**
 * Return the vector of ints that form the columns of the matrix K_height
 * @param height the size of each vector
 */
vector<vector<int>> columns_of_K(unsigned int height)
{
  vector<vector<int>> all_columns;
  for (unsigned int i = 0; i <= height; i++) {
    vector<vector<int>> sum_columns = columns_of_col_sum(height, i);
    all_columns.insert(all_columns.end(), sum_columns.begin(), sum_columns.end());
  }
  return all_columns;
}

/**
 * Return T_height - an upper triangular matrix of size height
 * @param height the number of rows in the returned matrix
 * REQUIRES: height > 0
 */
Matrix generate_T(unsigned int height)
{
  vector<vector<int>> all_columns;
  for (unsigned int i = 1; i <= height; i++) {
    vector<int> one_col;
    for (unsigned int j = 0; j < i; j++) {
      one_col.push_back(1);
    }
    for (unsigned int k = 0; k < height - 1; k++) {
      one_col.push_back(0);
    }
    all_columns.push_back(one_col);
  }
  return Matrix(1, all_columns);
}

/**
 * Return the result of performing the dot product between two vectors (sum of elementwise multiplcation)
 * @param a the first vector
 * @param b the second vector
 * REQUIRE: a.size() == b.size()
 */
int dot(vector<int>& a, vector<int>& b) {
  if (a.size() != b.size()) {
    std::cout << "Vectors need to be of same length to dot product.\n";
    return 0;
  }
  int ret = 0;
  for (unsigned int i = 0; i < a.size(); i++) {
    ret += a[i] * b[i];
  }
  return ret;
}

/**
 * Return true if list contains elem
 * @param elem the element to check for
 * @param list the list to check within
 */
template<class E>
bool list_contains(E &elem, vector<E> &list) {
  for (E e : list) {
    if (e == elem) return true;
  }
  return false;
}

/**
 * Return the matrix cross product between two matrices - 
 *  Let m and n be the number of columns in a and b, respectively. For each column in a, make n new columns 
 *  which are the elements of the column in a followed by the elements of the n columns in b.
 *  The result is a matrix with m*n columns.
 * @param a the first matrix
 * @param b the second matrix
 */
Matrix matrixCross(const Matrix& a, const Matrix& b) {
  vector<vector<int>> a_cols = a.getCols();
  vector<vector<int>> b_cols = b.getCols();
  vector<vector<int>> combined;
  for (unsigned int i = 0; i < a_cols.size(); i++) {
    vector<int> curr_col = a_cols[i];
    for (unsigned int j = 0; j < b_cols.size(); j++) {
      vector<int> col = curr_col;
      col.insert(col.end(), b_cols[j].begin(), b_cols[j].end());
      combined.push_back(col);
    }
  }
  return Matrix(1, combined);
}

/**
 * Given two matrices return a new matrix whose columns are the columns of a followed by the columns of b
 * @param a the first matrix
 * @param b the second matrix
 * REQUIRED: a and b have the same number of rows
 */
Matrix matrixCombine(const Matrix& a, const Matrix& b) {
  vector<vector<int>> a_cols = a.getCols();
  vector<vector<int>> b_cols = b.getCols();
  a_cols.insert(a_cols.end(), b_cols.begin(), b_cols.end());
  return Matrix(1, a_cols);
}

/**
 * Private helper for Avoid - generative recursion
 * Pop the last element in sel. Make a copy of rem. Add that last element to one of the copies of rem (call it rem_added).
 *  If rem_added avoids F, perform recursion on rem_added and sel, then add the matrix made from rem_added to the results.
 *  If not, do nothing with rem_added.
 *  Next, perform recursion on rem_skipped (the other copy of rem without the extra element), add these matrices
 *  to the vector of matrices from the previous step, then return the vector.
 * @param rem a vector of columns to choose from
 * @param sel a vector of columns already chosen
 * @param F the matrix to avoid
 */
vector<Matrix> Avoid(vector<vector<int>> rem, vector<vector<int>> sel, const Matrix &F) {
  if (rem.empty()) 
    return vector<Matrix>();
  vector<vector<int>> skip = sel;
  sel.push_back(rem.back());
  rem.pop_back();
  Matrix test(1, sel);
  vector<Matrix> ret;
  if (test.avoids(F)) {
    ret = Avoid(rem, sel, F);
    ret.push_back(test);
  }
  vector<Matrix> other = Avoid(rem, skip, F);
  ret.insert(ret.end(), other.begin(), other.end());
  return ret;
}

/**
 * Return a vector of m-rows, simple matrices such that for each vector A in the matrix, 
 *  A.avoids(F) is true.
 * @param m the number of rows each matrix should have
 * @param F the configuration to avoid
 */
vector<Matrix> Avoid(unsigned int m, const Matrix &F)
{
  vector<vector<int>> cols = columns_of_K(m);
  std::reverse(cols.begin(), cols.end());
  return Avoid(cols, vector<vector<int>>(), F);
}

/**
 * Given a vector of matrices, return the max of the number of columns from the matrices in the vector
 * @param list the vector of matrices
 */
unsigned int max_col_count(vector<Matrix> &list)
{
  unsigned int max = 0;
  for (unsigned int i = 0; i < list.size(); i++) {
    max = std::max(max, list[i].numCols());
  }
  return max;
}

/**
 * Return the maximum number of columns a matrix can have if it has m rows and avoids F.
 * @param m the number of rows each matrix should have
 * @param F the configuration to avoid
 */
unsigned int forb(unsigned int m, const Matrix &F)
{
  vector<Matrix> list = Avoid(m, F);
  return max_col_count(list);
}

/**
 * Return a vector of matrices that are in Avoid(m, F) which have the same number of columns as forb(m, F)
 * @param m the number of rows each matrix should have
 * @param F the configuration to avoid
 */
vector<Matrix> ext(unsigned int m, const Matrix &F)
{
  vector<Matrix> list = Avoid(m, F);
  unsigned int bound = max_col_count(list);
  return ext_match_helper(bound, list);
}

/**
 * Return a vector of matrices that are in Avoid(m, F) which have the same number of columns as bound
 * @param m the number of rows each matrix should have
 * @param F the configuration to avoid
 * @param bound the number of columns each matrix should have
 */
vector<Matrix> match(unsigned int m, const Matrix &F, unsigned int bound)
{
  vector<Matrix> list = Avoid(m, F);
  return ext_match_helper(bound, list);
}

/**
 * Return a subset of the list of matrices that have the same number of columns as given bound
 * @param bound the number of columns each matrix should have
 * @param list the vector of matrices
 */
vector<Matrix> ext_match_helper(unsigned int bound, vector<Matrix> &list)
{
  vector<Matrix> ret;
  for (unsigned int i = 0; i < list.size(); i++) {
    Matrix curr = list[i];
    if (curr.numCols() == bound) 
      ret.push_back(curr);
  }
  return ret;
}

/**
 * Given a vector of matrices, return a new vector of matrices such that none of the elements in the
 *  new vector are row permutations of one another
 * @param list the vector of matrices
 */
vector<Matrix> remove_row_perms(vector<Matrix> &list)
{
  vector<Matrix> ret;
  for (unsigned int i = 0; i < list.size(); i++) {
    Matrix curr = list[i];
    int outer = 0;
    for (unsigned int j = 0; j < ret.size(); j++) {
      if (curr.row_permutation_of(ret[j])) {
        outer = 1;
        break;
      }
    }
    if (outer) continue;
    ret.push_back(curr);
  }
  return ret;
}

/**
 * On each possible combination of k rows in A, tells you which column of K_k is missing
 *  1-indexed
 * @param k the number of rows to choose
 * @param A the matrix
 * REQUIRES: k <= A.numRows()
 */
void what_is_missing(unsigned int k, Matrix &A)
{
  int m = A.numRows();
  vector<vector<int>> k_tuples = columns_of_col_sum(m, k);
  for (vector<int> k_tuple : k_tuples) {
    vector<int> rows;
    for (unsigned int i = 0; i < k_tuple.size(); i++) {
      if (k_tuple[i]) 
        rows.push_back(i);
    }
    std::cout << "On rows";
    for (int row : rows) {
      std::cout << " " << row + 1;
    }
    std::cout << ":\n";
    vector<vector<int>> old_rows = A.getRows();
    vector<vector<int>> new_rows;
    for (int row : rows) {
      new_rows.push_back(old_rows[row]);
    }
    vector<vector<int>> cols_of_k = columns_of_K(k);
    vector<vector<int>> new_cols = rowColSwap(new_rows);
    vector<vector<int>> mat_cols;
    remove_duplicates(new_cols);
    for (vector<int> col : cols_of_k) {
      if (list_contains(col, new_cols))
        continue;
      mat_cols.push_back(col);
    }
    Matrix mat(1, mat_cols);
    std::cout << mat << std::endl;
  }
}

/**
 * Remove duplicates from a vector of vector of ints.
 * @param columns the vector of vector of ints
 */
void remove_duplicates(vector<vector<int>> &columns)
{
  vector<vector<int>> temp;
  for (vector<int> column : columns) {
    if (list_contains(column, temp))
      continue;
    temp.push_back(column);
  }
  columns = temp;
}

/**
 * Perform forb by looking for matrices with start number of columns and going down from there
 * @param m the number of rows
 * @param F the configuration to avoid
 * @param start the number of columns to start from
 */
unsigned int forbDown(unsigned int m, const Matrix &F, unsigned int start)
{
  while (start > 0) {
    vector<vector<int>> all_cols = columns_of_K(m);
    int total_possible = 1 << m;
    vector<vector<int>> set_of_possible_columns = columns_of_col_sum(total_possible, start);
    for (vector<int> columns : set_of_possible_columns) {
      vector<int> curr_cols;
      for (unsigned int i = 0; i < columns.size(); i++) {
        if (columns[i]) 
          curr_cols.push_back(i);
      }
      vector<vector<int>> mat_cols;
      for (int col : curr_cols) {
        mat_cols.push_back(all_cols[col]);
      }
      Matrix testMat(1, mat_cols);
      if (testMat.avoids(F)) {
        return start;
      }
    }
    start--;
  }
  return 0;
}

/**
 * Perform ext by looking for matrices with start number of columns and going down from there
 * @param m the number of rows
 * @param F the configuration to avoid
 * @param start the number of columns to start from
 */
vector<Matrix> extDown(unsigned int m, const Matrix &F, unsigned int start)
{
  while (start > 0) {
    vector<vector<int>> all_cols = columns_of_K(m);
    int total_possible = 1 << m;
    vector<vector<int>> set_of_possible_columns = columns_of_col_sum(total_possible, start);
    vector<Matrix> ret;
    for (vector<int> columns : set_of_possible_columns) {
      vector<int> curr_cols;
      for (unsigned int i = 0; i < columns.size(); i++) {
        if (columns[i]) 
          curr_cols.push_back(i);
      }
      vector<vector<int>> mat_cols;
      for (int col : curr_cols) {
        mat_cols.push_back(all_cols[col]);
      }
      Matrix testMat(1, mat_cols);
      if (testMat.avoids(F)) {
        ret.push_back(testMat);
      }
    }
    if (!ret.empty())
      return ret;
    start--;
  }
  return vector<Matrix>();
}
