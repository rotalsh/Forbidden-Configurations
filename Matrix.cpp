#include "Matrix.h"
#include <iostream>
#include <ostream>
#include <vector>

using std::vector;

unsigned int Matrix::numRows() const {
  return rows;
}

unsigned int Matrix::numCols() const {
  return cols;
}

/**
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

Matrix Matrix::dotWith(Matrix b) {
  return dot(*this, b);
}

void Matrix::printMatrix() {
  std::cout << "_" << std::endl;
  for (unsigned int i = 0; i < rows; i++) {
    std::cout << "| ";
    for (unsigned int j = 0; j < cols; j++) {
      std::cout << elems[i*cols + j] << " ";
    }
    std::cout << "|" << std::endl;
  }
}

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

vector<int> Matrix::getElems() const
{
  return elems;
}

Matrix Matrix::transposeOf() {
  return Matrix(0, getCols());
}

void Matrix::transpose() {
  Matrix A(0, getCols());
  rows = A.rows;
  cols = A.cols;
  elems = A.elems;
}

void Matrix::multiply(int a) {
  for (unsigned int i = 0; i < elems.size(); i++) {
    elems[i] *= a;
  }
}

/**
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

Matrix Matrix::operator*(int a)
{
  vector<int> nums;
  for (unsigned int i = 0; i < elems.size(); i++) {
    nums.push_back(elems[i] * a);
  }
  return Matrix(rows, cols, nums);
}

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
 * @todo don't forget to do this lul
 */
bool Matrix::avoid(const Matrix &b) const
{
  if (rows < b.rows || cols < b.cols) return true;
  if (*this == b) return false;
  return avoid(Node(makeVec(b.cols), vector<int>()), b);
}

/**
 * Private helper function for avoid, checking for different columns first
 */
bool Matrix::avoid(const Node& col_node, const Matrix &b) const
{
  if (col_node.rem.empty()) {
    return avoid(col_node.sel, Node(makeVec(b.rows), vector<int>()), b);
  }
  for (unsigned int i = 0; i < col_node.rem.size(); i++) {
    vector<int> curr_sel = col_node.sel;
    vector<int> curr_rem = col_node.rem;
    curr_sel.push_back(curr_rem[i]);
    curr_rem.erase(curr_rem.begin() + i);    
    bool attempt = avoid(Node(curr_rem, curr_sel), b);
    if (!attempt) return false;
  }
  return true;
}

/**
 * Private helper function for avoid, given a fixed column check different rows
 */
bool Matrix::avoid(const vector<int> &pot_cols, const Node &row_node, const Matrix &b) const
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
    std::cout << Matrix(0, mat) << std::endl;
    return !contains(Matrix(0, mat));
  }
  for (unsigned int i = 0; i < row_node.rem.size(); i++) {
    vector<int> curr_sel = row_node.sel;
    vector<int> curr_rem = row_node.rem;
    curr_sel.push_back(curr_rem[i]);
    curr_rem.erase(curr_rem.begin() + i);    
    bool attempt = avoid(pot_cols, Node(curr_rem, curr_sel), b);
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

bool Matrix::containedIn(const Matrix &b) const
{
  return b.contains(*this);
}

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

Matrix Matrix::operator+=(const Matrix &b)
{
  return *this = *this + b;
}

bool Matrix::operator==(const Matrix &b) const
{
  if (this == &b) return true;
  return rows == b.rows && cols == b.cols && elems == b.elems;
}

/**
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

Matrix operator*(int a, const Matrix &b)
{
  vector<int> nums;
  for (unsigned int i = 0; i < b.getElems().size(); i++) {
    nums.push_back(b.getElems()[i] * a);
  }
  return Matrix(b.numRows(), b.numCols(), nums);
}

std::ostream &operator<<(std::ostream &output, const vector<int> &a)
{
  for (unsigned int i = 0; i < a.size(); i++) {
    output << a[i] << (i == a.size()-1 ? "" : " ");
  }
  return output;
}

std::ostream& operator<<(std::ostream& output, const vector<vector<int>>& a) {
  for (unsigned int i = 0; i < a.size(); i++) {
    output << a[i] << (i == a.size()-1 ? "" : "\n");
  }
  return output;
}

std::ostream &operator<<(std::ostream &output, const Matrix &a)
{
  output << "_" << std::endl;
  unsigned int rows = a.numRows(), cols = a.numCols();
  vector<int> elems = a.getElems();
  for (unsigned int i = 0; i < rows; i++) {
    output << "| ";
    for (unsigned int j = 0; j < cols; j++) {
      output << elems[i*cols + j] << " ";
    }
    output << "|" << (i == rows-1 ? "" : "\n");
  }
  return output;
}

vector<int> makeVec(unsigned int i)
{
  vector<int> ret;
  for (unsigned int k = 0; k < i; k++) {
    ret.push_back(k);
  }
  return ret;
}

/**
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
