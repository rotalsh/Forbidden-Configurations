#include <iostream>
#include "Matrix.h"
#include <vector>

using std::vector;

void part1() {
  vector<int> kek = {0, 1, 1, 0};
  Matrix A(2, 2, kek);
  A.printMatrix();

  vector<int> lul = {0, 1, 1};
  Matrix B(3, 1, lul);
  Matrix C(1, 3, lul);
  B.printMatrix();
  C.printMatrix();

  Matrix D(1, {{1, 0, 0, 1}, {1, 0, 0, 0}, {1, 0, 1, 0}});
  D.printMatrix();
  vector<vector<int>> ayy = D.getCols();
  Matrix E(1, ayy);
  E.printMatrix();

  vector<vector<int>> byy = D.getRows();
  Matrix F(1, byy);
  F.printMatrix();

  Matrix G = dot(E, F);
  G.printMatrix();

  Matrix H = F.dotWith(E);
  H.printMatrix();

  D.transposeOf().printMatrix();

  D.transpose();
  D.printMatrix();
}

void part2() {
  Matrix A(2, 2, {1, 2, 0, 1});
  Matrix B = A;
  // B.printMatrix();
  // B.multiply(2);
  // B.printMatrix();

  Matrix C(2, 2, {1, 2, 1, 3});
  // (C * B).printMatrix();
  // (2 * C).printMatrix();
  // (C * 2).printMatrix();

  // (A + C).printMatrix();
  // A += C;
  // A.printMatrix();

  std::cout << (A == B) << std::endl;
  Matrix D(3, 2, {1, 1, 0, 1, 0, 1});
  Matrix E(2, 3, {1, 1, 0, 1, 0, 1});
  std::cout << (D == E) << std::endl;
  std::cout << (D == D) << std::endl;
}

void part3() {
  Matrix A(1, {{0, 0, 1}, {0, 1, 1}, {1, 0, 0}, {0, 1, 0}});
  // std::cout << A.isSimple() <<  std::endl;
  // A.printMatrix();
  Matrix B(2, 5, {1, 0, 1, 0, 0, 1, 1, 1, 0, 0});
  Matrix C(1, 2, {1, 0});
  Matrix D(2, 2, {1, 0, 1, 0});
  Matrix E(1, 3, {1, 1, 0});
  Matrix F(2, 2, {1, 0, 0, 1});
  B.printMatrix();
  std::cout << B.contains(C) << std::endl;
  std::cout << B.contains(D) << std::endl;
  std::cout << E.containedIn(B) << std::endl;
  std::cout << B.contains(F) << std::endl;
  std::cout << B.contains(Matrix(3, 2, {1, 0, 0, 1, 1, 0})) << std::endl;
}

void printTesting() {
  Matrix::Node A({0, 1}, {2});
  std::cout << A.rem << std::endl;
  std::cout << A.sel << std::endl;
  Matrix B(2, 3, {1, 2, 3, 4, 5, 6});
  vector<vector<int>> D = B.getRows();
  std::cout << D << std::endl;
  vector<vector<int>> C = rowColSwap(D);
  std::cout << C << std::endl;
  std::cout << B << std::endl;
}

// MORE TESTING NEEDED
void part4() {
  Matrix B(2, 5, {1, 0, 1, 0, 0, 1, 1, 1, 0, 0});
  B.printMatrix();
  Matrix F(2, 2, {1, 0, 0, 1});
  F.printMatrix();
  std::cout << B.avoid(F) << std::endl;
  std::cout << B.contains(F) << std::endl;
}

int main() {
  // part1();
  // part2();
  // part3();
  // printTesting();
  part4();
  return 0;
}