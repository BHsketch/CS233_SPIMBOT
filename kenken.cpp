#include "kenken.h"

#include <cassert>
#include <cstring>
#include <random>
#include <cstdio>

std::random_device rd;
std::mt19937 puzzle_rng(rd());

void swap_row(int i, int j, unsigned int **matrix) {
  unsigned int *p = *(matrix + j);
  *(matrix + j) = *(matrix + i);
  *(matrix + i) = p;
}

void swap_col(int i, int j, unsigned int **matrix, int N) {
  for (int k = 0; k < N; k++) {
    std::swap(matrix[k][i], matrix[k][j]);
  }
}

unsigned int **generate_matrix(int N) {
  assert(N >= 2);

  unsigned int **matrix = new unsigned int *[N];
  for (int i = 0; i < N; i++)
    matrix[i] = new unsigned int[N];

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int v = (i + j) % N + 1; // base cyclic Latin square
      // flip in "pair blocks" if both indices are odd
      if ((i % 2 == 1) && (j % 2 == 1)) {
        v = (i - j + N) % N + 1;
      }
      matrix[i][j] = v;
    }
  }

  std::uniform_int_distribution<> rand_int(0, N - 1);
  for (int i = 0; i < N; ++i) {
    int r1 = rand_int(puzzle_rng);
    int r2;
    while ((r2 = rand_int(puzzle_rng)) == r1)
      ;
    int c1 = rand_int(puzzle_rng);
    int c2;
    while ((c2 = rand_int(puzzle_rng)) == c1)
      ;
    // cout << r1 << r2 << c1 << c2 << endl;
    swap_col(c1, c2, matrix, N);
    swap_row(r1, r2, matrix);
  }
  // printmatrix(matrix, N);
  return matrix;
}

bool is_visitable(int row, int col, int N, int **visited) {
  if (row < 0 || col < 0 || row >= N || col >= N) {
    return false;
  }

  if (visited[row][col] == 1) {
    return false;
  }
  return true;
}

Kenken::Kenken(int N) {
  puzzle.size = N;
  puzzle.grid = new cell_t[N * N];

  // First, we want to randomly generate the numbers to fill the slots
  // each row/column should have numbers from 1 to N, no duplicates
  // cout << "generating cell values" << endl;
  unsigned int **matrix = generate_matrix(N);

  // Want to randomly generate cages, so number varies
  int **visited = new int *[N];
  for (int i = 0; i < N; i++) {
    visited[i] = new int[N];
    for (int j = 0; j < N; j++) {
      visited[i][j] = 0;
    }
  }
  std::vector<block_t> blocks;
  // cout << "generating cages" << endl;
  generate_cages(0, 0, cages, N, visited, blocks, matrix);

  // cage_t cage = {'-', 1, 2, new int[2] {0,1}};
  // cages.push_back(cage);
  // cout << "populating cell domain and cages" << endl;
  // cout << int(cages.size()) << " size\n";

  // first calculate domain
  // if it was a 4x4, then domain is 0xf
  // if it was a 8x8, then domain is 0xff
  int numberOfF = N / 4; // eg 9/4 = 2
  int numberOf1 = N % 4; // eg 9 % 4 = 1
  // eg 1 << (8+1), 1000000000
  int temp = 1 << (numberOfF * 4 + numberOf1);
  // 0111111111
  int domain = temp - 1;
  // cout << "domain is: " << domain << endl;

  // finally, populate the puzzle's grid's cell's domains and cages
  int sz = cages.size();
  // cout << sz << "\n";
  for (int i = 0; i < sz; i++) {
    for (int j = 0; j < cages[i].num_cell; j++) {
      /*
      cout << cages[i].num_cell;
      cout << "position is: " << cages[i].positions[j] << endl;
      */
      puzzle.grid[cages[i].positions[j]].domain = domain;
      puzzle.grid[cages[i].positions[j]].cage = &(cages[i]);
    }
  }

  // Freeing memory
  for (int i = 0; i < N; i++) {
    delete[] visited[i];
  }
  delete[] visited;
  for (int i = 0; i < N; i++) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

Kenken::~Kenken() {
  for (size_t i = 0; i < cages.size(); ++i) {
    delete[] cages[i].positions;
  }
  delete[] puzzle.grid;
}

// We randomly generate cages using a DFS method that randomly picks a direction
// to go in
// Randomly decide to make a cage
// If it is a 2 size cage, pick between subtraction and addition
// If it is other sized cage, pick addition.
// cage_t = {operation, target (sum), num_cell (size), positions}
// the positions cage_t expects are assuming it is a 1D array
// eg a 4*4 first row is positions {0,1,2,3}
// second row are positions {4,5,6,7}
// So a cage is possible with positions {0,4}
void Kenken::generate_cages(int row, int col, std::vector<cage_t> &cages, int N,
                            int **visited, std::vector<block_t> &blocks,
                            unsigned int **matrix) {
  /*    for (size_t i = 0; i < N; i++) {
          for (size_t j = 0; j < N; j++) {
              std::cout << visited[i][j] + " ";
          }
              std::cout << "\n";
      }
  */

  if (visited[row][col] == 1) {
    // if we have already visited this, don't go through it again!
    return;
  }
  visited[row][col] = 1;

  // of the four directios we can check next, this keeps track of the actual
  // visitable ones
  std::vector<block_t> visitable;
  // North
  if (is_visitable(row, col - 1, N, visited)) {
    block_t block = {row, col - 1};
    visitable.push_back(block);
  }
  // East
  if (is_visitable(row + 1, col, N, visited)) {
    block_t block = {row + 1, col};
    visitable.push_back(block);
  }
  // South
  if (is_visitable(row, col + 1, N, visited)) {
    block_t block = {row, col + 1};
    visitable.push_back(block);
  }
  // West
  if (is_visitable(row - 1, col, N, visited)) {
    block_t block = {row - 1, col};
    visitable.push_back(block);
  }

  // add current to blocks
  block_t curr_block = {row, col};
  blocks.push_back(curr_block);

  // dead end handling. Must form the cage!
  bool deadEnd = false;
  if (visitable.size() == 0) {
    deadEnd = true;
  }

  static std::uniform_int_distribution<> prob100(0, 99);
  // std::cout << row << ":" << col << std::endl;
  // std::cout << deadEnd << std::endl;
  // std::cout << visitable.size() << std::endl;
  // check the current, see if we should make a cage at this point and reset
  // blocks
  int size = blocks.size();
  if (size == 0) {
    // pass -- should not happen
  } else if (size == 1) {
    int chance = prob100(puzzle_rng); // 0 to 99
    if (deadEnd || chance < 25) {     // 25% chance
      int row = blocks[0].row;
      int col = blocks[0].col;
      // operation 0, target = value,  num_elements, {cell locations (in
      // 1D)}
      // std::cout << "one: " << row*N+col << std::endl;
      cage_t cage = {NONE, matrix[row][col], 1, new unsigned int[1]};
      cage.positions[0] = row * N + col;
      cages.push_back(cage);

      // reset blocks
      blocks.clear();
    }
  } else if (size == 2) {
    int chance = prob100(puzzle_rng);
    if (deadEnd || chance < 50) {
      int add_or_sub = prob100(puzzle_rng);
      if (add_or_sub < 50) {
        sub_cages(blocks, cages, N, matrix);
      } else {
        add_cages(blocks, cages, N, 2, matrix);
      }
      // reset blocks
      blocks.clear();
    }
  } else if (size == 3) {
    int chance = prob100(puzzle_rng);
    if (deadEnd || chance < 75) {
      add_cages(blocks, cages, N, 3, matrix);
      // reset blocks
      blocks.clear();
    }
  } else {
    add_cages(blocks, cages, N, 4, matrix);
    // reset blocks
    blocks.clear();
  }

  // Pick a random direction to travel, repeat until traveled N,S,E,W
  size = visitable.size();
  while (size > 0) {
    // randomly pick from 0 to size-1 (eg if size = 3, pick from 0,1,2)
    std::uniform_int_distribution<> rand_idx(0, size - 1);
    int index = rand_idx(puzzle_rng);      // index in range 0 to size-1
    block_t next_block = visitable[index]; // we already did bound checking,
                                           // so don't need to use at()
    // c.begin() refers to first element in vector. c.begin() + 1 is the
    // second element
    visitable.erase(visitable.begin() + index);

    // DFS
    generate_cages(next_block.row, next_block.col, cages, N, visited, blocks,
                   matrix);

    size = visitable.size();
  }
}

void Kenken::sub_cages(std::vector<block_t> &blocks, std::vector<cage_t> &cages,
                       int N, unsigned int **matrix) {
  int row0 = blocks[0].row;
  int col0 = blocks[0].col;
  int value0 = matrix[row0][col0];
  int index0 = row0 * N + col0;

  int row1 = blocks[1].row;
  int col1 = blocks[1].col;
  int value1 = matrix[row1][col1];
  int index1 = row1 * N + col1;

  unsigned int target = abs(value0 - value1);
  // std::cout << "cell_sub: " << index0 << " and " << index1 << std::endl;
  cage_t cage = {SUB, target, 2, new unsigned int[2]}; //  order for sub cages?
  cage.positions[0] = index0;
  cage.positions[1] = index1;
  cages.push_back(cage);
}

void Kenken::add_cages(std::vector<block_t> &blocks, std::vector<cage_t> &cages,
                       int N, unsigned int num_blocks, unsigned int **matrix) {
  unsigned int value = 0;
  int cell_locations[num_blocks];

  for (int i = 0; i < num_blocks; i++) {
    int row = blocks[i].row;
    int col = blocks[i].col;
    value += matrix[row][col];
    cell_locations[i] = (row * N + col);

    // std::cout << "cell: " << cell_locations[i] << std::endl;
  }
  // need to pass a pointer to a memory location for the array part
  unsigned int *cell_locs = new unsigned int[num_blocks];
  std::memcpy(cell_locs, cell_locations, num_blocks * sizeof(int));
  cage_t cage = {ADD, value, num_blocks, cell_locs};
  cages.push_back(cage);
}

namespace {} // namespace
