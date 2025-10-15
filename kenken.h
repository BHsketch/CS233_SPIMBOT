#ifndef KENKEN_H
#define KENKEN_H

#include <cstddef>
#include <vector>

#define PUZZLE_SIZE 4

struct block_t {
  int row;
  int col;
};

enum Operations : unsigned int {
  NONE = 0,
  ADD = 1,
  SUB = 2,
  MUL = 3,
  DIV = 4,
};

const unsigned UNSET = 0;

struct cage_t {
  // 0, 1, 2, 3, 4 for Operations
  Operations operation;
  // what we want to reach in the arithmetic
  unsigned target;
  // how many cells are in the cage
  unsigned num_cell;
  // an array of indices (row * size + col) for the cage
  unsigned *positions;
};

struct cell_t {
  // one-hot encoding for value (0 for EMPTY)
  unsigned int value;
  // all the possible values left (one-hot)
  unsigned int domain;
  // a pointer to the cage
  struct cage_t *cage;
};

struct puzzle_t {
  unsigned size;
  struct cell_t *grid;
};

// Represents a Puzzle for delivery to the bot. Will destroy itself correctly.
class Kenken {

public:
  // Create a Puzzle of dimensions N x N
  Kenken(int N);
  ~Kenken();

  const std::vector<cage_t> &get_cages() const { return cages; }
  puzzle_t &get_puzzle() { return puzzle; }
  const puzzle_t &get_puzzle() const { return puzzle; }

private:
  void generate_cages(int row, int col, std::vector<cage_t> &cages, int N,
                      int **visited, std::vector<block_t> &blocks,
                      unsigned int **matrix);
  void sub_cages(std::vector<block_t> &blocks, std::vector<cage_t> &cages,
                 int N, unsigned int **matrix);
  void add_cages(std::vector<block_t> &blocks, std::vector<cage_t> &cages,
                 int N, unsigned int num_blocks, unsigned int **matrix);

private:
  puzzle_t puzzle;
  std::vector<cage_t> cages;
};

#endif // #ifndef KENKEN_H
