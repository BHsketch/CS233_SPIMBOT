#include "kenken.h"
#include "solve.h"
#include <cmath>
#include <cstdio>

void print_board(const struct puzzle_t *p);
void print_mips_rep(const Kenken &kenken);

int main() {
  Kenken kenken(PUZZLE_SIZE);
  auto cages = kenken.get_cages();
  puzzle_t *puzzle = &kenken.get_puzzle();
  printf("Initial board: \n");
  print_board(puzzle);
  printf("\n Now solving puzzle... \n");
  solve(puzzle, cages); // cages.data returns the underlying contiguous storage from the vector
  print_board(puzzle);
}

int solve(struct puzzle_t *p, std::vector<cage_t> &cages) {
  int size = p->size;
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      rule1(p, i, j);
      fill_cell(p->grid + (i * p->size + j));
    }
  }
  printf("\n After rule 1: \n");
  print_board(p);
  
  // BHsketch ------------------------		

  CageOptions options_for_all_cages[MAX_NUM_CAGES] = {{{{0}},0}};

  unsigned num_cages = (unsigned)(cages.size());
  rule3_v2(p, cages.data(), num_cages, options_for_all_cages);	// pass pointer to the cage options array for rule3 to populate
  bool found = rule4_v2(p, options_for_all_cages, cages.data());	// rule4 uses the populated options_for_all_cages array to backtrack								

  //if(found) {
  //std::cerr<<"Found a solution!! The solution is:\n";

  //} else{
  //std::cerr<<"No solutions found :((\n";
  //}

  // ---------------------------------
  //int change;
  //struct cell_state *buffer = new cell_state[size * size];
  //while(!is_done(p)) {
    //change = 0;
    //for (int i = 0; i < size; ++i) {
      //for (int j = 0; j < size; ++j) {
        //if (p->grid[i * p->size + j].value) continue;
        //rule2(p, i, j);
        //change += fill_cell(p->grid + (i * p->size + j));
      //}
    //}

    //printf("\n After rule 2: \n");
    //print_board(p);
	
	
    //for (int i = 0; i < size; ++i) {
      //for (int j = 0; j < size; ++j) {
        //if (p->grid[i * p->size + j].value) continue;
        //rule3(p, i, j);
        //change += fill_cell(p->grid + (i * p->size + j));
      //}
    //}
    //printf("\n After rule 3: \n");
    //print_board(p);
    //printf("\n");

    //for (int i = 0; i < size; ++i) {
      //for (int j = 0; j < size; ++j) {
        //if (p->grid[i * p->size + j].value) continue;
        //rule4(p, i, j);
        //change += fill_cell(p->grid + (i * p->size + j));
      //}
    //}
    //printf("\n After rule 4: \n");
    //print_board(p);
    //printf("\n");
  
    //// if no change, there are multiple answers. Pick something and propagate changes
    //// Need to save state to backtrack if needed
    //if (change == 0) {
      //printf("No change\n");
      //for (int i = 0; i < size * size; ++i) {
        //if (p->grid[i].value) continue;
        //for (int k = 1; k < (p->size + 1); ++k) {
          //if (p->grid[i].domain & (1 << (k - 1))) {
            //// save state
            //for (int cell = 0; cell < size * size; ++cell) {
              //buffer[cell].domain = p->grid[cell].domain;
              //buffer[cell].value = p->grid[cell].value;
            //}
            //printf("Setting cell %d domain to %02x\n", i, (1 << (k - 1)));
            //p->grid[i].domain = (1 << (k - 1));
            //// try solving with forcing this value
            //if (solve(p)) {
              //// forced value actually works 
              //return 1;
            //} else {
              //// forced value didn't work, restore state, and go to next possible value
              //for (int cell = 0; cell < size * size; ++cell) {
                //p->grid[cell].domain = buffer[cell].domain;
                //p->grid[cell].value = buffer[cell].value;
              //}
            //}
          //}
        //}
      //}
      //// we have tried to force every value, clearly something went wrong
      //delete[] buffer;
      //return 0;
    //}
  //}
  //delete[] buffer;
  return 1;
}

void print_mips_rep(const Kenken &kenken) {
  const auto *p = &kenken.get_puzzle();
  const auto &c = kenken.get_cages();
  printf("puzzle:\t\t.word\t%d\tpuzzle_grid\n", p->size);
  printf("puzzle_grid:\n");

  for (int i = 0; i < p->size * p->size; ++i) {
    printf("cell_%02d:\t.word", i);
    printf("\t%d\t%d", p->grid[i].value, p->grid[i].domain);

    int j;
    for (j = 0; &c[j] != p->grid[i].cage; ++j)
      ;
    printf("\tcage_%02d\n", j);
  }

  for (int i = 0; i < c.size(); ++i) {
    const auto cc = c[i];
    printf("cage_%02d:\t.word", i);
    printf("\t%d\t%d\t%d", cc.operation, cc.target, cc.num_cell);
    printf("\tposition_%02d\n", i);

    printf("position_%02d:\t.word", i);
    for (int j = 0; j < cc.num_cell; ++j) {
      printf("\t%d", cc.positions[j]);
    }
    printf("\n");
  }
}

void _print_row_sep(const struct puzzle_t *p, int i) {
  printf("+");
  for (int j = 0; j < p->size; ++j) {
    if (i < 0 || i + 1 >= p->size ||
        p->grid[i * p->size + j].cage != p->grid[(i + 1) * p->size + j].cage)
      printf("-------+");
    else
      printf("       +");
  }
}

void _print_col_sep(const struct puzzle_t *p, int i, int j) {
  if (j < 0 || j + 1 >= p->size ||
      p->grid[i * p->size + j].cage != p->grid[i * p->size + (j + 1)].cage)
    printf("|");
  else
    printf(" ");
}

void _print_newline(const struct puzzle_t *p) { printf("\n"); }

void _print_cell(const struct puzzle_t *p, int i, int j) {
  cell_t *cell = &p->grid[i * p->size + j];
  cage_t *cage = cell->cage;

  int is_representative = 1;
  for (int k = 0; k < cage->num_cell; ++k) {
    int idx = cage->positions[k];
    int row = idx / p->size;
    int col = idx % p->size;
    if (row < i || col < j)
      is_representative = 0;
  }

  if (is_representative) {
    char ops[] = {' ', '+', '-', 'x', '/'};
    printf("%c%02x  ", ops[cage->operation], cage->target);
  } else {
    printf("     ");
  }
  if (cell->value) {
    printf("%02x", (int)std::log2(cell->value) + 1);
  } else {
    printf("%02X", cell->domain & 0xF);
  }
}

void print_board(const struct puzzle_t *p) {
  /**
   * Prints an ascii representation of the board.
   * Uses a bunch of helper functions to make it look nice.
   **/

  // schematic for ascii representation for grid:
  // each cell is a 1x5 grid (excluding seps)
  // the operation goes on the top left of each cage
  // the number goes on the right of each cell
  // the puzzle is in BLACK, the values are in GRAY
  // +-----+-----+-----+-----+
  // |6×  1|  3  |3÷  3|    1|
  // +-----+-----+-----+-----+
  // |    3|2÷  1|    2|    4|
  // +-----+-----+-----+-----+
  // |1–  2|  1  |6×  3|    2|
  // +-----+-----+-----+-----+
  // |3–  4|  1  |3+  1|    3|
  // +-----+-----+-----+-----+

  _print_row_sep(p, -1);
  for (int i = 0; i < p->size; ++i) {
    _print_newline(p);
    _print_col_sep(p, i, -1);
    for (int j = 0; j < p->size; ++j) {
      _print_cell(p, i, j);
      _print_col_sep(p, i, j);
    }
    _print_newline(p);
    _print_row_sep(p, i);
  }
}
