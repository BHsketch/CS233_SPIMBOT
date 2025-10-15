#include "kenken.h"
#include <cstddef>
#include <cstdio>

// BHsketch: useful when storing possibilities for each cage in rule3_v2 ---
const unsigned MAX_NUM_CAGES = 9;
const unsigned MAX_OPTIONS_PER_CAGE = 324;
const unsigned MAX_CAGE_SIZE = 5; 

struct CageOptions {                                            // struct created to track number of actual num options for each cell
    int options[MAX_OPTIONS_PER_CAGE][MAX_CAGE_SIZE];           // update values in this struct when determining options
    unsigned actual_num_options;
};
// -------------------------------------------------------------------------

struct rule_ll {
  void (*rule)(struct puzzle_t *, int, int);
  struct rule_ll *next;
};

struct recurse_t {
  unsigned position;
  struct cell_t *cell;
  struct recurse_t *next;
};

struct cell_state {
  unsigned value;
  unsigned domain;
};


int fill_cell(cell_t *cell);
void rule1(struct puzzle_t *p, int i, int j);
void rule2(struct puzzle_t *p, int i, int j);
void rule3(struct puzzle_t *p, int i, int j);
void rule4(struct puzzle_t *p, int i, int j);
void rule3_add(struct puzzle_t *p, int idx, struct cage_t *cage);
void rule3_sub(struct puzzle_t *p, int idx, struct cage_t *cage);
void rule3_mul(struct puzzle_t *p, int idx, struct cage_t *cage);
void rule3_div(struct puzzle_t *p, int idx, struct cage_t *cage);
bool rule3_add_recurse(struct puzzle_t *p, struct recurse_t *others, unsigned target);
bool rule3_mul_recurse(struct puzzle_t *p, struct recurse_t *others, unsigned target);
int solve(struct puzzle_t *p, std::vector<cage_t> &cages);
bool is_done(struct puzzle_t *p);
void remove_invalid(struct puzzle_t *p, struct rule_ll *rule);
void rule3_v2(struct puzzle_t *puzzle, struct cage_t *cages, unsigned num_cages, struct CageOptions *cageOptionsAll);
void recurse_over_cage(struct puzzle_t *puzzle, struct cage_t *cage, int cage_num, int cage_element, unsigned *hypothesis, struct CageOptions *cageOptionsAll);
bool check_cage(unsigned *hypothesis, unsigned num_cell, unsigned target, Operations op);
bool rule4_v2(struct puzzle_t *p, CageOptions *potential_solutions, struct cage_t *cages); 
int backtrack_cages(struct puzzle_t *p, unsigned cage_num, unsigned hypothesis[PUZZLE_SIZE*PUZZLE_SIZE], CageOptions *potential_solutions, struct cage_t *cages);
bool check_grid(struct puzzle_t *p, unsigned hypothesis[PUZZLE_SIZE*PUZZLE_SIZE]);
