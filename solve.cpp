#include "solve.h"
#include <cstddef>
#include <cstdio>
#include <assert.h>
#include <alloca.h>
#include<iostream>

int fill_cell(cell_t *cell) {
  /**
   * Returns whether it can write the value to the cell.
   **/

  if (cell->value != UNSET)
    return 0;
  if (cell->domain == 0)
    return 0;
  if (cell->domain & (cell->domain - 1)) // bit hack: zeroes everything from the least sig. bit to the first 1 in the number
    return 0;

  cell->value = cell->domain;
  return 1;
}

void rule1(struct puzzle_t *p, int i, int j) {
  /**
   * If the operation is NONE, then remove all values that aren't target.
   *
   * @param p, pointer to puzzle_t instance
   * @param i, row index
   * @param j, col index
   * @return none
   **/

  //printf("In rule 1\n");
  int idx = i * p->size + j;
  if (p->grid[idx].cage->operation == NONE) {
    p->grid[idx].domain = 1 << (p->grid[idx].cage->target - 1);
  }
}

void rule2(struct puzzle_t *p, int i, int j) {
  /**
   * Removes all values in the cage, row, col from domain.
   *
   * @param p, pointer to puzzle_t instance
   * @param i, row index
   * @param j, col index
   * @return none
   **/

   //printf("In rule 2\n");
  int size = p->size;
  int idx = i * size + j;

  //cage_t *cage = p->grid[idx].cage;

  // remove from the row and col
  for (int k = 0; k < size; ++k) {
    int idx2 = k * size + j;
    int idx3 = i * size + k;
    p->grid[idx].domain &= ~p->grid[idx2].value;
    p->grid[idx].domain &= ~p->grid[idx3].value;
  }

  // remove from the cage
  /*
  for (int k = 0; k < cage->num_cell; k++) {
    int idx2 = cage->positions[k];
    p->grid[idx].domain &= ~p->grid[idx2].value;
  }
    */
}

void rule3(struct puzzle_t *p, int i, int j) {
  /**
   * Performs operation of a cell
   *
   * @param p, pointer to puzzle_t instance
   * @param i, row index
   * @param j, col index
   * @return none
   **/
  int size = p->size;
  int idx = i * size + j;
  cage_t *cage = p->grid[idx].cage;
  Operations op = cage->operation;

  //printf("In rule 3\n");
  if (op == ADD) {
    rule3_add(p, idx, cage);
  } else if (op == SUB) {
    rule3_sub(p, idx, cage);
  } else if (op == MUL) {
    rule3_mul(p, idx, cage);
  } else if (op == DIV) {
    rule3_div(p, idx, cage);
  }
}

void rule3_add(struct puzzle_t *p, int idx, cage_t *cage) {
  // domain of current cell
  //printf("In rule 3 add\n");
  unsigned cur_domain = p->grid[idx].domain;
  unsigned target = cage->target;

  struct recurse_t *others = NULL;
  unsigned *saved_domains = new unsigned[cage->num_cell];

  // create a linked list of cells in the same cage as the current cell,
  // and save their domains in saved_domains
  for (int i = 0; i < cage->num_cell; ++i) {
    unsigned other_cell = cage->positions[i];
    if (other_cell != idx) {
      struct recurse_t *new_recurse = (struct recurse_t *) alloca(sizeof(struct recurse_t));
      new_recurse->position = other_cell;
      new_recurse->next = others;
      new_recurse->cell = &p->grid[other_cell];
      others = new_recurse;
    }
    saved_domains[i] = p->grid[other_cell].domain;
  }

  // for each possible value in this cell's domain,
  // (i) ensure that no cell in the same cage, which is also in the same row or column, has that value in its domain.
  // (ii) call rule3_add_recurse to try out remaining values for other cells
  // (iii) if this doesn't work, remove that possibility from the current cell's domain.
  // (iv) restore the domains of the other cells with older, saved, values
  for (int i = 1; i < (p->size + 1); ++i) {
    if (cur_domain & (1 << (i - 1))) {
      unsigned others_target = target - i;
      for (int j = 0; j < cage->num_cell; ++j) {
        unsigned other_cell = cage->positions[j];
        // Different cell on same row/column
        if (other_cell != idx && (other_cell / p->size == idx / p->size || other_cell % p->size == idx % p->size)) {
          // Temporarily remove the selected value from its domain
          p->grid[other_cell].domain &= ~(1 << (i - 1));
          printf("Removing value %d from cell %d, now has domain %02x\n", i, other_cell, p->grid[other_cell].domain);
        }
      }

	
      printf("Testing ADD: %d (other target = %d) for cell %d\n", i, others_target, idx);
      if (others_target && rule3_add_recurse(p, others, others_target)) {

      } else {
        printf("add: removing %d from %d\n", i, idx);
        cur_domain &= ~(1 << (i - 1));
      }

      for (int i = 0; i < cage->num_cell; ++i) {
        unsigned other_cell = cage->positions[i];
        p->grid[other_cell].domain = saved_domains[i];
      }
    }
  }
  delete[] saved_domains;
  p->grid[idx].domain = cur_domain;
}

bool rule3_add_recurse(struct puzzle_t *p, struct recurse_t *others, unsigned target) {
  //printf("In rule 3 add recurse\n");
  unsigned cur_domain = others->cell->domain;

  // Base case
  printf("In recurse for cell %d\n", others->position);
  if (others->next == NULL) {
    printf("Domain is %02x and target is %d\n", cur_domain, target);
    bool found = (cur_domain & (1 << (target - 1)));
    printf("found is %d\n", found);
    return found;
  }

  unsigned *saved_domains = new unsigned[others->cell->cage->num_cell];
  for (int i = 0; i < others->cell->cage->num_cell; ++i) {
    unsigned other_cell = others->cell->cage->positions[i];
    saved_domains[i] = p->grid[other_cell].domain;
  }

  for (int i = 1; i < (p->size + 1); ++i) {
    printf("Trying out value %d\n", i);
    if (cur_domain & (1 << (i - 1))) {
      unsigned others_target = target - i;
      for (int j = 0; j < others->cell->cage->num_cell; ++j) {
        unsigned idx = others->position;
        unsigned other_cell = others->cell->cage->positions[j];
        // Different cell on same row/column
        printf("Current cell: %d, and testing cell %d\n", idx, other_cell);
        if (other_cell != idx && (other_cell / p->size == idx / p->size || other_cell % p->size == idx % p->size)) {
          printf("Found cell %d\n", other_cell);
          // Temporarily remove the selected value from its domain
          p->grid[other_cell].domain &= ~(1 << (i - 1));
          printf("Removing value %d from cell %d, now has domain %02x\n", i, other_cell, p->grid[other_cell].domain);
        }
      }
      bool finish = others_target > 0 && rule3_add_recurse(p, others->next, others_target);
      for (int i = 0; i < others->cell->cage->num_cell; ++i) {
        unsigned other_cell = others->cell->cage->positions[i];
        p->grid[other_cell].domain = saved_domains[i];
      }
      if (finish) {
        printf("This value works!\n");
        delete[] saved_domains;
        return true;
      }
    }
  }
  delete[] saved_domains;
  return false;
}

void rule3_mul(struct puzzle_t *p, int idx, cage_t *cage) {
  // domain of current cell
  //printf("In rule 3 mul\n");
  unsigned cur_domain = p->grid[idx].domain;
  unsigned target = cage->target;

  struct recurse_t *others = NULL;
  unsigned *saved_domains = new unsigned[cage->num_cell];
  // create a linked list of cells in the same cage as the current cell,
  // and save their domains in saved_domains
  for (int i = 0; i < cage->num_cell; ++i) {
    unsigned other_cell = cage->positions[i];
    if (other_cell != idx) {
      struct recurse_t *new_recurse = (struct recurse_t *) alloca(sizeof(struct recurse_t));
      new_recurse->position = other_cell;
      new_recurse->next = others;
      new_recurse->cell = &p->grid[other_cell];
      others = new_recurse;
    }
    saved_domains[i] = p->grid[other_cell].domain;
  }

  for (int i = 1; i < (p->size + 1); ++i) {
    if (cur_domain & (1 << (i - 1))) {
      unsigned others_target = target / i;
      for (int j = 0; j < cage->num_cell; ++j) {
        unsigned other_cell = cage->positions[j];
        // Different cell on same row/column
        if (other_cell != idx && (other_cell / p->size == idx / p->size || other_cell % p->size == idx % p->size)) {
          // Temporarily remove the selected value from its domain
          p->grid[other_cell].domain &= ~(1 << (i - 1));
        }
      }

      printf("Testing MUL: %d (other target = %d) for cell %d\n", i, others_target, idx);
      if (target % i == 0 && rule3_mul_recurse(p, others, others_target)) {

      } else {
        printf("mul: removing %d from %d\n", i, idx);
        cur_domain &= ~(1 << (i - 1));
      }

      for (int i = 0; i < cage->num_cell; ++i) {
        unsigned other_cell = cage->positions[i];
        p->grid[other_cell].domain = saved_domains[i];
      }
    }
  }
  delete[] saved_domains;
  p->grid[idx].domain = cur_domain;
}

bool rule3_mul_recurse(struct puzzle_t *p, struct recurse_t *others, unsigned target) {
  //printf("In rule 3 mul recurse\n");
  unsigned cur_domain = others->cell->domain;

  // Base case
  printf("In recurse for cell %d\n", others->position);
  if (others->next == NULL) {
    printf("Domain is %02x and target is %d\n", cur_domain, target);
    bool found = (cur_domain & (1 << (target - 1)));
    printf("found is %d\n", found);
    return found;
  }

  unsigned *saved_domains = new unsigned[others->cell->cage->num_cell];
  for (int i = 0; i < others->cell->cage->num_cell; ++i) {
    unsigned other_cell = others->cell->cage->positions[i];
    saved_domains[i] = p->grid[other_cell].domain;
  }

  for (int i = 1; i < (p->size + 1); ++i) {
    printf("Trying out value %d\n", i);
    if (cur_domain & (1 << (i - 1))) {
      unsigned others_target = target - i;
      for (int j = 0; j < others->cell->cage->num_cell; ++j) {
        unsigned idx = others->position;
        unsigned other_cell = others->cell->cage->positions[j];
        // Different cell on same row/column
        printf("Current cell: %d, and testing cell %d\n", idx, other_cell);
        if (other_cell != idx && (other_cell / p->size == idx / p->size || other_cell % p->size == idx % p->size)) {
          printf("Found cell %d\n", other_cell);
          printf("Removing value %d from cell %d\n", i, other_cell);
          // Temporarily remove the selected value from its domain
          p->grid[other_cell].domain &= ~(1 << (i - 1));
        }
      }
      bool finish = target % i && rule3_mul_recurse(p, others->next, others_target);
      for (int i = 0; i < others->cell->cage->num_cell; ++i) {
        unsigned other_cell = others->cell->cage->positions[i];
        p->grid[other_cell].domain = saved_domains[i];
      }
      if (finish) {
        printf("This value works!\n");
        delete[] saved_domains;
        return true;
      }
    }
  }
  delete[] saved_domains;
  return false;
}

void rule3_div(struct puzzle_t *p, int cell_me, struct cage_t *cage) {
  //printf("In rule 3 div\n");
  // my domain
  unsigned me_domain = p->grid[cell_me].domain;
  unsigned target = cage->target;

  // find my partner (division involves exactly two cells in the cage)
  assert(cage->num_cell == 2);
  unsigned cell_other = (cage->positions[0] == cell_me) ? cage->positions[1] : cage->positions[0];
  unsigned other_domain = p->grid[cell_other].domain;

  // go through each bit of the domain and check if other has the complementary bit
  for (int i = 1 ; i < (p->size + 1) ; i ++) {
    if (me_domain & (1 << (i - 1))) {

      // DIV case 1: other / me = target -> other = target * me
      unsigned other_val1 = target * i;
      if (other_val1 <= p->size &&   // make sure the number fits in our range
	  other_domain & (1 << (other_val1 - 1))) {
	continue;  // don't need to remove this value from "me"
      }

      // SUB case 2: me / other = target -> other = me / target
      unsigned other_val2 = i / target;
      if ((other_val2 % target == 0) &&    // make sure it divides evenly
	  (other_domain & (1 << (other_val2 - 1)))) {
	continue;
      }

      // remove this bit because it doesn't work with "other"
      printf("div: removing %d from %d\n", i, cell_me);
      me_domain &= ~(1 << (i - 1));
    }
  }
  
  // put the updated "me_domain" back in the data structure
  p->grid[cell_me].domain = me_domain;
}

void rule3_sub(struct puzzle_t *p, int cell_me, struct cage_t *cage) {
  //printf("In rule 3 sub\n");
  // my domain
  unsigned me_domain = p->grid[cell_me].domain;
  unsigned target = cage->target;

  // find my partner (subtract involves exactly two cells in the cage)
  assert(cage->num_cell == 2);
  unsigned cell_other = (cage->positions[0] == cell_me) ? cage->positions[1] : cage->positions[0];
  unsigned other_domain = p->grid[cell_other].domain;

  // go through each bit of the domain and check if other has the complementary bit
  for (int i = 1 ; i < (p->size + 1); i ++) {
    if (me_domain & (1 << (i - 1))) {

      // SUB case 1: other - me = target -> other = target + me
      unsigned other_val1 = target + i;
      if (other_domain & (1 << (other_val1 - 1))) {
	continue;  // don't need to remove this value from "me"
      }

      // SUB case 2: me - other = target -> other = me - target
      unsigned other_val2 = i - target;
      if ((other_val2 > 0) &&
	  (other_domain & (1 << (other_val2 - 1)))) {
	continue;
      }

      // remove this bit because it doesn't work with "other"
      printf("sub: removing %d from %d\n", i, cell_me);
      me_domain &= ~(1 << (i - 1));
    }
  }
  
  // put the updated "me_domain" back in the data structure
  p->grid[cell_me].domain = me_domain;
}

void rule4(struct puzzle_t *p, int i, int j) {
  /**
   * If cell holds value in domain that is not available in any other cell in same column, row, then cell must have that value
   *
   * @param p, pointer to puzzle_t instance
   * @param i, row index
   * @param j, col index
   * @return none
   **/

   //printf("In rule 2\n");
  int size = p->size;
  int idx = i * size + j;
  unsigned int temp_domain = p->grid[idx].domain;
  // check the row
  for (int k = 0; k < size; ++k) {
    if (k == j) continue;
    int idx3 = i * size + k;
    temp_domain &= ~p->grid[idx3].domain;
  }

  if (temp_domain > 0) {
    p->grid[idx].domain = temp_domain;
    return;
  }

  // need to reset temp_domain
  temp_domain = p->grid[idx].domain;
  //check the column
  for (int k = 0; k < size; ++k) {
    if (k == i) continue;
    int idx2 = k * size + j;
    temp_domain &= ~p->grid[idx2].domain;
  }

  if (temp_domain > 0) {
    p->grid[idx].domain = temp_domain;
    return;
  }
}

bool is_done(struct puzzle_t *p) {
  // compute the number of cells on the board
  int size_squared = p->size * p->size;

  // walk through the cells, checking if they all have a single value
  for (int i = 0 ; i < size_squared ; i ++) {
    if (p->grid[i].value == UNSET) {
      // found a cell still with multiple values; we're not done
      return false;
    }
  }
  return true;
}

void remove_invalid(struct puzzle_t *p, struct rule_ll *rule) {
  /**
   * accesses a bunch of function pointers for the rules
   * say, this is a linked list
   * (ask the student to instantiate this themselves)
   **/
  if (rule == NULL)
    return;

  int size = p->size;
  for (int i = 1; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      int idx = i * size + j;
      if (p->grid[idx].value)
        continue;
      rule->rule(p, i, j);
      fill_cell(p->grid + (i * p->size + j));
    }
  }

  remove_invalid(p, rule->next);
}

bool rule4_v2(struct puzzle_t *p, CageOptions *potential_solutions, struct cage_t *cages) {      // assuming potential_solutions is filled in during rule3     
  //unsigned total_cells = p->size * p->size;
  unsigned hypothesis[PUZZLE_SIZE*PUZZLE_SIZE] = {0};
  int found = backtrack_cages(p, 0, hypothesis, potential_solutions, cages);
  if(found) {
	std::cerr<<"Found a solution!! Solution is: \n";
	
	for(int i=0; i<PUZZLE_SIZE; i++){
		for(int j=0; j<PUZZLE_SIZE; j++) {
			std::cerr<< hypothesis[i*PUZZLE_SIZE + j]<< " ";	
		}

		std::cerr<<"\n";
	}

	return true;
  } else{
	std::cerr<<"No solutions found :(( Hypothesis grid is: \n";
	for(int i=0; i<PUZZLE_SIZE; i++){
		for(int j=0; j<PUZZLE_SIZE; j++) {
			std::cerr<< hypothesis[i*PUZZLE_SIZE + j]<< " ";	
		}

		std::cerr<<"\n";
	}

	return false;
  }

  return false;
}

int backtrack_cages(struct puzzle_t *p, unsigned cage_num, unsigned hypothesis[PUZZLE_SIZE*PUZZLE_SIZE], CageOptions *potential_solutions, struct cage_t *cages) {
    if (cage_num >= MAX_NUM_CAGES) {                            // all cages assigned
        return 1;                                   
    }
    
    CageOptions this_cage_options = potential_solutions[cage_num];
    unsigned num_options = this_cage_options.actual_num_options;
    
    cage_t *this_cage = &cages[cage_num];   
    unsigned num_cells = this_cage->num_cell;
     
    for (int i = 0; i < num_options; i++) {
        for (int j = 0; j < num_cells; j++) {
            unsigned cell_position = this_cage->positions[j];
            hypothesis[cell_position] = this_cage_options.options[i][j];
        }
        int valid = 0;
        if (check_grid(p, hypothesis)) {                              // check_grid(hypothesis) is a function that checks rows and columns
            valid = backtrack_cages(p, cage_num + 1, hypothesis, potential_solutions, cages);
        }
        if (valid) {
            return 1;
        } else {
            for (int j = 0; j < num_cells; j++) {
                unsigned cell_position = this_cage->positions[j];
                hypothesis[cell_position] = 0;
            }
        }
    }
    return 0;
}

bool check_grid(struct puzzle_t *p, unsigned hypothesis[PUZZLE_SIZE*PUZZLE_SIZE]) {
	int size = p->size;

	int fail = 0;

	// check repetition in rows
	for(int i =0; i<size; i++) {
		for(int j = 0; j<size; j++) {
			if(hypothesis[i*PUZZLE_SIZE+j]!=0) {
				for(int k=j+1; k<size; k++) {
					
					if(hypothesis[i*PUZZLE_SIZE + j] == hypothesis[i*PUZZLE_SIZE + k]){
						fail = 1;
					}
				}
			}
		}
	}
	// check repetition in columns
	for(int j =0; j<size; j++) {
		for(int i = 0; i<size; i++) {
			if(hypothesis[i*PUZZLE_SIZE + j]!=0) {
				for(int k=i+1; k<size; k++) {
					
					if(hypothesis[i*PUZZLE_SIZE + j] == hypothesis[k*PUZZLE_SIZE + j]){
						fail = 1;
					}
				}
			}
		}
	}

	return (!fail);
}

void rule3_v2(struct puzzle_t *puzzle, struct cage_t *cages, unsigned num_cages, struct CageOptions *cageOptionsAll) {
	
	// for each cage in the puzzle
	for (int i = 0; i<num_cages; i++) {
		cage_t *cur_cage = &cages[i];
		unsigned hypothesis[MAX_CAGE_SIZE] = {0}; // array's size must be static
		std::cerr<<"\nprocessing cage with "<<cur_cage->num_cell<<" elements, operator is "<<cur_cage->operation<<" and target is "<<cur_cage->target<<"\n"; 
		recurse_over_cage(puzzle, cur_cage, i, 0, &hypothesis[0], cageOptionsAll);
		//if(!recurse_over_cage(puzzle, cur_cage, i, 0, &hypothesis[0], cageOptionsAll)) {
			//std::cerr<<"something went wrong. no possible sequences found for cage "<<i<<"\n";
		//}
	}
}


void recurse_over_cage(struct puzzle_t *puzzle, struct cage_t *cage, int cage_num, int cage_element, unsigned *hypothesis, struct CageOptions *cageOptionsAll) {
	// if this is the last element we're trying out: 
	// 		for any combo that works, store it in cageOptionsAll	
	if(cage_element >= cage->num_cell) {
		return;
	}

	if((unsigned)cage_element == (cage->num_cell - 1)) {
		std::cerr<<"	trying last element of cage "<<cage_num<<"\n";
		unsigned pos = (cage->positions)[cage_element];
		unsigned domain = (puzzle->grid)[pos].domain;
		for(int i = 1; i<=puzzle->size; i++) {
			// if that value is allowed
			if((domain & (1 << (i-1)))) {
				// enter that value in the cage
				hypothesis[cage_element] = i;
				std::cerr<<"		trying value "<<i<<" for cell "<<cage_element<<"\n";
				if(check_cage(&hypothesis[0], cage->num_cell, cage->target, cage->operation)) {
					// save possibility in cageOptionsAll
					std::cerr<<"found valid cage sequence\n";
					CageOptions *storage = &cageOptionsAll[cage_num];	
					unsigned processed_options = storage->actual_num_options;
					// if processed options is 5, index of the processed option will be 4, hence, next index to process will be 5
					for(int j = 0; j<cage->num_cell; j++) {
						(storage->options)[processed_options][j] = hypothesis[j];
					}
					(storage->actual_num_options) = processed_options+1;
				}
			}
		}
	}else {
		unsigned pos = (cage->positions)[cage_element];
		unsigned domain = (puzzle->grid)[pos].domain;
		for(int i = 1; i<=PUZZLE_SIZE; i++) {
			// if that value is allowed
			if((domain & (1 << (i-1)))) {
				// enter that value in the cage
				hypothesis[cage_element] = i;
				std::cerr<<"	trying value "<<i<<" for cell "<<cage_element<<"\n";
				recurse_over_cage(puzzle, cage, cage_num, cage_element+1, &hypothesis[0], cageOptionsAll);	
			}
		}

	}

}

bool check_cage(unsigned *hypothesis, unsigned num_cell, unsigned target, Operations op) {
	unsigned sum = 0;
	unsigned mul = 1;
	bool temp;
	switch(op){
		case ADD: 
			sum = 0;
				for(int i=0; i<num_cell; i++) {
					sum += hypothesis[i];
				}
			if(sum != target) {
				std::cerr<<"			check cage failed\n";
			}else{
				std::cerr<<"			check cage passed\n";
			}
			return (sum == target);
			break;
		case MUL:
			mul = 1;
				for(int i=0; i<num_cell; i++) {
					mul *= hypothesis[i];
				}
			if(mul != target) {
				std::cerr<<"			check cage failed\n";
			}else{
				
				std::cerr<<"			check cage passed\n";
			}
			return (mul == target);
		case SUB:
			 temp = ( (target == (hypothesis[0] - hypothesis[1])) || (target == (hypothesis[1] - hypothesis[0])) );
			 if(temp) {
				std::cerr<<"			check cage passed\n";
			 }else{
				std::cerr<<"			check cage failed\n";
			 }
			 return temp;
		case DIV:
			temp = ( (target == (hypothesis[0] / hypothesis[1])) || (target == (hypothesis[1] / hypothesis[0])) );
			if(temp) {
				std::cerr<<"			check cage passed\n";
			}else{
				std::cerr<<"			check cage failed\n";
			}

		case NONE:
			if(target == hypothesis[0]) {
				std::cerr<<"			check cage passed\n";
			}else{

				std::cerr<<"			check cage failed\n";
			}
			return (target == hypothesis[0]);
		default:
			return false;
	}
}

void print_hypothesis(unsigned hypothesis[PUZZLE_SIZE*PUZZLE_SIZE]) {
	std::cerr<<"\n";
	for(int i=0; i<PUZZLE_SIZE; i++){
		for(int j=0; j<PUZZLE_SIZE; j++) {
			std::cerr<< hypothesis[i*PUZZLE_SIZE + j]<< " ";	
		}

		std::cerr<<"\n";
	}

}

