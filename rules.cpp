#include"solve.h"
#include<iostream>
// initialize these earlier
//const unsigned MAX_NUM_CAGES = 9;
//const unsigned MAX_OPTIONS_PER_CAGE = 324;
//const unsigned MAX_CAGE_SIZE = 5; 

// TODO Write a get_cage(index) function
// TODO In first conditional of backtrack_cages, I think it should be cage_num >= num_cages_for_this_puzzle instead of cage_num >= MAX_NUM_CAGES
//
// Requirement: A pre-populated data structure that contains a list of cages for a given puzzle.

//struct CageOptions {                                            // struct created to track number of actual num options for each cell
    //int options[MAX_OPTIONS_PER_CAGE][MAX_CAGE_SIZE];           // update values in this struct when determining options
    //unsigned actual_num_options;
//};

void rule4(struct puzzle_t *p, CageOptions *potential_solutions, struct cage_t *cages) {      // assuming potential_solutions is filled in during rule3     
  unsigned total_cells = p->size * p->size;
  unsigned hypothesis[total_cells] = {0};
  backtrack_cages(p, 0, hypothesis, potential_solutions, cages);
}

int backtrack_cages(struct puzzle_t *p, unsigned cage_num, unsigned hypothesis[], CageOptions *potential_solutions, struct cage_t *cages) {
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

bool check_grid(struct puzzle_t *p, int hypothesis[]) {
	int size = p->size;

	int fail = 0;

	// check repetition in rows
	for(int i =0; i<size; i++) {
		for(int j = 0; j<size; j++) {
			if(hypothesis[i][j]!=0) {
				for(int k=j+1; k<size; k++) {
					
					if(hypothesis[i][j] == hypothesis[i][k]){
						fail = 1;
					}
				}
			}
		}
	}
	// check repetition in columns
	for(int j =0; j<size; j++) {
		for(int i = 0; i<size; i++) {
			if(hypothesis[i][j]!=0) {
				for(int k=i+1; k<size; k++) {
					
					if(hypothesis[i][j] == hypothesis[k][j]){
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
		if(!recurse_over_cage(puzzle, cur_cage, i, 0, &hypothesis[0], cageOptionsAll)) {
			std::cerr<<"something went wrong. no possible sequences found for cage "<<i<<"\n";
		}
	}
}


bool recurse_over_cage(struct puzzle_t *puzzle, struct cage_t *cage, int cage_num, int cage_element, unsigned *hypothesis, struct CageOptions *cageOptionsAll) {
	// if this is the last element we're trying out: 
	// 		for any combo that works, store it in cageOptionsAll	
	if((unsigned)cage_element == (cage->num_cell - 1)) {
		unsigned pos = (cage->positions)[cage_element];
		unsigned domain = (puzzle->grid)[pos]->domain;
		for(int i = 1; i<=puzzle->size; i++) {
			// if that value is allowed
			if((domain & (1 << (i-1)))) {
				// enter that value in the cage
				hypothesis[pos] = i;
				if(check_cage(&hypothesis[0], num_cell, cage->target, cage->operation)) {
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
	}
}

bool check_cage(unsigned *hypothesis, unsigned num_cell, unsigned target, Operations op) {
	switch(op){
		case ADD: 
			unsigned sum = 0
				for(int i=0; i<num_cell; i++) {
					sum += hypothesis[i];
				}
			return (sum == target);
			break;
		case MUL:
			unsigned mul = 0
				for(int i=0; i<num_cell; i++) {
					mul *= hypothesis[i];
				}
			return (mul == target);
		case SUB:
			return ( (target == (hypothesis[0] - hypothesis[1])) || (target == (hypothesis[1] - hypothesis[0])) );
		case DIV:
			return ( (target == (hypothesis[0] / hypothesis[1])) || (target == (hypothesis[1] / hypothesis[0])) );
		default:
			return (target == hypothesis[0]);
	}
}

/*
 *struct CageOptions {                                            // struct created to track number of actual num options for each cell
    int options[MAX_OPTIONS_PER_CAGE][MAX_CAGE_SIZE];           // update values in this struct when determining options
    unsigned actual_num_options;
};
 *Rule 3: 
1. For each cage in puzzle:
    1. cage_element_idx = 0 // element of the cage we're currently trying out.
    2. int[] hypothesis_cage = {} // new array of cage.num_elements size
    3. call recurse_over_cage(cage_ptr, cage_num, cage_element_idx, hypothesis_cage, cageOptionsAll):// inlined the function algorithm here
        1. if (cage_element == cage.num_elements -1 ) // if this is the last cage element
            1. for each value in domain of cage_element:
                1. works = try_value(hypothesis_cage)
                2. if (!works) return false
                3. else call `save(cage, hypothesis_grid)` // saves this possibility in our cage struct somewhere (or a separate data structure we provide)
        2. else 
            1. for each value in domain of cage_element:
                1. enter value into hypothesis grid
                2. call `recurse_over_cage(cage_ptr, cage_element_idx+1, hypothesis_cage)` // notice we incremented the element index

 *
 * */
