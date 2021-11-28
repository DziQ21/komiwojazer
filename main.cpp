#include "TSP.hpp"

#include <iostream>


int main() {
    cost_matrix_t cm = {{INF, 10, 8,   19, 12},
                      {10, INF, 20,  6,  3},
                      {8,   20, INF, 4,  2},
                      {19,  6,  4, INF,  7},
                      {12,  3,  2,   7, INF}};




//    StageState czesc(cm);
//    std::cout<<czesc.get_matrix();
//
//    auto b=czesc.reduce_cost_matrix();
//    std::cout<<czesc.get_matrix();
//    auto a=czesc.choose_new_vertex();
//    czesc.append_to_path(a.coordinates);
//    czesc.update_cost_matrix(a.coordinates);
//    std::cout<<czesc.get_matrix();
//
//
//    auto c=czesc.reduce_cost_matrix();
//    std::cout<<czesc.get_matrix();
//    auto d=czesc.choose_new_vertex();
//    czesc.append_to_path(d.coordinates);
//    czesc.update_cost_matrix(d.coordinates);
//    std::cout<<czesc.get_matrix();
//    std::cout<<"";
//    std::cout<<std::endl;
//
//    auto e=czesc.reduce_cost_matrix();
//    std::cout<<czesc.get_matrix();
//    auto f=czesc.choose_new_vertex();
//    czesc.append_to_path(f.coordinates);
//    czesc.update_cost_matrix(f.coordinates);
//    std::cout<<czesc.get_matrix();
//    auto lolek =czesc.get_path();
//    std::cout<<"";
//    std::cout<<std::endl;
    /* Rozwiązania:
 * 32 : 3 4 5 2 1
 * 32 : 2 5 4 3 1
 */
//    cost_matrix_t cm {
//            {INF, 12,   3,  45,   6},
//            {78, INF,  90,  21,   3},
//            { 5,  56, INF,  23,  98},
//            {12,   6,   8, INF,  34},
//            { 3,  98,   3,   2, INF}
//    };

    /* Rozwiązanie:
     * 30 : 4 3 2 0 1
    */

//    cost_matrix_t cm {
//            {INF,  3,  4,  2,  7},
//            {3,  INF,  4,  6,  3},
//            {4,  4,  INF,  5,  8},
//            {2,  6,  5,  INF,  6},
//            {7,  3,  8,  6,  INF},
//    };

    /* Rozwiązania:
     * 19 : 4 3 0 2 1
     * 19 : 1 2 0 3 4
    */
    printf("lol");
    tsp_solutions_t solutions = solve_tsp(cm);

    for (auto elem : solutions) {
        std::cout << elem.lower_bound << " : ";
        for (auto pelem : elem.path) {
            std::cout << pelem << " ";
        }
        std::cout << "\n";
    }

    return EXIT_SUCCESS;
}
