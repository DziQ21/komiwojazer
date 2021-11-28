#include "TSP.hpp"

#include <algorithm>
#include <stack>
#include <iostream>
#include <optional>

std::ostream& operator<<(std::ostream& os, const CostMatrix& cm) {
    for (std::size_t r = 0; r < cm.size(); ++r) {
        for (std::size_t c = 0; c < cm.size(); ++c) {
            const auto& elem = cm[r][c];
            os << (is_inf(elem) ? "INF" : std::to_string(elem)) << "\t";
        }
        os << "\n";
    }
    os << std::endl;

    return os;
}

/* PART 1 */

/**
 * Create path from unsorted path and last 2x2 cost matrix.
 * @return The vector of consecutive vertex.
 */
path_t StageState::get_path() {
    for(int i=0;i<unsorted_path_.size();i++)
        for(int j=0;j<unsorted_path_.size();j++)
            matrix_[unsorted_path_[i].col][unsorted_path_[j].row]=INF;
    int rest=0;
    for (int i = 0; i < matrix_.size(); ++i)
        for (int j = 0; j < matrix_.size(); ++j)
            if(matrix_[i][j]!=INF)
            {
                rest+=matrix_[i][j];
                unsorted_path_.push_back(vertex_t(i,j));
            }
    lower_bound_+=rest;
    path_t result;
    result.push_back(unsorted_path_[0].row);
    while(result.size()<unsorted_path_.size())
        for(int i=0;i<unsorted_path_.size();i++)
            if(unsorted_path_[i].row==result[result.size()-1])
                result.push_back(unsorted_path_[i].col);
    for(auto &el:result)
        el++;
    return result;
}

/**
 * Get minimum values from each row and returns them.
 * @return Vector of minimum values in row.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_rows() const {
    std::vector<cost_t> result(matrix_.size());
    for(int & i : result)
        i=INF;
    for(int i=0;i<result.size();i++)
        for(int j=0;j<result.size();j++)
            if(result[i]>matrix_[i][j])
                result[i]=matrix_[i][j];
    return result;
}

/**
 * Reduce rows so that in each row at least one zero value is present.
 * @return Sum of values reduced in rows.
 */
cost_t CostMatrix::reduce_rows() {
    auto rows=get_min_values_in_rows();
    for(auto &el:rows)
        if(el==INF)
            el=0;
    for(int i=0;i<rows.size();i++)
        for(int j=0;j<rows.size();j++)
            if(matrix_[i][j]!=INF)
                matrix_[i][j]-=rows[i];
    return std::accumulate(rows.begin(),rows.end(),0);
}

/**
 * Get minimum values from each column and returns them.
 * @return Vector of minimum values in columns.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_cols() const {
    std::vector<cost_t> result(matrix_.size());
    for(int & i : result)
        i=INF;
    for(int i=0;i<result.size();i++)
        for(int j=0;j<result.size();j++)
            if(result[j]>matrix_[i][j])
                result[j]=matrix_[i][j];
    return result;
}

/**
 * Reduces rows so that in each column at least one zero value is present.
 * @return Sum of values reduced in columns.
 */
cost_t CostMatrix::reduce_cols() {
    auto cols=get_min_values_in_cols();
    for(auto &el:cols)
        if(el==INF)
            el=0;
    for(int i=0;i<cols.size();i++)
        for(int j=0;j<cols.size();j++)
            if(matrix_[i][j]!=INF)
                matrix_[i][j]-=cols[j];
    return std::accumulate(cols.begin(),cols.end(),0);
}

/**
 * Get the cost of not visiting the vertex_t (@see: get_new_vertex())
 * @param row
 * @param col
 * @return The sum of minimal values in row and col, excluding the intersection value.
 */
 //dół prawo
cost_t CostMatrix::get_vertex_cost(std::size_t row, std::size_t col) const {
    int a=INF;
    int b=INF;
    for(int i=0;i<matrix_.size();i++)
        if(i!=row and matrix_[i][col]!=INF and matrix_[i][col]<a)
            a=matrix_[i][col];
     for(int i=0;i<matrix_.size();i++)
         if(i!=col and matrix_[row][i]!=INF and matrix_[row][i]<b)
             b=matrix_[row][i];
    return a+b;
}

/* PART 2 */

/**
 * Choose next vertex to visit:
 * - Look for vertex_t (pair row and column) with value 0 in the current cost matrix.
 * - Get the vertex_t cost (calls get_vertex_cost()).
 * - Choose the vertex_t with maximum cost and returns it.
 * @param cm
 * @return The coordinates of the next vertex.
 */
NewVertex StageState::choose_new_vertex() {
    std::vector<NewVertex> zeros;
    for(int i =0;i<matrix_.size();i++)
        for(int j =0;j<matrix_.size();j++)
            if(matrix_.get_matrix()[i][j]==0)
                zeros.push_back(NewVertex(vertex_t(i,j),matrix_.get_vertex_cost(i,j)));
    NewVertex result=zeros[0];
    for(int i =1;i<zeros.size();i++)
        if(zeros[i].cost>result.cost)
            result=zeros[i];
    return result;
}

/**
 * Update the cost matrix with the new vertex.
 * @param new_vertex
 */
void StageState::update_cost_matrix(vertex_t new_vertex) {
    for(int i=0;i<matrix_.size();i++)
    {
        matrix_[new_vertex.row][i]=INF;
        matrix_[i][new_vertex.col]=INF;
    }
    matrix_[new_vertex.col][new_vertex.row]=INF;

}

/**
 * Reduce the cost matrix.
 * @return The sum of reduced values.
 */
cost_t StageState::reduce_cost_matrix() {
    auto a=matrix_.reduce_rows();
    auto b=matrix_.reduce_cols();
    return a+b;
}

/**
 * Given the optimal path, return the optimal cost.
 * @param optimal_path
 * @param m
 * @return Cost of the path.
 */
cost_t get_optimal_cost(const path_t& optimal_path, const cost_matrix_t& m) {
    cost_t cost = 0;

    for (std::size_t idx = 1; idx < optimal_path.size(); ++idx) {
        cost += m[optimal_path[idx - 1] - 1][optimal_path[idx] - 1];
    }

    // Add the cost of returning from the last city to the initial one.
    cost += m[optimal_path[optimal_path.size() - 1] - 1][optimal_path[0] - 1];

    return cost;
}

/**
 * Create the right branch matrix with the chosen vertex forbidden and the new lower bound.
 * @param m
 * @param v
 * @param lb
 * @return New branch.
 */
StageState create_right_branch_matrix(cost_matrix_t m, vertex_t v, cost_t lb) {
    CostMatrix cm(m);
    cm[v.row][v.col] = INF;
    return StageState(cm, {}, lb);
}

/**
 * Retain only optimal ones (from all possible ones).
 * @param solutions
 * @return Vector of optimal solutions.
 */
tsp_solutions_t filter_solutions(tsp_solutions_t solutions) {
    cost_t optimal_cost = INF;
    for (const auto& s : solutions) {
        optimal_cost = (s.lower_bound < optimal_cost) ? s.lower_bound : optimal_cost;
    }

    tsp_solutions_t optimal_solutions;
    std::copy_if(solutions.begin(), solutions.end(),
                 std::back_inserter(optimal_solutions),
                 [&optimal_cost](const tsp_solution_t& s) { return s.lower_bound == optimal_cost; }
    );

    return optimal_solutions;
}

/**
 * Solve the TSP.
 * @param cm The cost matrix.
 * @return A list of optimal solutions.
 */
tsp_solutions_t solve_tsp(const cost_matrix_t& cm) {

    StageState left_branch(cm);

    // The branch & bound tree.
    std::stack<StageState> tree_lifo;

    // The number of levels determines the number of steps before obtaining
    // a 2x2 matrix.
    std::size_t n_levels = cm.size() - 2;

    tree_lifo.push(left_branch);   // Use the first cost matrix as the root.

    cost_t best_lb = INF;
    tsp_solutions_t solutions;

    while (!tree_lifo.empty()) {
        printf("w kolejce %d\n",tree_lifo.size());
        left_branch = tree_lifo.top();
        tree_lifo.pop();

        while (left_branch.get_level() != n_levels && left_branch.get_lower_bound() <= best_lb) {

            // Repeat until a 2x2 matrix is obtained or the lower bound is too high...

            if (left_branch.get_level() == 0) {
                left_branch.reset_lower_bound();
            }

            // 1. Reduce the matrix in rows and columns.
            auto new_cost= left_branch.reduce_cost_matrix(); // @TODO (KROK 1)

            // 2. Update the lower bound and check the break condition.
            left_branch.update_lower_bound(new_cost);
            if (left_branch.get_lower_bound() > best_lb) {
                break;
            }

            // 3. Get new vertex and the cost of not choosing it.

            NewVertex new_vertex = left_branch.choose_new_vertex().coordinates; // @TODO (KROK 2)

            // 4. @TODO Update the path - use append_to_path method.

            left_branch.append_to_path(new_vertex.coordinates);


            // 5. @TODO (KROK 3) Update the cost matrix of the left branch.
            left_branch.update_cost_matrix(new_vertex.coordinates);

            // 6. Update the right branch and push it to the LIFO.

            cost_t new_lower_bound = left_branch.get_lower_bound() + new_vertex.cost;
            tree_lifo.push(create_right_branch_matrix(cm, new_vertex.coordinates,
                                                      new_lower_bound));
            auto test=create_right_branch_matrix(cm, new_vertex.coordinates,
                                            new_lower_bound);
            printf("co jest ?");

        }

        if (left_branch.get_lower_bound() <= best_lb) {

            // If the new solution is at least as good as the previous one,
            // save its lower bound and its path.
            best_lb = left_branch.get_lower_bound();
            path_t new_path = left_branch.get_path();
            solutions.push_back({get_optimal_cost(new_path, cm), new_path});

        }
    }

    return filter_solutions(solutions); // Filter solutions to find only optimal ones.
}
