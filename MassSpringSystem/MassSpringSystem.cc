#include <SpringElement2DWithLength.hh>
#include "MassSpringSystem.hh"
#include <random>

namespace AOPT {

    bool MassSpringSystem::is_convex(const int _spring_element_type, const int _sparsity_type) {
        auto n_vertices = sg_.n_vertices();
        auto n_unknowns = 2*n_vertices;

        Vec points(n_unknowns);
        Mat h(n_unknowns, n_unknowns);
        SMat sh(n_unknowns, n_unknowns);

        int iter = 0;

        if(_spring_element_type == WITHOUT_LENGTH) {
            SpringElement2D se;
            if(_sparsity_type == DENSE) {
                MassSpringProblem2D msp(se, n_unknowns);
                for(size_t i=0; i<sg_.n_edges(); ++i)
                    msp.add_spring_element(sg_.edge(i).first, sg_.edge(i).second, sg_.coefficient(i), sg_.length(i));

                //------------------------------------------------------//
                //Todo: check the convexity of the function in SpringElement2D.hh
                //Hint: randomly set the coordinates of the vertices,
                //see if all the eigenvalues of the hessian matrix (Dense) are >=0

                for(int i = 0; i<int(n_unknowns); i++) {
                    points[i] = rng_.get_random_nd_vector(1)[0];
                }

                msp.eval_hessian(points, h);

                // std::cout << "Eigenvalues:";
                // std::cout << m.eigenvalues();
                // std::cout << std::endl;

                Eigen::VectorXcd eivals = h.eigenvalues();
                std::cout << eivals << std::endl;
                bool convex = true;

                for(int i = 0; i < eivals.rows(); i++) {
                    if(eivals[i].real() < 0.) {
                        convex = false;
                    }
                }

                if(convex) {
                    std::cout << "Problem is Convex" << std::endl;
                } else {
                    std::cout << "Problem is Non-convex" << std::endl;
                }

                //------------------------------------------------------//
            } else if(_sparsity_type == SPARSE) {
                //------------------------------------------------------//
                MassSpringProblem2DSparse msps(se, n_unknowns);
                for(size_t i=0; i<sg_.n_edges(); ++i)
                    msps.add_spring_element(sg_.edge(i).first, sg_.edge(i).second, sg_.coefficient(i), sg_.length(i));

                //check the gradient and hessian
                DerivativeChecker npd;
                npd.check_all(msps);

                //Todo: check the convexity of the function in SpringElement2D.hh
                //Hint: see if all the eigenvalues of the hessian matrix (Sparse) are >=0
                //This is the sparse version and the eigenvalues can be calculated with Spectra library
                for(int i = 0; i<int(n_unknowns); i++) {
                    points[i] = rng_.get_random_nd_vector(1)[0];
                }

                msps.eval_hessian(points, sh);
                //------------------------------------------------------//

            }
        } else if(_spring_element_type == WITH_LENGTH){
            SpringElement2DWithLength se2;
            if(_sparsity_type == DENSE) {
                MassSpringProblem2D msp(se2, n_unknowns);
                for(size_t i=0; i<sg_.n_edges(); ++i)
                    msp.add_spring_element(sg_.edge(i).first, sg_.edge(i).second, sg_.coefficient(i), sg_.length(i));

                //------------------------------------------------------//
                //Todo: check the convexity of the function in SpringElement2D.hh
                //Hint: randomly set the coordinates of the vertices,
                //see if all the eigenvalues of the hessian matrix (Dense) are >=0

                for(int i = 0; i<int(n_unknowns); i++) {
                    points[i] = rng_.get_random_nd_vector(1)[0];
                }

                msp.eval_hessian(points, h);

                // std::cout << "Eigenvalues:";
                // std::cout << m.eigenvalues();
                // std::cout << std::endl;

                Eigen::VectorXcd eivals = h.eigenvalues();
                std::cout << eivals << std::endl;
                bool convex = true;

                for(int i = 0; i < eivals.rows(); i++) {
                    if(eivals[i].real() < 0.) {
                        convex = false;
                    }
                }

                if(convex) {
                    std::cout << "Problem is Convex" << std::endl;
                } else {
                    std::cout << "Problem is Non-convex" << std::endl;
                }
                //------------------------------------------------------//
                //Todo: check the convexity of the function in SpringElement2DWithLength.hh
                //Hint: see if all the eigenvalues of the hessian matrix (Dense) are >=0

                //------------------------------------------------------//
            } else if(_sparsity_type == SPARSE) {
                MassSpringProblem2DSparse msps(se2, n_unknowns);
                for(size_t i=0; i<sg_.n_edges(); ++i)
                    msps.add_spring_element(sg_.edge(i).first, sg_.edge(i).second, sg_.coefficient(i), sg_.length(i));

                //check the gradient and hessian
                DerivativeChecker npd;
                npd.check_all(msps, 1e-5, 1e-2);

                //Todo: check the convexity of the function in SpringElement2DWithLength.hh
                //Hint: see if all the eigenvalues of the hessian matrix (Sparse) are >=0
                //This is the sparse version and the eigenvalues can be calculated with Spectra library

                //------------------------------------------------------//
            }
        }


        return false;
    }

    double MassSpringSystem::initial_system_energy(const int _spring_element_type, const int _sparsity_type) {
        auto n_vertices = sg_.n_vertices();
        auto n_unknowns = 2*n_vertices;

        Vec points(2*n_vertices);
        for(size_t i=0; i<n_vertices; ++i) {
            points[2*i] = sg_.point(i)[0];
            points[2*i+1] = sg_.point(i)[1];
        }


        if(_spring_element_type == WITHOUT_LENGTH) {
            SpringElement2D se;
            MassSpringProblem2D msp(se, n_unknowns);
            for(size_t i=0; i<sg_.n_edges(); ++i)
                msp.add_spring_element(sg_.edge(i).first, sg_.edge(i).second, sg_.coefficient(i), sg_.length(i));

            return msp.eval_f(points);

        } else if(_spring_element_type == WITH_LENGTH){
            SpringElement2DWithLength se2;
            MassSpringProblem2D msp(se2, n_unknowns);
            for(size_t i=0; i<sg_.n_edges(); ++i)
                msp.add_spring_element(sg_.edge(i).first, sg_.edge(i).second, sg_.coefficient(i), sg_.length(i));

            return msp.eval_f(points);
        }

        return -1.;
    }

    void MassSpringSystem::setup_spring_graph() {
        //------------------------------------------------------//
        //Todo: implement the function to setup the spring graph
        //hint: use the functions in SpringGraph.hh
        //first, add vertices
        //then, add edges

        // Adding points
        for(int i = 0; i <= n_grid_x_; i++) {
            for(int j = 0; j <= n_grid_y_; j++) {
                sg_.add_vertex(rng_.get_random_nd_vector(2));
            }
        }

        // Adding edges
        int plt; // Point left top
        int prt; // Point right top
        int plb; // Point left bottom
        int prb; // Point right bottom

        for(int i = 0; i < n_grid_x_; i++) { // For every column
            sg_.add_edge(i,i+1); // Edge to the top of column
            for (int j = 0; j < n_grid_y_; j++) { // For every row
                plt = (j * (n_grid_x_ + 1)) + i;
                prt = plt + 1;
                plb = plt + (n_grid_x_ + 1);
                prb = plb + 1;
                if(i == 0) {
                    sg_.add_edge(plt,plb); // Edge to the left of row
                }
                sg_.add_edge(plt,prb, 1., sqrt(2.));
                sg_.add_edge(prt,plb, 1., sqrt(2.));
                sg_.add_edge(prt,prb);
                sg_.add_edge(plb,prb);
            }
        }

        //------------------------------------------------------//
    }

    int MassSpringSystem::get_grid_index(const int _i, const int _j) const {
        assert(_i<=n_grid_x_ && _j<=n_grid_y_);
        return (n_grid_x_+1)*_j + _i;
    }
}
