#pragma once

#include <Eigen/Sparse>
#include "FunctionBase.hh"

//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================

    class MassSpringProblem2DSparse {
    public:
        using Vec = FunctionBase::Vec;
        using Mat = FunctionBase::Mat;
        using Edge = std::pair<int, int>;
        // sparse matrix type
        using SMat = Eigen::SparseMatrix<double>;
        // triplet
        using T = Eigen::Triplet<double>;

        MassSpringProblem2DSparse(FunctionBase& _spring, const int _n_unknowns) : func_(_spring), n_(_n_unknowns) {
            xe_.resize(func_.n_unknowns());
            ge_.resize(func_.n_unknowns());
            he_.resize(func_.n_unknowns(), func_.n_unknowns());
        }

        ~MassSpringProblem2DSparse() {}

        virtual int n_unknowns() {
            return n_;
        }

        virtual void initial_x(Vec &_x) {
            // use current solution
            for (unsigned int i = 0; i < n_; ++i)
                _x[i] = 0.;
        }

        virtual double eval_f(const Vec &_x) {
            double energy(0);
            //store k and l
            Vec coeff(2);
            
            //------------------------------------------------------//
            //TODO: assemble function values of all the elements (springs)
            for (int i = 0; i < elements_.size(); i++) {
                xe_.setZero();
                xe_[0] = _x[elements_[i].first*2];
                xe_[1] = _x[elements_[i].first*2+1];
                xe_[2] = _x[elements_[i].second*2];
                xe_[3] = _x[elements_[i].second*2+1];

                coeff[0] = ks_[i];
                coeff[1] = ls_[i];

                energy += func_.f(xe_, coeff);
            }
            //------------------------------------------------------//

            return energy;
        }

        virtual void eval_gradient(const Vec &_x, Vec &_g) {
            _g.resize(n_unknowns());
            _g.setZero();
            //store k and l
            Vec coeff(2);

            //------------------------------------------------------//
            //TODO: assemble local gradient vector to the global one
            for (int i = 0; i < elements_.size(); i++) {
                ge_.setZero();
                xe_.setZero();
                xe_[0] = _x[elements_[i].first*2];
                xe_[1] = _x[elements_[i].first*2+1];
                xe_[2] = _x[elements_[i].second*2];
                xe_[3] = _x[elements_[i].second*2+1];

                coeff[0] = ks_[i];
                coeff[1] = ls_[i];

                func_.grad_f(xe_, coeff, ge_);
                _g(elements_[i].first*2) += ge_(0);
                _g(elements_[i].first*2+1) += ge_(1);
                _g(elements_[i].second*2) += ge_(2);
                _g(elements_[i].second*2+1) += ge_(3);
            }
            //------------------------------------------------------//
        }

        //Hessian in dense matrix
        virtual void eval_hessian(const Vec &_x, SMat& _h) {
            //store k and l
            Vec coeff(2);

            //------------------------------------------------------//
            //TODO: assemble local hessian matrix to the global sparse hessian matrix
            //hint: use triplet in eigen to setup the sparse matrix
            std::vector<T> tripletList;
            for (int i = 0; i < elements_.size(); i++) {
                xe_.setZero();
                he_.setZero();
                xe_[0] = _x[elements_[i].first*2];
                xe_[1] = _x[elements_[i].first*2+1];
                xe_[2] = _x[elements_[i].second*2];
                xe_[3] = _x[elements_[i].second*2+1];

                Vec indexes(4);
                indexes[0] = elements_[i].first*2;
                indexes[1] = elements_[i].first*2+1;
                indexes[2] = elements_[i].second*2;
                indexes[3] = elements_[i].second*2+1;

                coeff[0] = ks_[i];
                coeff[1] = ls_[i];

                func_.hess_f(xe_, coeff, he_);

                for(int i = 0; i <= 3; i++) {
                    for(int j = 0; j <= 3; j++) {
                        tripletList.push_back(T(indexes[i],indexes[j],he_(i,j)));
                    }
                }

            }

            _h.setFromTriplets(tripletList.begin(), tripletList.end());
            //------------------------------------------------------//
        }


        void add_spring_element(const int _v_idx0, const int _v_idx1, const double _k = 1., const double _l = 1.)
        {
            if(2*_v_idx0 > (int)n_ || _v_idx0 < 0 || 2*_v_idx1 >= (int)n_ || _v_idx1 < 0 )
                std::cout << "Warning: invalid spring element was added... " << _v_idx0 << " " << _v_idx1 << std::endl;
            else {
                elements_.push_back(Edge(_v_idx0, _v_idx1));
                ks_.push_back(_k);
                ls_.push_back(_l);
            }
        }


    private:
        int n_;
        std::vector<Edge> elements_;
        FunctionBase& func_;
        //vector of constants
        std::vector<double> ks_;
        std::vector<double> ls_;

        // coordinates of each element's nodes
        Vec xe_;
        // gradient for each element
        Vec ge_;
        // hessian for each element
        Mat he_;
    };

//=============================================================================
}


