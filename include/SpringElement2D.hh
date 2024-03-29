#pragma once

#include "FunctionBase.hh"

//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================

    class SpringElement2D : public FunctionBase {
    public:
        // E_ab(x) = 1/2 * k * ((x[0] - x[2])^2 + (x[1] - x[3])^2)
        // constructor
        SpringElement2D() : FunctionBase() {}

        // number of unknowns
        inline virtual int n_unknowns() override { return 4; }

        // funcion evaluation
        inline virtual double f(const Vec &_x) override {
            return 0.;
        }

        // k and l are stored in _coeffs, i.e._coeffs[0] is k, _coeffs[1] is l
        inline virtual double f(const Vec &_x, const Vec &_coeffs) override {
            //------------------------------------------------------//
            //Todo: implement the function f(x) = 1/2 * k * ((x[0] - x[2])^2 + (x[1] - x[3])^2)
            double result = 1./2. * _coeffs[0] * ((_x[0] - _x[2])*(_x[0] - _x[2]) + (_x[1] - _x[3])*(_x[1] - _x[3]));
            return result;
            //------------------------------------------------------//
            return 0.;
        }

        // gradient evaluation
        inline virtual void grad_f(const Vec &_x, Vec &_g) override {}

        inline virtual void grad_f(const Vec &_x, const Vec &_coeffs, Vec &_g) override {
            //------------------------------------------------------//
            //Todo: implement the gradient and store in _g
            _g[0] = _coeffs[0] * (_x[0]-_x[2]);
            _g[1] = _coeffs[0] * (_x[1]-_x[3]);
            _g[2] = -_coeffs[0] * (_x[0]-_x[2]);
            _g[3] = -_coeffs[0] * (_x[1]-_x[3]);
            //------------------------------------------------------//
        }

        // hessian matrix evaluation
        inline virtual void hess_f(const Vec &_x, Mat &_H) override {}

        inline virtual void hess_f(const Vec &_x, const Vec &_coeffs, Mat &_H) override {
            //------------------------------------------------------//
            //Todo: implement the hessian matrix and store in _H
            for(int i = 0; i <= 3; i++) {
                for(int j = 0; j <= 3; j++) {
                    int sign = ((i < 2 && j<2) ||(i>1 && j>1)) ? 1 : -1;
                    double entry = (i%2 == j%2)? (sign * _coeffs[0]) : 0;
                    _H(i,j) = entry;
                }
            }
            //------------------------------------------------------//

        }

        // initial point for optimization
        inline virtual void initial_x(Vec &_x) override {}

        // callback that is called after optimization finished
        inline virtual void handle_result(const Vec &_x) override {
            std::cout << "solution converged to x = " << std::endl << _x << std::endl;
        }
    };

//=============================================================================
}


