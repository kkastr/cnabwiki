// Copyright (c) 2009-2017 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.



#ifndef __EVALUATOR_NANOPORE_EFIELD_H__
#define __EVALUATOR_NANOPORE_EFIELD_H__

#ifndef NVCC
#include <string>
#endif

#include <math.h>
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"


// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

//! Class for evaluating the analytic electric field through a nanopore

class EvaluatorNanoporeEfield
    {
    public:

        //! type of parameters this external potential accepts
       //!typedef struct param{} param_type;
        typedef Scalar3 param_type;
        typedef Scalar3 field_type;

        //! Constructs the constraint evaluator
        /*! \param X position of particle
            \param box box dimensions
            \param params per-type parameters of external potential
        */
        DEVICE EvaluatorNanoporeEfield(Scalar3 X, const BoxDim& box, const param_type& params, const field_type& field)
            : m_pos(X),
              m_box(box), p_center_pos(field), f_params(params)
            {
            }

        
        DEVICE static bool needsDiameter() { return false; }
        //! Accept the optional diameter value
        /*! \param di Diameter of particle i
        */
        DEVICE void setDiameter(Scalar di) { }

       
        DEVICE static bool needsCharge() { return true; }

        /*! \param qi Charge of particle i
        */
        DEVICE void setCharge(Scalar qi) { m_qi = qi; }

        //! Declares additional virial cotribututions are needed for the external field
        /*! No contribution
        */
        DEVICE static bool requestFieldVirialTerm() { return true; }

        //! Evaluate the force, energy and virial
        /*! \param F force vector
            \param energy value of the energy
            \param virial array of six scalars for the upper triangular virial tensor
        */
        DEVICE void evalForceEnergyAndVirial(Scalar3& F, Scalar& energy, Scalar* virial)
            {

            // cx,cy,cz are the coordinates of the center of the nanopore
            // a,c are the radii of the ellipse in oblate spheroidal coordinates
            // for a circular pore a = c
            // V0 = effective voltage drop
            Scalar pi = M_PI;

            Scalar rx = m_pos.x - p_center_pos.x;
            Scalar ry = m_pos.y - p_center_pos.y;
            Scalar rz = m_pos.z - p_center_pos.z;



            Scalar V0 = f_params.x;
            Scalar a =  f_params.y;
            Scalar c =  f_params.z;

          
            Scalar mu,nu,phi;
            Scalar rho,d1,d2;
            Scalar pref;
            Scalar factor;
            Scalar rescale_coeff;

            Scalar3 E;

            E.x = E.y =E.z = Scalar(0);

            rho=fast::sqrt(rx*rx + ry*ry);
            d1=fast::sqrt( (rho+c)*(rho+c) + rz*rz);
            d2=fast::sqrt( (rho-c)*(rho-c) + rz*rz);


            mu=fabs(acosh( (d1+d2)/(2*c)));
            nu=fast::acos((d1-d2)/(2*c));

            phi=atan2(ry,rx);


            rescale_coeff = pi*a*cosh(mu)*fast::sqrt( sinh(mu)*sinh(mu) + fast::sin(nu)*fast::sin(nu) );
            
            factor=V0/rescale_coeff;
            
            // this automatically computer the inverse square root
            pref= fast::rsqrt( sinh(mu)*sinh(mu) + fast::sin(nu)*fast::sin(nu) );
            if(rz<0){
                  E.x=Scalar(-factor*pref*sinh(mu)*fast::cos(nu)*fast::cos(phi));
                  E.y=Scalar(-factor*pref*sinh(mu)*fast::cos(nu)*fast::sin(phi));
            }else{
                  E.x=Scalar(factor*pref*sinh(mu)*fast::cos(nu)*fast::cos(phi));
                  E.y=Scalar(factor*pref*sinh(mu)*fast::cos(nu)*fast::sin(phi));
            }
            E.z=Scalar(factor*pref*cosh(mu)*fast::sin(nu));
            
            /*
            Instead of using doubles (which causes slowdown), check if
            the resulting field is NAN in any of the components.
            This is done to avoid the weird issue caused by passing floats close to zero
            through so many trig functions.
            */ 

            if (isnan(E.x)) {
                E.x = Scalar(0.0);
            }
            
            if (isnan(E.y)) {
                E.y = Scalar(0.0);
            }

            if (isnan(E.z)) {
                E.z = Scalar(0.0);
            }



            F.x = m_qi*E.x;
            F.y = m_qi*E.y;
            F.z = m_qi*E.z;
            

      
            }

        #ifndef NVCC
        //! Get the name of this potential
        /*! \returns The potential name. Short and all lowercase.
        */
        static std::string getName()
            {
            return std::string("np_efield");
            }
        #endif

    protected:
        Scalar3 m_pos;                //!< particle position
        BoxDim m_box;                 //!< box dimensions
        Scalar m_qi;                  //!< particle charge
        Scalar3 p_center_pos;         //!< contains position of pore center
        Scalar3 f_params;             //!< contains V0,a,c
        

   };


#endif 
