// Copyright (c) 2009-2017 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.



#ifndef __EVALUATOR_THICK_NANOPORE_EFIELD_H__
#define __EVALUATOR_THICK_NANOPORE_EFIELD_H__

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

class EvaluatorThickNanoporeEfield
    {
    public:

        //! type of parameters this external potential accepts
       //!typedef struct param{} param_type;
        typedef Scalar4 param_type;
        typedef Scalar3 field_type;

        //! Constructs the constraint evaluator
        /*! \param X position of particle
            \param box box dimensions
            \param params per-type parameters of external potential
        */
        DEVICE EvaluatorThickNanoporeEfield(Scalar3 X, const BoxDim& box, const param_type& params, const field_type& field)
            : m_pos(X),
              m_box(box),p_center_pos(field), f_params(params)
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
            double pi = M_PI;

            double rx =(double) m_pos.x - p_center_pos.x;
            double ry =(double) m_pos.y - p_center_pos.y;
            double rz =(double) m_pos.z - p_center_pos.z;

            double V0 = (double) f_params.x;
            double a = (double)  f_params.y;
            double c = (double) f_params.z;
            double htp = (double) f_params.w/2.0;

            double sgn_z = rz/abs(rz);

          
            double mu,nu,phi;
            double rho,d1,d2;
            double pref;
            double factor;
            double Vpart;
            double rescale_coeff;

            Scalar3 E;

            E.x = E.y =E.z = Scalar(0);

            rho=sqrt(rx*rx + ry*ry);

            if (abs(rz) > htp){

                
                d1=sqrt( (rho+c)*(rho+c) + (rz-sgn_z*htp)*(rz-sgn_z*htp));
                d2=sqrt( (rho-c)*(rho-c) + (rz-sgn_z*htp)*(rz-sgn_z*htp));


                mu=fabs(acosh( (d1+d2)/(2*c)));
                nu=acos((d1-d2)/(2*c));

                phi=atan2(ry,rx);

                Vpart = V0 / (2.0 + 4*htp/(pi*a));

                rescale_coeff = pi*a*cosh(mu)*sqrt( sinh(mu)*sinh(mu) + sin(nu)*sin(nu) );
                factor=2*Vpart/rescale_coeff;

                pref=(double) 1.0/sqrt( sinh(mu)*sinh(mu) + sin(nu)*sin(nu) );
                
                E.x=Scalar(sgn_z*factor*pref*sinh(mu)*cos(nu)*cos(phi));
                E.y=Scalar(sgn_z*factor*pref*sinh(mu)*cos(nu)*sin(phi));
                E.z=Scalar(factor*pref*cosh(mu)*sin(nu));
            } else if (abs(rz) <= htp) {

                Vpart = V0 / (1.0 + pi*a/(2*htp));
                E.x = E.y = Scalar(0);
                E.z = Scalar(Vpart/(2*htp));

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
            return std::string("thick_np_efield");
            }
        #endif

    protected:
        Scalar3 m_pos;                //!< particle position
        BoxDim m_box;                 //!< box dimensions
        Scalar m_qi;                  //!< particle charge
        Scalar3 p_center_pos;         //!< contains x,y,z position of the pore center
        Scalar4 f_params;             //!< contains V0,a,c,t_p
   };


#endif 
