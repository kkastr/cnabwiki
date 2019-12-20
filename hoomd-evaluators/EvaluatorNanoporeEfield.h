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
              m_box(box),
              m_field(field), cx(params.x), cy(params.y), cz(params.z)
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

            Scalar rx = m_pos.x - cx;
            Scalar ry = m_pos.y - cy;
            Scalar rz = m_pos.z - cz;

            Scalar V0 =  m_field.x;
            Scalar a = m_field.y;
            Scalar c = m_field.z;

          
            Scalar mu,nu,phi;
            Scalar rho,d1,d2;
            Scalar3 E;
            Scalar pref;
            Scalar factor;

            rho=sqrt(rx*rx + ry*ry);
            d1=sqrt( (rho+c)*(rho+c) + rz*rz);
            d2=sqrt( (rho-c)*(rho-c) + rz*rz);


            mu=fabs(acosh( (d1+d2)/(2*c)));
            nu=acos( (d1-d2)/(2*c));

            phi=atan2(ry,rx);


            Scalar denom = pi*a*cosh(mu)*sqrt( sinh(mu)*sinh(mu) + sin(nu)*sin(nu) );
            factor=V0/denom;

            pref=1.0/sqrt( sinh(mu)*sinh(mu) + sin(nu)*sin(nu) );
            if(rz<0){
                  E.x=-factor*pref*sinh(mu)*cos(nu)*cos(phi);
                  E.y=-factor*pref*sinh(mu)*cos(nu)*sin(phi);
            }else{
                  E.x=factor*pref*sinh(mu)*cos(nu)*cos(phi);
                  E.y=factor*pref*sinh(mu)*cos(nu)*sin(phi);
            }
            E.z=factor*pref*cosh(mu)*sin(nu);


            F = m_qi*E;

            Scalar phi_el = (V0/pi)*atan(sinh(mu));
            energy = -m_qi*phi_el;

            virial[0] = F.x*m_pos.x;
            virial[1] = F.x*m_pos.y;
            virial[2] = F.x*m_pos.z;
            virial[3] = F.y*m_pos.y;
            virial[4] = F.y*m_pos.z;
            virial[5] = F.z*m_pos.z;
      
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
        Scalar3 m_field;              //!< contains V0,a,c
        Scalar cx;                    //!< x position of the pore center
        Scalar cy;                    //!< y position of the pore center
        Scalar cz;                    //!< z position of the pore center
   };


#endif 